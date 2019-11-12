"""Wrappers to the hilde calculate functions"""
import numpy as np

from hilde.helpers.converters import dict2atoms
from hilde.settings import TaskSettings, Settings
from hilde.tasks.calculate import calculate_socket, calculate

TIME_SUMMARY_LINE = "          Detailed time accounting                     :  max(cpu_time)    wall_clock(cpu1)\n"
E_F_INCONSISTENCY = "  ** Inconsistency of forces<->energy above specified tolerance.\n"


def check_if_failure_ok(lines, walltime):
    """Checks if the FHI-aims calculation finished"""
    line_sum = np.where(lines == TIME_SUMMARY_LINE)[0]
    time_used = float(lines[line_sum[0] + 1].split(":")[1].split("s")[1])
    sum_present = line_sum.size > 0

    if walltime and sum_present and time_used / walltime > 0.95:
        return True

    if E_F_INCONSISTENCY in lines:
        return True

    return False


def wrap_calc_socket(
    atoms_dict_to_calculate,
    calc_dict,
    metadata,
    phonon_times=None,
    mem_use=None,
    trajectory="trajectory.son",
    workdir=".",
    backup_folder="backups",
    walltime=None,
    **kwargs,
):
    """Wrapper for the clalculate_socket function

    Parameters
    ----------
    atoms_dict_to_calculate:list of dicts
        A list of dicts representing the cellsto calculate the forces on
    calc_dict:dict
        A dictionary representation of the ASE Calculator used to calculatethe Forces
    metadata:dict
        metadata for the force trajectory file
    phonon_times:list
        List of all the phonon calculation times
    trajectory:str
        file name for the trajectory file
    workdir:str
        work directory for the force calculations
    backup_folder:str
        Directory to store backups
    walltime:int
        number of seconds to run the calculation for

    Returns
    -------
    bool
        True if all the calculations completed

    Raises
    ------
    RuntimeError
        If the calculation fails
    """
    atoms_to_calculate = []
    if calc_dict["calculator"].lower() == "aims":
        settings = TaskSettings(name=None, settings=Settings(settings_file=None))
        if "species_dir" in calc_dict["calculator_parameters"]:
            from os import path

            species_type = calc_dict["calculator_parameters"]["species_dir"].split("/")[
                -1
            ]
            calc_dict["calculator_parameters"]["species_dir"] = path.join(
                settings.machine.basissetloc, species_type
            )
        calc_dict["command"] = settings.machine.aims_command
        if walltime:
            calc_dict["calculator_parameters"]["walltime"] = walltime - 180

    for at_dict in atoms_dict_to_calculate:
        atoms_to_calculate.append(dict2atoms(at_dict, calc_dict, False))
    calculator = dict2atoms(atoms_dict_to_calculate[0], calc_dict, False).calc
    try:
        return calculate_socket(
            atoms_to_calculate,
            calculator,
            metadata=metadata,
            trajectory=trajectory,
            workdir=workdir,
            backup_folder=backup_folder,
            **kwargs,
        )
    except RuntimeError:
        if calc_dict["calculator"].lower() == "aims":
            lines = np.array(open(workdir + "calculations/aims.out").readlines())
            failure_okay = check_if_failure_ok(lines, walltime)
            if not failure_okay:
                raise RuntimeError(
                    "FHI-aims failed to converge, and it is not a walltime issue"
                )
            return True

        raise RuntimeError("The calculation failed")


def wrap_calculate(atoms, calc, workdir=".", walltime=1800):
    """Wrapper for the clalculate_socket function

    Parameters
    ----------
    atoms:Atoms
        Structure.
    calculator:calculator
        Calculator.
    workdir:folder
        Folder to perform calculation in.
    walltime:int
        number of seconds to run the calculation for

    Returns
    -------
    bool
        True if all the calculations completed

    Raises
    ------
    RuntimeError
        If the calculation fails
"""
    calc.parameters["walltime"] = walltime
    try:
        return calculate(atoms, calc, workdir)
    except RuntimeError:
        if calc.name.lower() == "aims":
            lines = np.array(open(workdir + "/aims.out").readlines())
            failure_okay = check_if_failure_ok(lines, walltime)

            if failure_okay:
                return atoms

            raise RuntimeError(
                "FHI-aims failed to converge, and it is not a walltime issue"
            )
        raise RuntimeError("The calculation failed")
