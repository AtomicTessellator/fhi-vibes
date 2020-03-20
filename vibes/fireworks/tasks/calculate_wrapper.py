"""Wrappers to the vibes calculate functions"""
from pathlib import Path

import numpy as np

from vibes.filenames import filenames
from vibes.helpers.converters import dict2atoms
from vibes.helpers.hash import hash_dict
from vibes.settings import Settings, TaskSettings
from vibes.tasks.calculate import calculate, calculate_socket

T_S_LINE = (
    "          Detailed time accounting                     : "
    " max(cpu_time)    wall_clock(cpu1)\n"
)
E_F_INCONSISTENCY = "  ** Inconsistency of forces<->energy above specified tolerance.\n"


def check_if_failure_ok(lines, walltime):
    """Checks if the FHI-aims calculation finished"""
    line_sum = np.where(lines == T_S_LINE)[0]
    time_used = float(lines[line_sum[0] + 1].split(":")[1].split("s")[1])
    sum_present = line_sum.size > 0

    if walltime and sum_present and time_used / walltime > 0.95:
        return True

    if E_F_INCONSISTENCY in lines:
        return True

    return False


def wrap_calc_socket(
    atoms_dict_to_calculate,
    calculator_dict,
    metadata,
    phonon_times=None,
    mem_use=None,
    trajectory_file=filenames.trajectory,
    workdir=".",
    backup_folder="backups",
    walltime=None,
    fw_settings=None,
    **kwargs,
):
    """Wrapper for the clalculate_socket function

    Parameters
    ----------
    atoms_dict_to_calculate:list of dicts
        A list of dicts representing the cellsto calculate the forces on
    calculator_dict:dict
        A dictionary representation of the ASE Calculator used to calculatethe Forces
    metadata:dict
        metadata for the force trajectory file
    phonon_times:list
        List of all the phonon calculation times
    trajectory_file:str
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
    if calculator_dict["calculator"].lower() == "aims":
        settings = TaskSettings(name=None, settings=Settings(settings_file=None))
        if "species_dir" in calculator_dict["calculator_parameters"]:
            from os import path

            species_type = calculator_dict["calculator_parameters"][
                "species_dir"
            ].split("/")[-1]
            calculator_dict["calculator_parameters"]["species_dir"] = path.join(
                settings.machine.basissetloc, species_type
            )
        calculator_dict["command"] = settings.machine.aims_command
        if walltime:
            calculator_dict["calculator_parameters"]["walltime"] = walltime - 180

    for at_dict in atoms_dict_to_calculate:
        atoms_to_calculate.append(dict2atoms(at_dict, calculator_dict, False))
    calculator = dict2atoms(atoms_dict_to_calculate[0], calculator_dict, False).calc
    if "use_pimd_wrapper" in calculator.parameters:
        if calculator.parameters["use_pimd_wrapper"][0][:5] == "UNIX:":
            atoms_hash = hash_dict({"to_calc": atoms_dict_to_calculate})
            name = atoms_to_calculate[0].get_chemical_formula()
            calculator.parameters["use_pimd_wrapper"][
                0
            ] += f"cm_{name}_{atoms_hash[:15]}"
    try:
        return calculate_socket(
            atoms_to_calculate,
            calculator,
            metadata=metadata,
            trajectory_file=trajectory_file,
            workdir=workdir,
            backup_folder=backup_folder,
            **kwargs,
        )
    except RuntimeError:
        if calculator_dict["calculator"].lower() == "aims":
            path = Path(workdir) / "calculations"
            lines = np.array(open(path / filenames.output.aims).readlines())
            failure_okay = check_if_failure_ok(lines, walltime)
            if not failure_okay:
                raise RuntimeError(
                    "FHI-aims failed to converge, and it is not a walltime issue"
                )
            return True

        raise RuntimeError("The calculation failed")


def wrap_calculate(atoms, calculator, workdir=".", walltime=1800, fw_settings=None):
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
    calculator.parameters["walltime"] = walltime
    calculator.parameters.pop("use_pimd_wrapper", None)
    try:
        return calculate(atoms, calculator, workdir)
    except RuntimeError:
        if calculator.name.lower() == "aims":
            path = Path(workdir)
            lines = np.array(open(path / filenames.output.aims).readlines())
            failure_okay = check_if_failure_ok(lines, walltime)

            if failure_okay:
                return atoms

            raise RuntimeError(
                "FHI-aims failed to converge, and it is not a walltime issue"
            )
        raise RuntimeError("The calculation failed")
