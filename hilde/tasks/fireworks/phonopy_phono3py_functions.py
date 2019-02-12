'''Functions used to wrap around HiLDe Phonopy/Phono3py functions'''
from ase.md.velocitydistribution import PhononHarmonics
from ase import units as u
from pathlib import Path

import numpy as np

from hilde.helpers.converters import dict2atoms
from hilde.phonon_db.row import PhononRow
from hilde.phonopy import displacement_id_str
from hilde.phonopy.workflow import bootstrap
from hilde.settings import Settings, AttributeDict
from hilde.structure.convert import to_Atoms
from hilde.trajectory import step2file, metadata2file
from hilde.tasks.calculate import calculate_socket


def bootstrap_phonon(
    atoms,
    calc,
    kpt_density=None,
    phonopy_settings=None,
    phono3py_settings=None,
    fw_settings=None,
):
    '''
    Creates a Settings object and passes it to the bootstrap function
    Args:
        atoms (ASE Atoms Object): Atoms object of the primitive cell
        calc (ASE Calculator): Calculator for the force calculations
        kpt_density (float): k-point density for the MP-Grid
        phonopy_settings (dict): kwargs for phonopy setup
        phono3py_settings (dict): kwargs for phono3py setup
        fw_settings (dict): FireWork specific settings
    Returns (dict): The output of hilde.phonopy.workflow.bootstrap for phonopy and phono3py
    '''
    settings = Settings(settings_file=None)
    settings.atoms = atoms
    if kpt_density:
        settings["control_kpt"] = AttributeDict({"density": kpt_density})
    kwargs_boot = {}
    if calc.name.lower() != "aims":
        kwargs_boot["calculator"] = calc
    else:
        settings["control"] = calc.parameters.copy()
        if "species_dir" in settings.control:
            sd = settings["control"].pop("species_dir")
            settings["basisset"] = AttributeDict({"type": sd.split("/")[-1]})

        if (phonopy_settings and "use_pimd_wrapper" in phonopy_settings) or (
            phono3py_settings and "use_pimd_wrapper" in phono3py_settings
        ):
            settings["socketio"] = AttributeDict(
                {"port": phonopy_settings.pop("use_pimd_wrapper")}
            )

        if "aims_command" in settings.control:
            del settings.control["aims_command"]

    outputs = {}
    if phonopy_settings:
        settings["phonopy"] = phonopy_settings.copy()
        if "serial" in settings.phonopy:
            del settings.phonopy["serial"]
        outputs["phonopy"] = bootstrap(name="phonopy", settings=settings, **kwargs_boot)
    if phono3py_settings:
        settings["phono3py"] = phono3py_settings.copy()
        if "serial" in settings.phono3py:
            del settings.phono3py["serial"]
        outputs["phono3py"] = bootstrap(
            name="phono3py", settings=settings, **kwargs_boot
        )

    return outputs

def wrap_calc_socket(
    atoms_dict_to_calculate,
    calc_dict,
    metadata,
    trajectory="trajectory.yaml",
    workdir=".",
    backup_folder="backups",
    walltime=1800,
    **kwargs,
):
    '''
    Wrapper for the clalculate_socket function
    Args:
        atoms_dict_to_calculate (list of dicts): A list of dicts representing the cells
                                                 to calculate the forces on
        calc_dict (dict): A dictionary representation of the ASE Calculator used to calculate
                          the Forces
        metadata (dict): metadata for the force trajectory file
        trajectory (str): file name for the trajectory file
        workdir (str): work directory for the force calculations
        backup_folder (str): Directory to store backups
        walltime (int): number of seconds to run the calculation for
    Returns (bool): True if all the calculations completed
    '''
    atoms_to_calculate = []
    if calc_dict["calculator"].lower() == "aims":
        settings = Settings(settings_file=None)
        if "species_dir" in calc_dict["calculator_parameters"]:
            from os import path

            species_type = calc_dict["calculator_parameters"]["species_dir"].split("/")[
                -1
            ]
            calc_dict["calculator_parameters"]["species_dir"] = path.join(
                settings.machine.basissetloc, species_type
            )
        calc_dict["command"] = settings.machine.aims_command

    for at_dict in atoms_dict_to_calculate:
        for key, val in calc_dict.items():
            at_dict[key] = val
        atoms_to_calculate.append(dict2atoms(at_dict))
    calculator = dict2atoms(atoms_dict_to_calculate[0]).calc
    return calculate_socket(
        atoms_to_calculate,
        calculator,
        metadata=metadata,
        trajectory=trajectory,
        workdir=workdir,
        backup_folder=backup_folder,
        walltime=walltime,
        **kwargs,
    )

def collect_to_trajectory(trajectory, calculated_atoms, metadata):
    '''
    Collects forces to a single trajectory file
    Args:
        workdir (str): Path to the work directory
        trajectory (str): file name for the trajectory file
        calculated_atoms (list of ASE Atoms): Results of the force calculations
        metadata (dict): metadata for the phonon calculations
    '''
    traj = Path(trajectory)
    traj.parent.mkdir(exist_ok=True, parents=True)
    if "Phonopy" in metadata:
        for el in metadata["Phonopy"]["displacement_dataset"]["first_atoms"]:
            el["number"] = int(el["number"])

    if "Phono3py" in metadata:
        for el1 in metadata["Phono3py"]["displacement_dataset"]["first_atoms"]:
            el1["number"] = int(el1["number"])
            for el2 in el1["second_atoms"]:
                el2["number"] = int(el2["number"])

    metadata2file(metadata, trajectory)
    if isinstance(calculated_atoms[0], dict):
        temp_atoms = [dict2atoms(cell) for cell in calculated_atoms]
    else:
        temp_atoms = calculated_atoms.copy()
    calculated_atoms = sorted(
        temp_atoms,
        key=lambda x: x.info[displacement_id_str] if x else len(calculated_atoms) + 1,
    )
    for atoms in calculated_atoms:
        if atoms:
            step2file(atoms, atoms.calc, trajectory)

def setup_harmonic_analysis(
    atoms,
    calc,
    phonon_dict,
    temperatures,
    n_samples=1,
    deterministic=False,
    fw_settings=None,
    **kwargs,
):
    '''
    Initializes harmonic analysis functions
    Args:
    Return:
    '''
    phonon_dict["forces_2"] = np.array(phonon_dict["forces_2"])
    ph = PhononRow(dct=phonon_dict).to_phonon()
    n_atoms = ph.get_supercell().get_number_of_atoms()
    force_constants = (
        ph.get_force_constants().swapaxes(1, 2).reshape(2 * (3 * n_atoms,))
    )
    sc = to_Atoms(ph.get_supercell())
    atoms_to_calc = []
    for temp in temperatures:
        atoms_to_calc.append({"temperature": temp, "calc_atoms": []})
        for ii in range(n_samples):
            atoms = sc.copy()
            PhononHarmonics(
                atoms,
                force_constants,
                quantum=False,
                temp=temp * u.kB,
                plus_minus=deterministic,
            )
            atoms_to_calc[-1]["calc_atoms"].append(atoms)
    return atoms_to_calc
