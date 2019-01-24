from hilde.phonopy.wrapper import preprocess as ph_preprocess
from hilde.phono3py.wrapper import preprocess as ph3_preprocess
from pathlib import Path

from hilde.tasks.fireworks.general_py_task import get_func
from hilde.helpers.converters import dict2atoms
from hilde.phonopy import displacement_id_str
from hilde.trajectory import step2file, metadata2file

def run_serial_phonon(atoms, calc, kpt_density, func_path, phonon_settings):
    settings = Settings(settings_file=None)
    settings["atoms"] = dict2atoms(atoms)
    settings["control_kpt"] = {"density": kpt_density}

    if calc["calculator"] != "Aims":
        raise ValueError("The calculator has to be aims")

    settings["control"] = calc["calculator_parameters"]
    if "species_dir" in settings.control:
        settings.control["species_type"] = settings.control.pop("species_dir").split("/")[-1]
    if "use_pimd_wrapper" in phonon_settings:
        settings["socketio"] = {"port": phonon_settings.pop(use_pimd_wrapper)}
    if "aims_command" in settings.control:
        del(settings.control["aims_command"])

    name = func_path.split(".")[-2]
    settings[name] = phonon_settings
    func = get_func(func_path)
    return func(settings=settings)


def preprocess(atoms, calc, kpt_density=None, phonopy_settings=None, phono3py_settings=None):
    '''
    Function to preprocess a phonon calculation
    Args:
        atoms (ASE Atoms object): Structure of the material whose phonons are being calculated
        calc (ASE Calculator): Calculator to calculate the forces
        phonopy_settings (dict): settings for the Phonopy Calculation
        phono3py_settings (dict): settings for the Phono3py Calculation
    Return: tuple
        phonopy preprocess outputs, phono3py preprocess outputs
    '''
    atoms.set_calculator(calc)
    if phonopy_settings:
        phonopy_preprocess = ph_preprocess(atoms, **phonopy_settings)
        if kpt_density is not None:
            update_k_grid(phonopy_preprocess[1], calc, kpt_density)
    else:
        phonopy_preprocess = None
    if phono3py_settings:
        phono3py_preprocess = ph3_preprocess(atoms, **phono3py_settings)
        if kpt_density is not None:
            update_k_grid(phono3py_preprocess[1], calc, kpt_density)
    else:
        phono3py_preprocess = None
    return phonopy_preprocess, phono3py_preprocess

def collect_forces_to_trajectory(
    trajectory,
    calculated_atoms,
    metadata,
):
    Path(trajectory).parents[0].mkdir(exist_ok=True, parents=True)
    for el in metadata["Phonopy"]["displacement_dataset"]["first_atoms"]:
        el["number"] = int(el["number"])
    metadata2file(metadata, trajectory)
    if isinstance(calculated_atoms[0], dict):
        temp_atoms = [dict2atoms(cell) for cell in calculated_atoms]
    else:
        temp_atoms = calculated_atoms.copy()
    calculated_atoms = sorted(
        temp_atoms,
        key=lambda x: x.info[displacement_id_str] if x else len(calculated_atoms) + 1,
    )
    for nn, atoms in enumerate(calculated_atoms):
        if atoms:
            step2file(atoms, atoms.calc, trajectory)