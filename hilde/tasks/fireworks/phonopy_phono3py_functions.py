from hilde.phonopy.wrapper import preprocess as ph_preprocess
from hilde.phono3py.wrapper import preprocess as ph3_preprocess
from pathlib import Path

from hilde.tasks.fireworks.general_py_task import get_func
from hilde.helpers.converters import dict2atoms, dict2results, atoms2dict
from hilde.phonopy import displacement_id_str
from hilde.phonopy.workflow import bootstrap
from hilde.settings import Settings, AttributeDict
from hilde.trajectory import step2file, metadata2file
from hilde.tasks.calculate import calculate_socket

def bootstrap_phonon(atoms, calc, kpt_density, phonopy_settings=None, phono3py_settings=None, fw_settings=None):
    settings = Settings(settings_file=None)
    settings["atoms"] = atoms
    settings["control_kpt"] = AttributeDict({"density": kpt_density})
    kwargs_boot = {}
    if calc.name.lower() != "aims":
        kwargs_boot["calculator"] = calc
    else:
        settings["control"] = calc.parameters.copy()
        if "species_dir" in settings.control:
            sd = settings["control"].pop("species_dir")
            settings["basisset"] = AttributeDict({"type": sd.split("/")[-1]})

        if(
            (phonopy_settings and "use_pimd_wrapper" in phonopy_settings) or
            (phono3py_settings and "use_pimd_wrapper" in phono3py_settings)
        ):
            settings["socketio"] = AttributeDict({"port": phonopy_settings.pop(use_pimd_wrapper)})

        if "aims_command" in settings.control:
            del(settings.control["aims_command"])

    outputs = {}
    if phonopy_settings:
        settings["phonopy"] = phonopy_settings.copy()
        if "serial" in settings.phonopy:
            del(settings.phonopy['serial'])
        outputs["phonopy"] = bootstrap(name="phonopy", settings=settings, **kwargs_boot)
    if phono3py_settings:
        settings["phono3py"] = phono3py_settings.copy()
        if "serial" in settings.phono3py:
            del(settings.phono3py['serial'])
        outputs["phono3py"] = bootstrap(name="phono3py", settings=settings, **kwargs_boot)

    return outputs

def wrap_calc_socket(atoms_dict_to_calculate, calc_dict, metadata, trajectory="trajectory.yaml", workdir=".", backup_folder="backups", walltime=1800, **kwargs):
    atoms_to_calculate = []
    if calc_dict["calculator"].lower() == "aims":
        settings = Settings(settings=None)
        if "species_dir" in calc_dict["calculator_parameters"]:
            from os import path
            species_type = calc_dict["calculator_parameters"]["species_dir"].split("/")[-1]
            calc_dict["calculator_parameters"]["species_dir"] = path.join(settings.machine.basissetloc, species_type)
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
# def preprocess(atoms, calc, kpt_density=None, phonopy_settings=None, phono3py_settings=None):
#     '''
#     Function to preprocess a phonon calculation
#     Args:
#         atoms (ASE Atoms object): Structure of the material whose phonons are being calculated
#         calc (ASE Calculator): Calculator to calculate the forces
#         phonopy_settings (dict): settings for the Phonopy Calculation
#         phono3py_settings (dict): settings for the Phono3py Calculation
#     Return: tuple
#         phonopy preprocess outputs, phono3py preprocess outputs
#     '''
#     atoms.set_calculator(calc)
#     if phonopy_settings:
#         phonopy_preprocess = ph_preprocess(atoms, **phonopy_settings)
#         if kpt_density is not None:
#             update_k_grid(phonopy_preprocess[1], calc, kpt_density)
#     else:
#         phonopy_preprocess = None
#     if phono3py_settings:
#         phono3py_preprocess = ph3_preprocess(atoms, **phono3py_settings)
#         if kpt_density is not None:
#             update_k_grid(phono3py_preprocess[1], calc, kpt_density)
#     else:
#         phono3py_preprocess = None
#     return phonopy_preprocess, phono3py_preprocess

def collect_to_trajectory(
    workdir,
    trajectory,
    calculated_atoms,
    metadata,
):
    trajectory = Path(workdir) / trajectory
    Path(workdir).mkdir(exist_ok=True, parents=True)
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
    for nn, atoms in enumerate(calculated_atoms):
        if atoms:
            step2file(atoms, atoms.calc, trajectory)