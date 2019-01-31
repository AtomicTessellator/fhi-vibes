"""Generate FWActions after setting Phonon Calculations"""
from ase.symbols import Symbols
from fireworks import FWAction

from phonopy import Phonopy

import numpy as np

from hilde.helpers.converters import calc2dict, atoms2dict
from hilde.helpers.supercell import make_cubic_supercell
from hilde.fireworks.workflow_generator import generate_firework
from hilde.fireworks.workflows.phonon_wf import generate_phonon_fw, generate_phonon_postprocess_fw
from hilde.materials_fp.material_fingerprint import get_phonon_dos_fingerprint_phononpy
from hilde.structure.convert import to_Atoms


def post_bootstrap(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    Generate an FWAction to add either the serial or parallel phonon workflow calculations
    Args:
        atoms (dict): Dictionary representing the ASE Atoms object of theprimitive cell
        calc (dict): Dictionary representing the ASE Calculator for the force calculations
        outputs (dict): The outputs generated from the bootstrap function
        func (str): path to the bootstrap function
        func_fw_out (str): path to this function
        func_kwargs (dict): kwargs used in the bootstrapping
        func_fw_kwargs (dict): kwargs used for the force calculations
        fw_settings (dict): FireWorks settings
    Returns (FWAction): The updated spec and detours to calculate the phonon
                        supercells with displacements
    """
    detours = []
    update_spec = {}
    if "phonopy" in outputs:
        ph_fw_set = fw_settings.copy()
        ph_outputs = outputs["phonopy"]
        ph_settings = func_fw_kwargs["phonopy_settings"].copy()
        update_spec["ph_metadata"] = ph_outputs["metadata"]
        if "spec" in ph_fw_set:
            ph_fw_set["spec"].update(update_spec)
        else:
            ph_fw_set["spec"] = update_spec.copy()
        ph_fw_set["metadata_spec"] = "ph_metadata"
        ph_fw_set["mod_spec_add"] = "ph_forces"
        calc_dict = calc2dict(ph_outputs["calculator"])
        if ph_settings["serial"]:
            update_spec["ph_calculated_atoms"] = [
                atoms2dict(at) for at in ph_outputs["atoms_to_calculate"]
            ]
            update_spec["ph_calculator"] = calc_dict
            ph_fw_set["spec"].update(update_spec)
            ph_fw_set["calc_atoms_spec"] = "ph_calculated_atoms"
            ph_fw_set["calc_spec"] = "ph_calculator"

            detours = add_socket_calc_to_detours(detours, ph_settings, ph_fw_set, "ph")
        else:
            detours = add_single_calc_to_detours(
                detours,
                ph_settings,
                atoms,
                ph_outputs["atoms_to_calculate"],
                calc_dict,
                ph_fw_set,
                "ph",
            )
    if "phono3py" in outputs:
        ph3_fw_set = fw_settings.copy()
        ph3_outputs = outputs["phono3py"]
        ph3_settings = func_fw_kwargs["phono3py_settings"].copy()
        update_spec["ph3_metadata"] = ph3_outputs["metadata"]
        if "spec" in ph3_fw_set:
            ph3_fw_set["spec"].update(update_spec)
        else:
            ph3_fw_set["spec"] = update_spec.copy()
        ph3_fw_set["mod_spec_add"] = "ph3_forces"
        ph3_fw_set["metadata_spec"] = "ph3_metadata"
        calc_dict = calc2dict(ph3_outputs["calculator"])
        if ph3_settings["serial"]:
            update_spec["ph3_calculated_atoms"] = [
                atoms2dict(at) for at in ph3_outputs["atoms_to_calculate"]
            ]
            update_spec["ph3_calculator"] = calc_dict
            ph3_fw_set["spec"].update(update_spec)
            ph3_fw_set["calc_atoms_spec"] = "ph_calculated_atoms"
            ph3_fw_set["calc_spec"] = "ph3_calculator"
            detours = add_socket_calc_to_detours(
                detours, ph3_settings, ph3_fw_set, "ph3"
            )
        else:
            detours = add_single_calc_to_detours(
                detours,
                ph3_settings,
                atoms,
                ph3_outputs["atoms_to_calculate"],
                calc_dict,
                ph3_fw_set,
                "ph3",
            )
    return FWAction(update_spec=update_spec, detours=detours)


def add_socket_calc_to_detours(detours, func_kwargs, fw_settings, prefix):
    """
    Generates a Firework to run a socket calculator and adds it to the detours
    Args:
        detours (list of Fireworks): Current list of detours
        func_kwargs (dict): kwargs needed to do the socket I/O calculation
        fw_settings (dict): FireWorks settings
        prefix (str): ph for phonopy and ph3 for phono3py calculations
    Returns (list of Fireworks): an updated detours list
    """
    calc_kwargs = {}
    calc_keys = ["trajectory", "workdir", "backup_folder", "walltime"]
    for key in calc_keys:
        if key in func_kwargs:
            calc_kwargs[key] = func_kwargs[key]
    fw = generate_firework(
        func="hilde.tasks.fireworks.phonopy_phono3py_functions.wrap_calc_socket",
        func_fw_out="hilde.tasks.fireworks.fw_out.calculate.socket_calc_check",
        func_kwargs=calc_kwargs,
        atoms_calc_from_spec=False,
        inputs=[
            prefix + "_calculated_atoms",
            prefix + "_calculator",
            prefix + "_metadata",
        ],
        fw_settings=fw_settings,
    )
    detours.append(fw)
    return detours


def add_single_calc_to_detours(
    detours, func_fw_kwargs, atoms, atoms_list, calc_dict, fw_settings, prefix
):
    """
    Adds a group of Fireworks to run as single calculations
    Args:
        detours (list of Fireworks): Current list of detours
        func_kwargs (dict): kwargs needed to do the socket I/O calculation
        atoms (dict): Dictionary representing the ASE Atoms object of theprimitive cell
        atosm_list (list of Atoms): List of supercells to perform force calculations on
        calc_dict (dict): Dictionary representing the ASE Calculator for the force calculations
        fw_settings (dict): FireWorks settings
        prefix (str): ph for phonopy and ph3 for phono3py calculations
    Returns (list of Fireworks): an updated detours list
    """
    for i, sc in enumerate(atoms_list):
        if not sc:
            continue
        fw_settings = fw_settings.copy()
        fw_settings["from_db"] = False
        if "kpoint_density_spec" in fw_settings:
            del fw_settings["kpoint_density_spec"]
        sc.info["displacement_id"] = i
        sc_dict = atoms2dict(sc)
        for key, val in calc_dict.items():
            sc_dict[key] = val
        calc_kwargs = {"workdir": func_fw_kwargs["workdir"] + f"/{i:05d}"}
        fw_settings["fw_name"] = (
            prefix + f"forces_{Symbols(atoms['numbers']).get_chemical_formula()}_{i}"
        )
        detours.append(
            generate_firework(
                func="hilde.tasks.calculate.calculate",
                func_fw_out="hilde.tasks.fireworks.fw_out.calculate.mod_spec_add",
                func_kwargs=calc_kwargs,
                atoms=sc_dict,
                calc=calc_dict,
                atoms_calc_from_spec=False,
                fw_settings=fw_settings,
            )
        )
    return detours

def converge_phonons(
    func, func_fw_out, *args, fw_settings=None, **kwargs
):
    ph = kwargs["outputs"]
    prev_dos_fp = None
    if isinstance(ph, Phonopy):
        if "prev_dos_fp" in kwargs:
            prev_dos_fp = kwargs["prev_dos_fp"].copy()
            dos_fp = get_phonon_dos_fingerprint_phononpy(
                ph,
                min_e=np.min(prev_dos_fp[0]),
                max_e=np.max(prev_dos_fp[0]),
                nbins=prev_dos_fp[-1],
            )

        if prev_dos_fp is not None and check_phonon_conv(dos_fp, prev_dos_fp, kwargs["conv_crit"]):
            return FWAction()

        # If Not Converged update phonons
        pc = to_Atoms(ph.get_primitive())
        _, sc_mat = make_cubic_supercell(
            pc,
            len(pc.numbers)*np.linalg.det(ph.get_supercell_matrix())+50
        )

        func_kwargs = {
            "type" : "phonopy"
            "displacement": ph._displacement_dataset['first_atoms'][0]['displacement'][0],
            "supercell_matrix": sc_mat,
        }

        if "spec" in fw_settings:
            fw_settings["spec"]["prev_dos_fp"] = dos_fp
        else:
            fw_settings["spec"] = {"prev_dos_fp": dos_fp}

        if "queueadapter" in fw_settings:
            qadapter = fw_settings["queueadapter"]
        else:
            qadapter = None

        if kwargs["init_wd"].split("/")[-1][9] == "sc_natoms_":
            wd = "/".join(kwargs["init_wd"].split("/")[:-1]) + f"/sc_natoms_{np.linalg.det(sc_mat)*len(pc.numbers)}"
        else:
            wd = kwargs["wd"] + f"/sc_natoms_{np.linalg.det(sc_mat)*len(pc.numbers)}"

        init_fw = generate_phonon_fw(
            pc, wd, fw_settings, qadapter, func_kwargs, update_in_spec=False
        )

        if kwargs["wd"].split("/")[-1][9] == "sc_natoms_":
            wd = "/".join(kwargs["wd"].split("/")[:-1]) + f"/sc_natoms_{np.linalg.det(sc_mat)*len(pc.numbers)}"
        else:
            wd = kwargs["wd"] + f"/sc_natoms_{np.linalg.det(sc_mat)*len(pc.numbers)}"
        kwargs["prev_dos_fp"] = dos_fp
        analysis_fw = generate_phonon_postprocess_fw(
            pc, wd, fw_settings, func_kwargs
        )
        analysis_fw.parents = [init_fw]
        detours = [init_fw, analysis_fw]
        return FWAction(detours=detours, update_spec={"prev_dos_fp": dos_fp})

    from phono3py.phonon3 import Phono3py
    return FWAction()

def check_phonon_conv(dos_fp, prev_dos_fp, conv_crit):
    from hilde.materials_fp.material_fingerprint import scalar_product
    return scalar_product(dos_fp, prev_dos_fp, col=1, normalize=True) >= conv_crit