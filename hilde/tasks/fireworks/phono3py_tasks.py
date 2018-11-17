""" Define FireTasks for electronic structure calculations """
import numpy as np
from fireworks import FWAction
import importlib as il

from hilde.helpers.brillouinzone import get_bands_and_labels
from hilde.phonon_db.phonon_db import connect
from hilde.phonon_db.row import phonon_to_dict, phonon3_to_dict, PhononRow
from hilde.phonopy import phono as ph, displacement_id_str
from hilde.helpers.input_exchange import patoms2dict, dict2patoms, pAtoms
from hilde.tasks.calculate import setup_multiple
from .utility_tasks import add_phonon_to_db

ph3 = il.import_module("hilde.phono3py.phono3")

module_name = __name__

def initialize_phono3py(phono3py_settings, workdir, atoms_ideal, symprec=1e-5):
    """
    A wrapper function to initialize all phonopy calculations and add the new
    FireWorks as detours in the workflow
    Args:
        atoms_ideal: dict generated from patoms2dict
            The dictionary representation of the atoms or pAtoms object of the ndisplaced atoms
        smatrix: list
            The supercell matrix of the calculation
        workdir: str
            The base work directory for the calculation
        symprec: float
            The precision used to determine the system's symmetry
    Returns: FWAction
        A FWAction that will update the spec of the FireWorks to include the proper atom_dicts
        and workdirs
    """
    print("Initializing phono3py calculations")
    atoms = dict2patoms(atoms_ideal)
    print(atoms.calc)
    phono3py_settings["atoms"] = atoms
    phonon3, sc2, sc3, scs2, scs3 = ph3.preprocess(**phono3py_settings)
    print(len(scs2), len(scs3))
    scs, workdirs = setup_multiple(scs2 + scs3, atoms.get_calculator(), workdir, mkdir=False)
    for i, cell in enumerate(scs):
        if cell:
            cell.info[displacement_id_str] = i
    atom_dicts = [patoms2dict(cell) for cell in scs2 + scs3]
    return FWAction(update_spec={"atom_dicts": atom_dicts, "workdirs": workdirs})


def calc_phono3py_force_constants(phono3py_settings, atoms_ideal, calc_atoms, symprec=1e-5):
    """
    A wrapper function to calculate 2nd order force constants with phonopy where
    the displacement cells were calculated in its parent FireWork and adds the
    results to a phonon_db
    Args:
        atoms_ideal: dict generated from patoms2dict
            The dictionary representation of the atoms or pAtoms object of the undisplayed atoms
        smatrix: list
            The supercell matrix of the calculation
        db_path: str
            Path to the phonon_db where these calculations would be included
        calc_atoms: str
            The key in the spec where the displaced atoms with calculated forces are stored
        symprec: float
            The precision used to determine the system's symmetry
    Returns: FWAction
        A FWAction that will store the calculated phonopy objects in the MongoDB for FireWorks
    """
    print("Calculating second and third order force constants")
    atoms = dict2patoms(atoms_ideal)
    phono3py_settings["atoms"] = atoms
    disp_cells = [dict2patoms(ca) for ca in calc_atoms]
    disp_cells = sorted(
        disp_cells, key=lambda x: x.info[displacement_id_str] if x else len(disp_cells) + 1
    )
    phonon3, sc2, sc3, scs2, scs3 = ph3.preprocess(**phono3py_settings)
    fc2_forces = ph3.get_forces(disp_cells[: len(scs2)])
    used_forces = len(scs2)
    fc3_cells = []
    for sc in scs3:
        if sc:
            fc3_cells.append(disp_cells[used_forces])
            used_forces += 1
        else:
            fc3_cells.append(None)
    fc3_forces = ph3.get_forces(fc3_cells)
    phonon3.produce_fc2(fc2_forces)
    phonon3.produce_fc3(fc3_forces)
    return phonon3


def add_keys(phonon3_dict, fw_dict):
    """
    Checks the phonon_dictionary and adds any missing keys
    Args:
        phonon_dict: dict
            The original phonon dictionary
        fw_dict: dict
            The new phonon dictionary
    Returns: dict
        The new phonon dictionary with any keys from the old one that were missing
    """
    for key in phonon_dict:
        if key not in fw_dict:
            fw_dict[key] = phonon_dict[key]
    return fw_dict


def calc_phono3py_kappa(phonon3, temps, mesh=[11, 11, 11]):
    """
    Calculates the phonon dispersion using the phonopy force constants
    Args:
        phonon_dict: dict
            A dictionary representation of the phonopy object
    Returns:
        updated spec for the phonon_dict including the calculated properties
    """
    print("Calculating harmonic band structure")
    phonon3.run_thermal_conductivity(write_kappa=True, temperatures=temps)
    to_fw = add_keys(phonon_dict, phonon_to_dict(phonon, True))
    return FWAction(update_spec={"phonon_dict": to_fw})


def analyze_phono3py(
    db_path,
    phono3py_settings,
    atoms_ideal,
    calc_atoms,
    symprec=1e-5,
    temps=None,
    mesh=None,
    calc_type="calc",
    **kwargs,
):
    phonon3 = calc_phono3py_force_constants(phono3py_settings, atoms_ideal, calc_atoms, symprec)
    if temps is not None:
        if mesh is None:
            calc_phono3py_kappa(phonon3, temps)
        else:
            calc_phono3py_kappa(phonon3, temps, mesh)
    phonon_dict = phonon3_to_dict(phonon3)
    print(phonon_dict["fc_3"].shape)
    # del(phonon_dict["fc_3"])
    db_args = [db_path, atoms_ideal, phonon_dict, calc_type, symprec]
    return add_phonon_to_db(*db_args, **kwargs)


initialize_phono3py.name = f"{module_name}.{initialize_phono3py.__name__}"
calc_phono3py_force_constants.name = f"{module_name}.{calc_phono3py_force_constants.__name__}"
calc_phono3py_kappa.name = f"{module_name}.{calc_phono3py_kappa.__name__}"
analyze_phono3py.name = f"{module_name}.{analyze_phono3py.__name__}"
