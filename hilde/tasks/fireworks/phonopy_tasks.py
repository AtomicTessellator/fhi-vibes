""" Define FireTasks for electronic structure calculations """
import numpy as np
from fireworks import FWAction

from hilde.helpers.brillouinzone import get_bands_and_labels
from hilde.phonon_db.phonon_db import connect
from hilde.phonon_db.row import phonon_to_dict, PhononRow
from hilde.phonopy import wrapper as ph, displacement_id_str
from hilde.helpers.converters import atoms2dict, dict2atoms, pAtoms
from hilde.tasks.calculate import setup_multiple

module_name = __name__


def initialize_phonopy(smatrix, workdir, atoms_ideal, symprec=1e-5):
    """
    A wrapper function to initialize all phonopy calculations and add the new
    FireWorks as detours in the workflow
    Args:
        atoms_ideal: dict generated from atoms2dict
            The dictionary representation of the atoms or Atoms object of the ndisplaced atoms
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
    print("Initializing phonopy calculations")
    atoms = dict2atoms(atoms_ideal)
    smatrix = np.array(smatrix).reshape(3, 3)
    _, _, supercells_with_disps = ph.preprocess(atoms, smatrix.T, symprec=symprec)
    supercells_with_disps, workdirs = setup_multiple(
        supercells_with_disps, atoms.get_calculator(), workdir, mkdir=False
    )
    for i, cell in enumerate(supercells_with_disps):
        cell.info[displacement_id_str] = i
    atom_dicts = [atoms2dict(cell) for cell in supercells_with_disps]
    return FWAction(update_spec={"atom_dicts": atom_dicts, "workdirs": workdirs})


def calc_phonopy_force_constants(smatrix, atoms_ideal, calc_atoms, symprec=1e-5):
    """
    A wrapper function to calculate 2nd order force constants with phonopy where
    the displacement cells were calculated in its parent FireWork and adds the
    results to a phonon_db
    Args:
        atoms_ideal: dict generated from atoms2dict
            The dictionary representation of the atoms or Atoms object of the undisplayed atoms
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
    print("Calculating harmonic force constants")
    atoms = dict2atoms(atoms_ideal)
    disp_cells = [dict2atoms(ca) for ca in calc_atoms]
    disp_cells = sorted(disp_cells, key=lambda x: x.info[displacement_id_str])
    smatrix = np.array(smatrix).reshape(3, 3)
    phonon, _, _ = ph.preprocess(atoms, smatrix.T, symprec=symprec)
    phonon.set_forces([cell.get_forces() for cell in disp_cells])
    phonon.produce_force_constants()
    return FWAction(update_spec={"phonon_dict": phonon_to_dict(phonon)})


def add_keys(phonon_dict, fw_dict):
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


def calc_phonopy_band_structure(phonon_dict):
    """
    Calculates the phonon dispersion using the phonopy force constants
    Args:
        phonon_dict: dict
            A dictionary representation of the phonopy object
    Returns:
        updated spec for the phonon_dict including the calculated properties
    """
    print("Calculating harmonic band structure")
    phonon = PhononRow(phonon_dict).to_phonon()
    supercell = pAtoms(phonopy_atoms=phonon.get_supercell())
    bands, _ = get_bands_and_labels(supercell)
    phonon.set_band_structure(bands)
    to_fw = add_keys(phonon_dict, phonon_to_dict(phonon, True))
    return FWAction(update_spec={"phonon_dict": to_fw})


def calc_phonopy_dos(mesh, phonon_dict):
    """
    Calculates the phonon density of states using the phonopy force constants
    Args:
        mesh: list of ints (size=3)
            The k-grid mesh sizes
        phonon_dict: dict
            A dictionary representation of the phonopy object
    Returns:
        updated spec for the phonon_dict including the calculated properties
    """
    print("Calculating harmonic DOS")
    phonon = PhononRow(phonon_dict).to_phonon()
    phonon.set_mesh(mesh)
    phonon.set_total_DOS(freq_pitch=0.1, tetrahedron_method=True)
    to_fw = add_keys(phonon_dict, phonon_to_dict(phonon))
    return FWAction(update_spec={"phonon_dict": to_fw})


def calc_phonopy_thermal_prop(mesh, temps, phonon_dict):
    """
    Calculates the thermal properties of a material using the phonopy force constants
    Args:
        mesh: list of ints (size=3)
            The k-grid mesh sizes
        temps: list of floats
            Temperatures to calculate the thermal properties at in Kelvin
        phonon_dict: dict
            A dictionary representation of the phonopy object
    Returns:
        updated spec for the phonon_dict including the calculated properties
    """
    print("Calculating harmonic thermal properties")
    phonon = PhononRow(phonon_dict).to_phonon()
    phonon.set_mesh(mesh)
    phonon.set_thermal_properties(temperatures=temps)
    to_fw = add_keys(phonon_dict, phonon_to_dict(phonon))
    return FWAction(update_spec={"phonon_dict": to_fw})


initialize_phonopy.name = f"{module_name}.{initialize_phonopy.__name__}"
calc_phonopy_force_constants.name = f"{module_name}.{calc_phonopy_force_constants.__name__}"
calc_phonopy_band_structure.name = f"{module_name}.{calc_phonopy_band_structure.__name__}"
calc_phonopy_dos.name = f"{module_name}.{calc_phonopy_dos.__name__}"
calc_phonopy_thermal_prop.name = f"{module_name}.{calc_phonopy_thermal_prop.__name__}"
