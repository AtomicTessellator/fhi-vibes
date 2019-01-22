""" Module containing wrapper functions to work with Phonopy """

import numpy as np
from phonopy.structure.atoms import PhonopyAtoms
from hilde.structure.convert import to_Atoms
from hilde.helpers.fileformats import last_from_yaml
from hilde.trajectory import input2dict


displacement_id_str = "displacement_id"


def last_calculation_id(trajectory):
    """ return the id of the last computed supercell """
    disp_id = -1

    try:
        dct = last_from_yaml(trajectory)
        disp_id = dct["info"][displacement_id_str]
    except (FileNotFoundError, KeyError):
        pass

    return disp_id


def to_phonopy_atoms(atoms):
    phonopy_atoms = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.get_cell(),
        masses=atoms.get_masses(),
        positions=atoms.get_positions(wrap=True),
    )
    return phonopy_atoms


def enumerate_displacements(cells, info_str=displacement_id_str):
    """ Assign a displacemt id to every atoms obect in cells.

    Args:
        cells (list): atoms objects created by, e.g., phonopy
        info_str (str): how to name the child

    Returns:
        list: cells with id attached to atoms.info (inplace)

    """
    for nn, scell in enumerate(cells):
        if scell is None:
            continue
        scell.info[info_str] = nn


def metadata2dict(atoms, calc, obj):
    """ convert metadata information to plain dict """

    prim_data = input2dict(atoms)

    obj_dict = {
        "version": obj.get_version(),
        "primitive": prim_data["atoms"],
        "supercell_matrix": obj.get_supercell_matrix().astype(int).tolist(),
        "symprec": float(obj._symprec),
        "displacement_dataset": obj.get_displacement_dataset(),
    }

    try:
        displacements = obj.get_displacements()
        obj_dict.update({"displacements": displacements})
    except AttributeError:
        pass

    supercell = to_Atoms(obj.get_supercell())
    supercell_data = input2dict(supercell, calc)

    return {str(obj.__class__.__name__): obj_dict, **supercell_data}
