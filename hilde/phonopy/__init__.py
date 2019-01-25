""" Module containing wrapper functions to work with Phonopy """

import numpy as np
from phonopy.structure.atoms import PhonopyAtoms
from hilde.structure.convert import to_Atoms
from hilde.helpers.fileformats import last_from_yaml
from hilde.helpers.converters import input2dict
from hilde.helpers.attribute_dict import AttributeDict as adict


displacement_id_str = "displacement_id"

defaults = adict(
    {
        # for phono3py compatibility
        "displacement": 0.03,
        "symprec": 1e-5,
        "trigonal": False,
        "is_diagonal": False,
        "q_mesh": [11, 11, 11],
    }
)


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


def get_supercells_with_displacements(phonon):
    """ Create a phonopy object and supercells etc. """

    supercell = to_Atoms(
        phonon.get_supercell(),
        info={
            "supercell": True,
            "supercell_matrix": phonon.get_supercell_matrix().T.flatten().tolist(),
        },
    )

    scells = phonon.get_supercells_with_displacements()

    supercells_with_disps = [to_Atoms(cell) for cell in scells]

    enumerate_displacements(supercells_with_disps)

    return phonon, supercell, supercells_with_disps


def metadata2dict(phonon, calculator):
    """ convert metadata information to plain dict """

    atoms = to_Atoms(phonon.get_primitive())

    prim_data = input2dict(atoms)

    phonon_dict = {
        "version": phonon.get_version(),
        "primitive": prim_data["atoms"],
        "supercell_matrix": phonon.get_supercell_matrix().T.astype(int).tolist(),
        "symprec": float(phonon.get_symmetry().get_symmetry_tolerance()),
        "displacement_dataset": phonon.get_displacement_dataset(),
    }

    try:
        displacements = phonon.get_displacements()
        phonon_dict.update({"displacements": displacements})
    except AttributeError:
        pass

    supercell = to_Atoms(phonon.get_supercell())
    supercell_data = input2dict(supercell, calculator)

    return {str(phonon.__class__.__name__): phonon_dict, **supercell_data}


from .workflow import run_phonopy
from .postprocess import postprocess
