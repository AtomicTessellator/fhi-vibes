""" tools for storing phonopy supercells

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import numpy as np
from . import input2dict, results2dict
from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml
from hilde.structure.convert import to_Atoms
from hilde.phonopy import displacement_id_str
from hilde.helpers.hash import hash_atoms


def step2file(atoms, calc, displacement_id, file="phonopy_trajectory.yaml"):
    """ Save the current state of MD to file """

    step = step2dict(atoms, calc, displacement_id)

    to_yaml(step, file)


def metadata2file(atoms, calc, obj, file="phonopy_trajectory.yaml"):
    """ save opt metadata to file """

    metadata = metadata2dict(atoms, calc, obj)

    to_yaml(metadata, file, mode="w")


def step2dict(atoms, calc, displacement_id):
    """ extract information from opt step and convet to plain dict """

    atoms_hash = hash_atoms(atoms)
    obj_dict = {displacement_id_str: displacement_id, "hash": atoms_hash}

    return {"info": obj_dict, **results2dict(atoms, calc, append_cell=False)}


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