""" tools for storing phonopy supercells 

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import numpy as np
from . import input2dict, results2dict
from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml
from hilde.structure.convert import to_Atoms


def step2file(atoms, calc, displacement_id, file="phonopy_trajectory.yaml"):
    """ Save the current state of MD to file """

    to_yaml(step2dict(atoms, calc, displacement_id), file)


def metadata2file(atoms, calc, obj, file="phonopy_trajectory.yaml"):
    """ save opt metadata to file """

    metadata = metadata2dict(atoms, calc, obj)

    to_yaml(metadata, file, mode="w")


def step2dict(atoms, calc, displacement_id):
    """ extract information from opt step and convet to plain dict """

    # info from opt algorithm
    obj_dict = {"id": displacement_id}

    return {"phonopy": obj_dict, **results2dict(atoms, calc, append_cell=False)}


def metadata2dict(atoms, calc, obj):
    """ convert metadata information to plain dict """

    prim_data = input2dict(atoms)

    obj_dict = {
        "supercell_matrix": obj.get_supercell_matrix().astype(int).tolist(),
        "symprec": float(obj._symprec),
        "displacements": obj._displacements,
        "displacement_dataset": obj.get_displacement_dataset(),
        "primitive": prim_data["atoms"],
    }

    supercell = to_Atoms(obj.get_supercell())
    supercell_data = input2dict(supercell, calc)

    return {"phonopy": obj_dict, **supercell_data}
