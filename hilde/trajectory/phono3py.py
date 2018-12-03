""" tools for storing phono3py supercells

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import numpy as np
from . import input2dict, results2dict
from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml
from hilde.structure.convert import to_Atoms


def step2file(atoms, calc, displacement_id, file="phono3py_trajectory.yaml"):
    """ Save the current state of MD to file """

    to_yaml(step2dict(atoms, calc, displacement_id), file)


def metadata2file(atoms, calc, obj, file="phono3py_trajectory.yaml"):
    """ save opt metadata to file """

    metadata = metadata2dict(atoms, calc, obj)

    to_yaml(metadata, file, mode="w")


def step2dict(atoms, calc, displacement_id):
    """ extract information from opt step and convet to plain dict """

    # info from opt algorithm
    obj_dict = {"id": displacement_id}

    return {"phono3py": obj_dict, **results2dict(atoms, calc, append_cell=False)}


def metadata2dict(atoms, calc, obj):
    """ convert metadata information to plain dict """

    prim_data = input2dict(atoms)

    obj_dict = {
        "supercell_matrix_2": obj._phonon_supercell_matrix.astype(int).tolist(),
        "supercell_matrix_3": obj.get_supercell_matrix().astype(int).tolist(),
        "symprec": float(obj._symprec),
        # "displacements": obj.get_displacements(),
        "displacement_dataset_2": obj.get_phonon_displacement_dataset(),
        "displacement_dataset_3": obj.get_displacement_dataset(),
        "primitive": prim_data["atoms"],
    }

    supercell_2 = to_Atoms(obj.get_phonon_supercell())
    supercell_3 = to_Atoms(obj.get_supercell())
    supercell_data = {
        "supercell_2": input2dict(supercell_2, calc)["atoms"],
        "supercell_3": input2dict(supercell_3, calc)["atoms"],
        "atoms": input2dict(supercell_3, calc)["atoms"],
        "calculator": input2dict(supercell_3, calc)["calculator"],
    }

    return {"phono3py": obj_dict, **supercell_data}
