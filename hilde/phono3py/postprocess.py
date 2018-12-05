""" Provide a full highlevel phonopy workflow """

from pathlib import Path
import pickle
import numpy as np

from hilde.helpers.converters import dict2atoms
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import PhononRow
import hilde.phonopy.wrapper as ph
from hilde.phonopy import displacement_id_str
from hilde.structure.convert import to_Atoms
from hilde.trajectory import reader as traj_reader


def postprocess(
    phonon3,
    trajectory="trajectory.yaml",
    force_constants_file="force_constants.dat",
    pickle_file="phono3py.pick",
    **kwargs,
):
    """ Phonopy postprocess """

    trajectory = Path(workdir) / trajectory

    if trajectory.is_file():
        calculated_atoms = traj_reader(trajectory)
    else:
        raise ValueError("Either calculated_atoms or trajectory must be defined")


    fc3_cells = []
    used_forces = 0
    for ii, cell in phonon3.get_supercells_with_displacements():
        if cell is not None and ii == :
            fc3_cells.append(calculated_atoms[used_forces])
            used_forces += 1
        else:
            fc3_forces.append(None)
    fc3_forces = ph3.get_forces(fc3_cells)
    phonon3.produce_fc3(fc3_forces)



    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    # compute and save force constants
    phonon3.set_fc3(force_sets)

    with (Path(workdir) / pickle_file).open("wb") as fp:
        pickle.dump(phonon, fp)

    exit(
        f"*** Force Sets provided in {pickle_file}, "
        "full postprocess not yet supported"
    )
