""" Provide a full highlevel phonopy workflow """

from pathlib import Path
import pickle
import numpy as np

from hilde.helpers.converters import dict2atoms
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import PhononRow
import hilde.phono3py.wrapper as ph3
from hilde.phonopy import displacement_id_str
from hilde.phono3py.wrapper import defaults
from hilde.structure.convert import to_Atoms
from hilde.trajectory import reader as traj_reader


def postprocess(
    phonon3,
    calculated_atoms=None,
    trajectory="trajectory.yaml",
    force_constants_file="force_constants.dat",
    pickle_file="phono3py.pick",
    **kwargs,
):
    """ Phonopy postprocess """
    if "workdir" in kwargs:
        workdir = kwargs["workdir"]
    else:
        workdir = "."

    if "fireworks" in kwargs and kwargs["fireworks"]:
        if isinstance(phonon3, dict):
            phonon3 = PhononRow(dct=phonon3).to_phonon3()
        else:
            phonon3 = PhononRow(phonon3=phonon3).to_phonon3()

        if "displacement" in kwargs:
            displacement = kwargs["displacement"]
        else:
            displacement = ph3.defaults["displacement"]

        if "cutoff_pair_distance" in kwargs:
            cutoff_pair_distance = kwargs["cutoff_pair_distance"]
        else:
            cutoff_pair_distance = ph3.defaults["cutoff_pair_distance"]

        phonon3.generate_displacements(
            distance=displacement,
            cutoff_pair_distance=cutoff_pair_distance,
            is_plusminus="auto",
            is_diagonal=True,
        )

    trajectory = Path(workdir) / trajectory
    if calculated_atoms:
        if "fireworks" in kwargs and kwargs["fireworks"]:
            temp_atoms = [dict2atoms(cell) for cell in calculated_atoms]
        else:
            temp_atoms = calculated_atoms.copy()
        calculated_atoms = sorted(
            temp_atoms,
            key=lambda x: x.info[displacement_id_str] if x else len(disp_cells) + 1,
        )
        fc3_cells = []
        used_forces = 0
        for cell in phonon3.get_supercells_with_displacements():
            if cell is not None:
                fc3_cells.append(calculated_atoms[used_forces])
                used_forces += 1
            else:
                fc3_forces.append(None)
    elif trajectory.is_file():
        fc3_cells = traj_reader(trajectory)
    else:
        raise ValueError("Either calculated_atoms or trajectory must be defined")

    fc3_forces = ph3.get_forces(fc3_cells)
    phonon3.produce_fc3(fc3_forces)

    with (Path(workdir) / pickle_file).open("wb") as fp:
        pickle.dump(phonon3, fp)

    print(
        f"*** Force Sets provided in {pickle_file}, "
        "full postprocess not yet supported"
    )
