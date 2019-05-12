""" Provide a full highlevel phonopy workflow """

from pathlib import Path
import numpy as np
from hilde.helpers.converters import dict2results
from hilde.phonopy import displacement_id_str
from hilde.trajectory import reader as traj_reader
from hilde.helpers.pickle import psave

from hilde.phono3py.wrapper import prepare_phono3py
from hilde.phonopy.postprocess import postprocess as postprocess2


def postprocess(
    trajectory="phono3py/trajectory.yaml",
    trajectory_fc2="phonopy/trajectory.yaml",
    pickle_file="phonon3.pick",
    **kwargs,
):
    """ Phono3py postprocess """

    trajectory3 = Path(trajectory)

    # first run phonopy postprocess
    try:
        phonon = postprocess2(trajectory=trajectory_fc2)
    except FileNotFoundError:
        phonon = None

    # read the third order trajectory
    calculated_atoms, metadata_full = traj_reader(trajectory3, True)
    metadata = metadata_full["Phono3py"]

    phono3py_settings = {
        "atoms": dict2results(metadata["primitive"]),
        "supercell_matrix": metadata["supercell_matrix"],
        "phonon_supercell_matrix": phonon.get_supercell_matrix() if phonon else None,
        "fc2": phonon.get_force_constants() if phonon else None,
        "cutoff_pair_distance": metadata["displacement_dataset"]["cutoff_distance"],
        "symprec": metadata["symprec"],
        "displacement_dataset": metadata["displacement_dataset"],
        **kwargs
    }

    phonon3 = prepare_phono3py(**phono3py_settings)
    zero_force = np.zeros([len(calculated_atoms[0]), 3])

    # collect the forces and put zeros where no supercell was created
    force_sets = []
    disp_scells = phonon3.get_supercells_with_displacements()
    for nn, scell in enumerate(disp_scells):
        if scell:
            atoms = calculated_atoms.pop(0)
            assert atoms.info[displacement_id_str] == nn
            force_sets.append(atoms.get_forces())
        else:
            force_sets.append(zero_force)

    phonon3.produce_fc3(force_sets)

    if pickle_file:
        psave(phonon3, trajectory3.parent / pickle_file)

    return phonon3
