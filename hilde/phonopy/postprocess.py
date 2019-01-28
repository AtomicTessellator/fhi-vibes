""" Provide a full highlevel phonopy workflow """
from pathlib import Path

from hilde.helpers.converters import dict2results
from hilde.phonopy.wrapper import prepare_phonopy
from hilde.trajectory import reader
from hilde.helpers.pickle import psave


def postprocess(
    trajectory="phonopy/trajectory.yaml", pickle_file="phonon.pick", **kwargs
):
    """ Phonopy postprocess """
    trajectory = Path(trajectory)

    calculated_atoms, metadata = reader(trajectory, True)

    primitive = dict2results(metadata["Phonopy"]["primitive"])
    supercell_matrix = metadata["Phonopy"]["supercell_matrix"]
    symprec = metadata["Phonopy"]["symprec"]

    phonon = prepare_phonopy(primitive, supercell_matrix, symprec=symprec)

    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    phonon.produce_force_constants(force_sets)

    if pickle_file:
        psave(phonon, trajectory.parent / pickle_file)

    return phonon
