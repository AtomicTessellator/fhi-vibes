""" Provide a full highlevel phonopy workflow """
from pathlib import Path
import pickle

from hilde.helpers.converters import dict2results
from hilde.phonopy.wrapper import prepare_phonopy
from hilde.trajectory import reader


def postprocess(
    workdir=".", trajectory="trajectory.yaml", pickle_file="phonon.pick", **kwargs
):
    """ Phonopy postprocess """
    trajectory = Path(workdir) / trajectory

    calculated_atoms, metadata = reader(trajectory, True)

    primitive = dict2results(metadata["Phonopy"]["primitive"])
    supercell_matrix = metadata["Phonopy"]["supercell_matrix"]
    symprec = metadata["Phonopy"]["symprec"]

    phonon = prepare_phonopy(primitive, supercell_matrix, symprec=symprec)

    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    phonon.produce_force_constants(force_sets)

    if pickle_file:
        with (Path(workdir) / pickle_file).open("wb") as file:
            pickle.dump(phonon, file)

    return phonon
