""" Provide a full highlevel phonopy workflow """
from pathlib import Path

import numpy as np

from hilde.helpers.converters import dict2results
from hilde.phonopy.wrapper import prepare_phonopy, get_force_constants
from hilde.trajectory import reader
from hilde.helpers.pickle import psave
from hilde.io import write


def postprocess(
    trajectory="phonopy/trajectory.yaml", pickle_file="phonon.pick", **kwargs
):
    """ Phonopy postprocess """
    trajectory = Path(trajectory)

    calculated_atoms, metadata = reader(trajectory, True)

    primitive = dict2results(metadata["Phonopy"]["primitive"])
    supercell = dict2results(metadata["atoms"])
    supercell_matrix = metadata["Phonopy"]["supercell_matrix"]
    supercell.info = {"supercell_matrix": str(supercell_matrix)}
    symprec = metadata["Phonopy"]["symprec"]

    phonon = prepare_phonopy(primitive, supercell_matrix, symprec=symprec)
    phonon._displacement_dataset = metadata["Phonopy"]["displacement_dataset"].copy()

    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    phonon.produce_force_constants(force_sets)
    if pickle_file:
        psave(phonon, trajectory.parent / pickle_file)

    # Save the supercell
    fname = "geometry.in.supercell"
    write(supercell, fname)
    print(f"Supercell written to {fname}")

    force_constants = get_force_constants(phonon)
    fname = "force_constants.dat"
    np.savetxt(fname, force_costants)
    print(f"Force constants saved to {fname}.")

    return phonon
