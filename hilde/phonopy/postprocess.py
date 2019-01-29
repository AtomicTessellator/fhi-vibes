""" Provide a full highlevel phonopy workflow """
from pathlib import Path

from hilde.helpers.converters import dict2results
from hilde.phonopy.wrapper import prepare_phonopy
from hilde.trajectory import reader
from hilde.helpers.pickle import psave

def collect_forces_to_trajectory(
    trajectory,
    calculated_atoms,
    metadata,
):
    Path(trajectory).parents[0].mkdir(exist_ok=True, parents=True)
    for el in metadata["Phonopy"]["displacement_dataset"]["first_atoms"]:
        el["number"] = int(el["number"])

    to_yaml(metadata, trajectory, mode="w")

    if isinstance(calculated_atoms[0], dict):
        temp_atoms = [dict2atoms(cell) for cell in calculated_atoms]
    else:
        temp_atoms = calculated_atoms.copy()
    calculated_atoms = sorted(
        temp_atoms,
        key=lambda x: x.info[displacement_id_str] if x else len(calculated_atoms) + 1,
    )
    for nn, atoms in enumerate(calculated_atoms):
        if atoms:
            step2file(atoms, atoms.calc, trajectory)


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
    print(trajectory.parent, pickle_file)
    if pickle_file:
        psave(phonon, trajectory.parent / pickle_file)

    return phonon
