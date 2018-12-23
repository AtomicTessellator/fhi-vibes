""" Provide a full highlevel phonopy workflow """
from pathlib import Path
import pickle
import numpy as np

from phonopy import Phonopy

from hilde.helpers.converters import dict2atoms, dict2results
from hilde import konstanten as const
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import PhononRow
import hilde.phonopy.wrapper as ph
from hilde.phonopy import displacement_id_str
from hilde.structure.convert import to_Atoms, to_phonopy_atoms
from hilde.trajectory import reader as traj_reader
from hilde.trajectory.phonons import step2file, to_yaml

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
            step2file(atoms, atoms.calc, nn, trajectory)


def postprocess(
    # phonon,
    metadata=None,
    calculated_atoms=None,
    trajectory="trajectory.yaml",
    workdir=".",
    force_constants_file="force_constants.dat",
    displacement=0.01,
    fireworks=False,
    pickle_file="phonon.pick",
    db_kwargs=None,
    **kwargs,
):
    """ Phonopy postprocess """
    trajectory = Path(workdir) / trajectory

    if fireworks:
        collect_forces_to_trajectory(trajectory, calculated_atoms, metadata)

    calculated_atoms, metadata = traj_reader(trajectory, True)
    ph_atoms = to_phonopy_atoms(dict2results(metadata["Phonopy"]["primitive"]), wrap=True)
    phonon = Phonopy(
        ph_atoms,
        supercell_matrix=np.array(metadata["Phonopy"]["supercell_matrix"]).reshape(3,3),
        is_symmetry=True,
        factor=const.omega_to_THz,
        **kwargs
    )
    phonon.set_displacement_dataset(metadata["Phonopy"]['displacement_dataset'])

    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    # compute and save force constants
    force_constants = ph.get_force_constants(phonon, force_sets)
    np.savetxt(Path(workdir) / force_constants_file, force_constants)

    # set the (default) bandstructure
    ph.get_bandstructure(phonon)

    with (Path(workdir) / pickle_file).open("wb") as fp:
        pickle.dump(phonon, fp)
    if db_kwargs is not None:
        db_kwargs = db_kwargs.copy()
        db_path = db_kwargs.pop("db_path")
        update_phonon_db(
            db_path,
            to_Atoms(phonon.get_unitcell()),
            phonon,
            symprec=phonon._symprec,
            sc_matrix_2=list(phonon.get_supercell_matrix().flatten()),
            **db_kwargs
        )