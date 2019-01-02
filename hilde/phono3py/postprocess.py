""" Provide a full highlevel phonopy workflow """

from pathlib import Path
from phono3py.phonon3 import Phono3py
import pickle
import numpy as np
import json
from phono3py.file_IO import write_fc3_to_hdf5
from hilde.helpers.converters import dict2atoms, dict2results
from hilde import konstanten as const
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import PhononRow
import hilde.phono3py.wrapper as ph3
from hilde.phonopy import displacement_id_str
from hilde.phono3py.wrapper import defaults
from hilde.structure.convert import to_Atoms, to_phonopy_atoms
from hilde.trajectory import reader as traj_reader
from hilde.trajectory.phonons import step2file, to_yaml

def collect_forces_to_trajectory(
    trajectory,
    calculated_atoms,
    metadata,
):
    Path(trajectory).parents[0].mkdir(exist_ok=True, parents=True)
    for el in metadata["Phono3py"]["displacement_dataset"]["first_atoms"]:
        el["number"] = int(el["number"])
        for d, dd in enumerate(el["direction"]):
            el['direction'][d] = int(dd)
        for el2 in el["second_atoms"]:
            el2["number"] = int(el2["number"])

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
    # phonon3,
    metadata=None,
    calculated_atoms=None,
    symmetrize_fc3=False,
    workdir=".",
    trajectory="trajectory.yaml",
    force_constants_file="force_constants_3.hdf5",
    pickle_file="phono3py.pick",
    db_kwargs=None,
    fireworks=False,
    **kwargs,
):
    trajectory = Path(workdir) / trajectory
    force_constants_file = Path(workdir) / force_constants_file
    if fireworks:
        collect_forces_to_trajectory(trajectory, calculated_atoms, metadata)

    calculated_atoms, metadata = traj_reader(trajectory, True)
    # if not phonon3:
    ph_atoms = to_phonopy_atoms(dict2results(metadata["Phono3py"]["primitive"]), wrap=True)
    phonon3_kwargs = kwargs.copy()
    if "fc2" in phonon3_kwargs:
        del(phonon3_kwargs['fc2'])
    phonon3 = Phono3py(
        ph_atoms,
        supercell_matrix=np.array(metadata["Phono3py"]["supercell_matrix"]).reshape(3,3),
        is_symmetry=True,
        frequency_factor_to_THz=const.omega_to_THz,
        **phonon3_kwargs,
    )
    phonon3.set_displacement_dataset(metadata["Phono3py"]['displacement_dataset'])
    if "fc2" in kwargs:
        phonon3.set_fc2(fc2)
    forces_shape = calculated_atoms[0].get_forces().shape

    force_sets_fc3 = []
    used_forces = 0

    for ii, cell in enumerate(phonon3.get_supercells_with_displacements()):
        if cell is not None:
            ref_atoms = calculated_atoms[used_forces]
            if ref_atoms.info[displacement_id_str] == ii:
                force_sets_fc3.append(ref_atoms.get_forces())
                used_forces += 1
        else:
            force_sets_fc3.append(np.zeros(forces_shape))

    # compute force constants
    if symmetrize_fc3:
        phonon3.produce_fc3(
            force_sets_fc3,
            is_translational_symmetry=True,
            is_permutation_symmetry=True,
            is_permutation_symmetry_fc2=True,
        )
    else:
        phonon3.produce_fc3(force_sets_fc3)

    write_fc3_to_hdf5(phonon3.get_fc3(), filename=force_constants_file)

    if db_kwargs is not None:
        print(db_kwargs)
        db_path = db_kwargs.pop("db_path")
        update_phonon_db(
            db_path,
            to_Atoms(phonon3.get_unitcell()),
            phonon3,
            symprec=phonon3._symprec,
            sc_matrix_3=list(phonon3.get_supercell_matrix().flatten()),
            **db_kwargs
        )

