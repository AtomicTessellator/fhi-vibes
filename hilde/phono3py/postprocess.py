""" Provide a full highlevel phonopy workflow """

from pathlib import Path
import pickle
import numpy as np

from hilde.helpers.converters import dict2atoms
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import PhononRow
import hilde.phono3py.wrapper as ph3
from hilde.phonopy import displacement_id_str
from hilde.structure.convert import to_Atoms
from hilde.trajectory import reader as traj_reader


def postprocess(
    phonon3,
    calculated_atoms=None,
    trajectory="trajectory.yaml",
    workdir=".",
    force_constants_file_2="second_order_force_constants.dat",
    force_constants_file_3="third_order_force_constants.dat",
    displacement=0.03,
    q_mesh=[11,11,11],
    cutoff_pair_distance=10.0,
    symprec=1e-5,
    log_level=2,
    pickle_file="phonon3.pick",
    db_path=None,
    fireworks=False,
    **kwargs,
):
    """ Phonopy postprocess """

    if fireworks:
        if isinstance(phonon3, dict):
            phonon3 = PhononRow(dct=phonon3).to_phonon3()
        else:
            phonon3 = PhononRow(phonon3=phonon3).to_phonon3()
        phonon3.generate_displacements(
            distance=displacement,
            cutoff_pair_distance=cutoff_pair_distance,
            is_plusminus='auto',
            is_diagonal=True
        )
        if not phonon3._mesh:
            phonon3._mesh = np.array(q_mesh, dtype='intc')
    if calculated_atoms:
        if fireworks:
            temp_atoms = [dict2atoms(cell) for cell in calculated_atoms]
        else:
            temp_atoms = calculated_atoms.copy()
        calculated_atoms = sorted(
            temp_atoms, key=lambda x: x.info[displacement_id_str] if x else len(disp_cells) + 1
        )

    elif Path(trajectory).is_file():
        calculated_atoms = traj_reader(trajectory)
    else:
        raise ValueError("Either calculated_atoms or trajectory must be defined")

    fc2_forces = ph3.get_forces(calculated_atoms[:len(phonon3.get_phonon_supercells_with_displacements())])
    fc3_cells = []
    used_forces = len(fc2_forces)
    for cell in phonon3.get_supercells_with_displacements():
        if cell is not None:
            fc3_cells.append(calculated_atoms[used_forces])
            used_forces += 1
        else:
            fc3_forces.append(None)
    fc3_forces = ph3.get_forces(fc3_cells)
    phonon3.produce_fc2(fc2_forces)
    phonon3.produce_fc3(fc3_forces)

    phonon3.run_thermal_conductivity(write_kappa=True)

    # compute and save force constants
    n_atoms = phonon3.get_phonon_supercell().get_number_of_atoms()
    fc2 = phonon3.get_fc2().swapaxes(1, 2).reshape(2 * (3 * n_atoms,))
    np.savetxt(Path(workdir) / force_constants_file_2, fc2)
    n_atoms = phonon3.get_supercell().get_number_of_atoms()
    fc3 = phonon3.get_fc3().swapaxes(4, 3).swapaxes(4, 2).swapaxes(2,1).reshape(3 * (3 * n_atoms,))
    with open(str(Path(workdir) / force_constants_file_3), 'w') as outfile:
        for i,slice in enumerate(fc3):
            outfile.write(f"# New Slice Number {i}\n")
            np.savetxt(outfile, slice)
    with (Path(workdir) / pickle_file).open("wb") as fp:
        pickle.dump(phonon3, fp)
    if db_path:
        update_phonon3_db(
            db_path,
            to_Atoms(phonon3.get_unitcell()),
            phonon3,
            calc_type="phonon3",
            symprec=phonon3._symprec,
            sc_matrix_2=list(phonon3.get_phonon_supercell_matrix().flatten()),
            sc_matrix_3=list(phonon3.get_supercell_matrix().flatten()),
            **kwargs
        )