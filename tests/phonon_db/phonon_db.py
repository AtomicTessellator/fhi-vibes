""" Test for the phonon database """

def test_phonon_db():
    from pathlib import Path

    import numpy as np

    from hilde.helpers.hash import hash_traj

    from hilde.phonon_db.phonon_db import connect
    from hilde.phonopy.postprocess import postprocess

    from hilde.trajectory import reader

    # Write the initial structure with the phonopy object
    phonon = postprocess(trajectory="trajectory.son")

    traj_hash, meta_hash = hash_traj(*reader("trajectory.son", True), True)
    hashes = {
        "traj_hash": traj_hash,
        "meta_hash": meta_hash,
    }

    # Get the row from the database
    for db_path in [Path("test.json")]:
        db = connect(db_path)
        db.write(phonon, hashes)
        row = list(
            db.select(
                sc_matrix_2=phonon.get_supercell_matrix(),
                traj_hash=hashes["traj_hash"],
                meta_hash=hashes["meta_hash"],
                has_fc2=True,
                columns=["id", "fc_2", "force_2"],
            )
        )[0]

        assert np.max(np.abs(row.fc_2[:] - phonon.get_force_constants()[:])) < 1e-12

        db_path.unlink()

if __name__ == "__main__":
    test_phonon_db()
