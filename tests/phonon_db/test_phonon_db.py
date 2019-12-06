""" Test for the phonon database """


def test_phonon_db():
    from pathlib import Path

    import numpy as np

    from vibes.helpers.hash import hash_traj

    from vibes.phonon_db.phonon_db import connect
    from vibes.phonopy.postprocess import postprocess

    from vibes.trajectory import reader

    # Write the initial structure with the phonopy object
    trajectory = Path(__file__).parent / "trajectory.son"
    phonon = postprocess(trajectory=trajectory)

    traj_hash, meta_hash = hash_traj(*reader(trajectory, True), True)
    hashes = {"traj_hash": traj_hash, "meta_hash": meta_hash}

    # Get the row from the database
    db_path = Path(__file__).parent / "test.json"
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
