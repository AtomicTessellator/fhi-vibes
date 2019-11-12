"""Tests the phonon_db interface"""


def test_phonon_db_interface():
    from pathlib import Path

    import numpy as np

    from hilde.phonon_db.database_interface import traj_to_database, from_database
    from hilde.phonopy.postprocess import postprocess

    db_path = Path("test.json")

    trajectory = Path(__file__).parent / "trajectory.son"

    hashes = traj_to_database(db_path, trajectory, True)

    selection = [(key, "=", val) for key, val in hashes.items()]
    ph = from_database(db_path, selection, get_phonon=True)
    ph_traj = postprocess(trajectory)

    assert (
        np.max(np.abs(ph.get_force_constants() - ph_traj.get_force_constants())) < 1e-12
    )

    db_path.unlink()


if __name__ == "__main__":
    test_phonon_db_interface()
