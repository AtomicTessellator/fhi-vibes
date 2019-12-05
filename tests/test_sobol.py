import numpy as np
from vibes.helpers.sobol import RandomState


def test_sobol():
    rng = RandomState(dimension=3, seed=4)

    try:
        rng.rand(1)
        raise RuntimeError("this should not happen")
    except ValueError:
        pass

    rand = rng.rand(1, 3)
    norm = np.linalg.norm(rand - np.array([0.95140271, 0.63109021, 0.34202771]))

    assert norm < 1e-8, norm

if __name__ == "__main__":
    test_sobol()
