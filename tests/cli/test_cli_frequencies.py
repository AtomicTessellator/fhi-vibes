import subprocess as sp
from pathlib import Path

import numpy as np

_file = "frequencies.dat"
parent = Path(__file__).parent


cmd = r"""vibes utils fc frequencies"""


def test_output():
    """check created frequencies vs reference"""
    sp.run(cmd.split(), cwd=parent)
    file = parent / _file
    frequencies = np.loadtxt(file)
    reference = np.loadtxt(parent / "ref" / _file)

    for ii, (f1, f2) in enumerate(zip(frequencies, reference)):
        if abs(f1) < 1e-5:
            assert np.allclose(f1, f2, rtol=1), (f"Frequency {ii}: ", f1, f2)
        else:
            assert np.allclose(f1, f2), (f"Frequency {ii}: ", f1, f2)

    file.unlink()


if __name__ == "__main__":
    test_output()
