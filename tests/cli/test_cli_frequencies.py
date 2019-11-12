from pathlib import Path
import subprocess as sp


import numpy as np
from ase.io import read

_filename = "frequencies.dat"
parent = Path(__file__).parent


cmd = r"""hilde utils fc frequencies"""


def test_run_cmd():
    """create samples with the cli tool"""
    sp.run(cmd.split(), cwd=parent)


def test_output():
    """check created frequencies vs reference"""
    file = parent / _filename
    frequencies = np.loadtxt(file)
    reference = np.loadtxt(parent / "ref" / _filename)

    for ii, (f1, f2) in enumerate(zip(frequencies, reference)):
        if abs(f1) < 1e-5:
            assert np.allclose(f1, f2, rtol=0.01), (f"Frequency {ii}: ", f1, f2)
        else:
            assert np.allclose(f1, f2), (f"Frequency {ii}: ", f1, f2)

    file.unlink()


if __name__ == "__main__":
    test_run_cmd()
    test_output()
