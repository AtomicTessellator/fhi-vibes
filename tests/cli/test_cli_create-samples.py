from pathlib import Path
import subprocess as sp

from jinja2 import Template

import numpy as np
from ase.io import read

parent = Path(__file__).parent


template = r"""hilde utils create-samples geometry.in.supercell
-fc FORCE_CONSTANTS_remapped
--zacharias
-T {{ temperature }}
-n {{ n_samples }}
-seed {{ seed }}
--propagate {{ propagate }}
"""

args = {"temperature": 600, "n_samples": 2, "seed": 4, "propagate": 1}

cmd = Template(template).render(args)


def test_run_cmd():
    """create samples with the cli tool"""
    sp.run(cmd.split(), cwd=parent)


def test_output():
    """check created samples vs reference"""
    files = parent.glob("geometry.in.supercell.*")
    for file in files:
        reference = parent / "ref" / file.name

        atoms = read(file, format="aims")
        ref_atoms = read(reference, format="aims")

        # positions
        pos1, pos2 = atoms.get_positions(), ref_atoms.get_positions()
        for ii, (p1, p2) in enumerate(zip(pos1, pos2)):
            assert np.allclose(p1, p2), (f"Position {ii}: ", p1, p2)

        # velocities
        vs1, vs2 = atoms.get_velocities(), ref_atoms.get_velocities()
        for ii, (p1, p2) in enumerate(zip(vs1, vs2)):
            assert np.allclose(p1, p2), (f"Velocity {ii}: ", p1, p2)

        # header
        idx = slice(4, 17)
        h1 = file.read_text().split("\n")[idx]
        h2 = reference.read_text().split("\n")[idx]

        for ii, (line1, line2) in enumerate(zip(h1, h2)):
            assert line1 == line2, (f"Header, line {ii}: ", line1, line2)

    for file in files:
        file.unlink()


if __name__ == "__main__":
    test_run_cmd()
    test_output()
