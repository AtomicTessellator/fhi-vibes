"""test the MD workflow"""

import shutil
from pathlib import Path

import numpy as np
from ase.calculators.emt import EMT
from ase.io import read

from vibes.helpers import cwd
from vibes.relaxation.context import RelaxationContext
from vibes.settings import Settings

parent = Path(__file__).parent
workdir = parent / "relaxation"
logfile = workdir / "relaxation.log"
settings_file = parent / "relaxation.in"
ref_structure = read(parent / ".Al.in.ref", format="aims")


def test_relaxation(workdir=workdir):
    with cwd(parent):

        settings = Settings(settings_file=settings_file)

        if workdir.exists():
            shutil.rmtree(workdir)

        ctx = RelaxationContext(settings, workdir=workdir)
        ctx.atoms.rattle(0.01, seed=4)
        ctx.calc = EMT()

        ctx.run()
        ctx.run()

        # compare the relaxation path:
        log = logfile.read_text()
        ref_log = (parent / ".relaxation.log.ref").read_text()
        assert len(log) == len(ref_log), (log, ref_log)

        # compare final structure
        for pos, refpos in zip(ctx.atoms.positions, ref_structure.positions):
            assert np.allclose(pos, refpos)


if __name__ == "__main__":
    test_relaxation()
