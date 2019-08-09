"""test the MD workflow"""

import shutil
from pathlib import Path
import pytest
from hilde.settings import Settings
from hilde.helpers import cwd
from hilde.molecular_dynamics.context import MDContext


parent = Path(__file__).parent


@pytest.mark.filterwarnings("ignore:Subprocess")
def test_md():
    with cwd(parent):

        ctx = MDContext(Settings(settings_file="md.in"))

        workdir = ctx.workdir / "calculations"

        if workdir.exists():
            shutil.rmtree(workdir)

        try:
            ctx.run()
        except RuntimeError:
            assert (workdir / "control.in").exists()
            shutil.rmtree(ctx.workdir)


if __name__ == "__main__":
    test_md()
