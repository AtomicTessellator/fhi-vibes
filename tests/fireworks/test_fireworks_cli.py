import shutil
import subprocess as sp
from pathlib import Path

from vibes.fireworks.launchpad import LaunchPad
from vibes.helpers.paths import cwd

parent = Path(__file__).parent


def test_fireworks_cli():
    lp = LaunchPad(strm_lvl="INFO")
    lp.reset("", require_password=False)

    commands = (
        f"vibes fireworks add_wf -w {parent}/workflow_C.in",
        "vibes fireworks rlaunch rapidfire --max_loops 3",
    )

    with cwd(parent):
        sp.run(commands[0].split())
        with cwd("fireworks_launchers", mkdir=True):
            sp.run(commands[1].split())

        assert lp.get_wf_by_fw_id(1).state == "COMPLETED"

        shutil.rmtree("test_run/")
        shutil.rmtree("fireworks_launchers/")


if __name__ == "__main__":
    test_fireworks_cli()
