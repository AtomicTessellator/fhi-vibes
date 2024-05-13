import shutil
import subprocess as sp
from pathlib import Path

import yaml

from vibes.fireworks.launchpad import LaunchPad
from vibes.helpers.paths import cwd

parent = Path(__file__).parent


def test_fireworks_cli():
    lp = LaunchPad(strm_lvl="INFO")

    commands = (
        f"vibes fireworks add_wf -l my_launchpad.yaml -w {parent}/workflow_C.in",
        "vibes fireworks rlaunch -l ../my_launchpad.yaml rapidfire --max_loops 3",
    )

    with cwd(parent):
        with open("my_launchpad.yaml", "w") as lp_file:
            yaml.dump(lp.as_dict(), lp_file)

        sp.run(commands[0].split(), check=False)
        with cwd("fireworks_launchers", mkdir=True):
            sp.run(commands[1].split(), check=False)

        assert lp.get_wf_by_fw_id(2).state == "COMPLETED"

        shutil.rmtree("test_run/")
        shutil.rmtree("fireworks_launchers/")
        Path("my_launchpad.yaml").unlink()


if __name__ == "__main__":
    test_fireworks_cli()
