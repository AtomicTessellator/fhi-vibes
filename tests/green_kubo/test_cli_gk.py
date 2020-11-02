import subprocess as sp
from pathlib import Path

import pytest


parent = Path(__file__).parent


commands = ["vibes info vdos test.nc -p"]  # , "vibes info gk -p"]
commands_fail = ["vibes output gk test.nc"]


@pytest.mark.parametrize("cmd", commands)
def test_cmd(cmd):
    sp.run(cmd.split(), cwd=parent, check=True)


@pytest.mark.parametrize("cmd", commands_fail)
def test_cmd_fail(cmd):
    with pytest.raises(sp.CalledProcessError):
        sp.run(cmd.split(), cwd=parent, check=True)


if __name__ == "__main__":
    for cmd in commands:
        test_cmd(cmd)
