import subprocess as sp
from pathlib import Path

import pytest

parent = Path(__file__).parent


commands = ["vibes aigk vdos test.nc -p", "vibes info hf test.nc -p"]


@pytest.mark.parametrize("cmd", commands)
def test_cmd(cmd):
    sp.run(cmd.split(), cwd=parent, check=True)


if __name__ == "__main__":
    for cmd in commands:
        test_cmd(cmd)
