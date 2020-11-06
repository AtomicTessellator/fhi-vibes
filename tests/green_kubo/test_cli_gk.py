import subprocess as sp
from pathlib import Path

import pytest


parent = Path(__file__).parent


commands = [
    "vibes info vdos test.nc -p",
    "vibes output gk test.nc --filter_threshold 0.8",
    "vibes info gk -p",
]


@pytest.mark.parametrize("cmd", commands)
def test_cmd(cmd):
    sp.run(cmd.split(), cwd=parent, check=True)


if __name__ == "__main__":
    for cmd in commands:
        test_cmd(cmd)
