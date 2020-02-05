import subprocess as sp
from pathlib import Path

import pytest

parent = Path(__file__).parent


commands = (
    "vibes --help",
    "vibes output md",
    "vibes output phon phonopy.son --full --q_mesh 5 5 5",
    "vibes template aims",
    "vibes info geometry geometry.in.primitive",
    "vibes utils suggest_k_grid geometry.in.primitive",
    "vibes utils trajectory 2db",
    "vibes utils trajectory 2tdep",
    "vibes utils trajectory pick_sample -n 1",
)


@pytest.mark.parametrize("cmd", commands)
def test_cmd(cmd):
    sp.run(cmd.split(), cwd=parent, check=True)


if __name__ == "__main__":
    for cmd in commands:
        test_cmd(cmd)
