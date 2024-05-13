import subprocess as sp
from pathlib import Path

import pytest

parent = Path(__file__).parent


commands = (
    "vibes info anharmonicity trajectory.nc",
    "vibes info anharmonicity trajectory.nc --per_sample",
    "vibes info anharmonicity trajectory.nc --per_mode --dry",
)


@pytest.mark.parametrize("cmd", commands)
def test_cmd(cmd):
    sp.run(cmd.split(), cwd=parent, check=True)


if __name__ == "__main__":
    for cmd in commands:
        test_cmd(cmd)
