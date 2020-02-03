import subprocess as sp
from pathlib import Path

parent = Path(__file__).parent


def test_vibes_help():
    cmd = "vibes --help"
    sp.run(cmd.split(), cwd=parent, check=True)


def test_vibes_output():
    cmd = "vibes output phon trajectory.son --full --q_mesh 5 5 5"
    sp.run(cmd.split(), cwd=parent, check=True)


def test_vibes_template():
    cmd = "vibes template aims"
    sp.run(cmd.split(), cwd=parent, check=True)


if __name__ == "__main__":
    test_vibes_help()
    test_vibes_output()
    test_vibes_template()
