import subprocess as sp
from pathlib import Path

parent = Path(__file__).parent


def test_vibes_help():
    cmd = "vibes --help"
    sp.run(cmd.split(), cwd=parent, check=True)


def test_vibes_md_output():
    cmd = "vibes output md md.son"
    sp.run(cmd.split(), cwd=parent, check=True)


def test_vibes_phonopy_output():
    cmd = "vibes output phon phonopy.son --full --q_mesh 5 5 5"
    sp.run(cmd.split(), cwd=parent, check=True)


def test_vibes_template():
    cmd = "vibes template aims"
    sp.run(cmd.split(), cwd=parent, check=True)


if __name__ == "__main__":
    test_vibes_help()
    test_vibes_md_output()
    test_vibes_phonopy_output()
    test_vibes_template()
