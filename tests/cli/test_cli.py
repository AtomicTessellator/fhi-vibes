import subprocess as sp
from pathlib import Path

import pytest


# from vibes.cli.scripts.nomad_upload import upload_folder_dry as nomad_upload_folder

parent = Path(__file__).parent


commands = (
    "vibes --help",
    "vibes info geometry geometry.in.primitive",
    "vibes utils geometry suggest-k-grid geometry.in.primitive",
    "vibes info phonopy",
    "vibes info relaxation relaxation.son -v",
)

gp = "geometry.in.primitive"
gs = "geometry.in.supercell"
commands_files = [
    ["vibes run singlepoint lj.in", "lj/trajectory.son"],
    ["vibes utils backup calculations", "backups/backup.00002.19636B4A.tgz"],
    ["vibes output phonopy phonopy.son --full --q_mesh 5 5 5", "output"],
    ["vibes utils trajectory 2db", "trajectory.db"],
    ["vibes utils trajectory 2tdep trajectory.son", "tdep"],
    ["vibes utils trajectory 2csv", "trajectory.csv"],
    ["vibes utils trajectory 2xyz", "trajectory.xyz"],
    ["vibes utils trajectory pick -n 1", "geometry.in.00001"],
    ["vibes utils hash trajectory.son", "hash.toml"],
    ["vibes utils fc frequencies", "frequencies.dat"],
]

commands_files_depend = [
    [
        [f"vibes utils geometry get-deformation {gp} {gs}", "deformation.dat"],
        [f"vibes utils geometry apply-deformation {gp}", f"{gp}.deformed"],
    ],
    [
        ["vibes output md", "trajectory.nc"],
        ["vibes output md", "trajectory.nc"],
        ["vibes output md --force", "trajectory.nc"],
    ],
]


def _run(cmd, cwd=parent):
    sp.run(cmd.split(), cwd=cwd, check=True)


def _exists(file, cwd=True):
    if cwd:
        assert (parent / file).exists()
    else:
        assert Path(file).exists()


@pytest.mark.parametrize("cmd", commands)
def test_cmd(cmd):
    """check command"""
    _run(cmd)


@pytest.mark.parametrize("cmd,file", commands_files)
def test_cmd_and_files(cmd, file):
    """check command and see if output file exists"""
    _run(cmd)
    _exists(file)


@pytest.mark.parametrize("group", commands_files_depend)
def test_cmd_and_files_depend(group):
    """check command and see if output file exists with dependencies"""
    for (c, f) in group:
        test_cmd_and_files(c, f)


if __name__ == "__main__":
    for cmd in commands:
        test_cmd(cmd)
    for (cmd, file) in commands_files:
        test_cmd_and_files(cmd, file)
