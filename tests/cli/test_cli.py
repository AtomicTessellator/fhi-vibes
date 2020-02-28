import subprocess as sp
from pathlib import Path

import pytest

from vibes.cli.scripts.nomad_upload import upload_folder_dry as nomad_upload_folder

parent = Path(__file__).parent


commands = (
    "vibes --help",
    "vibes info geometry geometry.in.primitive",
    "vibes utils suggest_k_grid geometry.in.primitive",
    "vibes info phonopy",
)

commands_files = [
    ["vibes utils backup calculations", "backups/backup.00002.19636B4A.tgz"],
    ["vibes output md", "trajectory.nc"],
    ["vibes output phon phonopy.son --full --q_mesh 5 5 5", "output"],
    ["vibes utils trajectory 2db", "trajectory.db"],
    ["vibes utils trajectory 2tdep", "tdep"],
    ["vibes utils trajectory pick_sample -n 1", "geometry.in.1"],
    ["vibes utils nomad upload calculations --token test --dry", nomad_upload_folder],
    ["vibes utils hash trajectory.son", "hash.toml"],
    ["vibes utils fc frequencies", "frequencies.dat"],
]

commands_tmpdir = ["vibes template aims"]


def _run(cmd, cwd=parent):
    sp.run(cmd.split(), cwd=cwd, check=True)


def _exists(file, cwd=True):
    if cwd:
        assert (parent / file).exists()
    else:
        assert Path(file).exists()


@pytest.mark.parametrize("cmd", commands)
def test_cmd(cmd):
    _run(cmd)


@pytest.mark.parametrize("cmd", commands_tmpdir)
def test_cmd_tmpdir(cmd, tmp_path):
    _run(cmd, cwd=tmp_path)


@pytest.mark.parametrize("cmd,file", commands_files)
def test_cmd_and_files(cmd, file):
    _run(cmd)
    _exists(file)


if __name__ == "__main__":
    for cmd in commands:
        test_cmd(cmd)
    for (cmd, file) in commands_files:
        test_cmd_and_files(cmd, file)
