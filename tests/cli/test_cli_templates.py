import subprocess as sp

import pytest

command = "vibes template "
templates = ["calculator aims", "calculator lj", "phonopy", "md", "relaxation"]


def _run(cmd, cwd):
    sp.run(cmd.split(), cwd=cwd, check=True)


@pytest.mark.parametrize("template", templates)
def test_cmd_tmpdir(template, tmp_path):
    _run(command + template, cwd=tmp_path)
