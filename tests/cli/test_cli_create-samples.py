from pathlib import Path
import subprocess as sp

from jinja2 import Template

parent = Path(__file__).parent


template = r"""hilde utils create-samples geometry.in.supercell
-fc FORCE_CONSTANTS_remapped
--zacharias
-T {{ temperature }}
-n {{ n_samples }}
-seed {{ seed }}
--propagate {{ propagate }}
"""

args = {"temperature": 600, "n_samples": 2, "seed": 4, "propagate": 1}

cmd = Template(template).render(args)


def test_run_cmd():
    """create samples with the cli tool"""
    sp.run(cmd.split(), cwd=parent)


def test_output():
    """check created samples vs reference"""
    for file in parent.glob("geometry.in.supercell.*"):
        reference = parent / "ref" / file.name

        t1, t2 = file.read_text().split("\n")[4:], reference.read_text().split("\n")[4:]

        assert t1 == t2, (t1, t2)
        file.unlink()


if __name__ == "__main__":
    test_run_cmd()
    test_output()
