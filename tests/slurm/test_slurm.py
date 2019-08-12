from pathlib import Path

from hilde import Settings
from hilde.slurm.submit import submit

parent = Path(__file__).parent


ref_submit_script = """#!/bin/bash -l

#SBATCH -J test|h
#SBATCH -o log.test.%j
#SBATCH -e log.test.%j
#SBATCH -D ./
#SBATCH --mail-type=all
#SBATCH --mail-user=knoop@fhi-berlin.mpg.de
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-core=1
#SBATCH -t 0:30:00
#SBATCH --partition=express


hilde info geometry
"""


def test_submit():
    s = Settings(parent / "slurm.in").slurm

    logfile = parent / "test.log"
    submit_script = parent / "submit.sh"
    submit(s, submit_command="bash", log=logfile, file=submit_script)

    assert logfile.exists()
    assert submit_script.exists()
    assert submit_script.read_text() == ref_submit_script


if __name__ == "__main__":
    test_submit()
