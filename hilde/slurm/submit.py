import time
import subprocess as sp

from hilde.helpers import talk
from .generate import generate_jobscript


def submit(
    dct, command=None, submit_command="sbatch", file="submit.sh", log=".submit.log"
):
    """submit the job described in dct"""

    if command:
        dct.update({"command": command})

    # write jobscribt to file
    generate_jobscript(dct, file=file)

    cmd = [submit_command, file]

    submit_output = sp.run(cmd, text=True, capture_output=True)

    if submit_output == "":
        submit_output = "empty (e.g. local computation)"

    talk(submit_output.stdout)

    try:
        timestr = time.strftime("%Y/%m/%d_%H:%M:%S")
        with open(log, "a") as f:
            f.write(f"{timestr}: {submit_output.stdout}\n")
            if submit_output.stderr:
                f.write(f"{timestr} [STDERR]: \n{submit_output.stderr}\n")
    except (IndexError, ValueError):
        print("Error during slurm submission: {:s}".format(submit_output))
