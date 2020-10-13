import subprocess as sp
import time
from pathlib import Path

from .generate import generate_jobscript


def submit(
    dct,
    command=None,
    submit_command="sbatch",
    file="submit.sh",
    submit_log=".submit.log",
    log_folder="log",
    dry=False,
):
    """submit the job described in dct"""
    Path(log_folder).mkdir(exist_ok=True)

    if command:
        dct.update({"command": command})

    dct.update({"logfile": str(Path(log_folder) / dct["name"])})

    # write jobscribt to file
    generate_jobscript(dct, file=file)

    if dry:
        print(f"DRY RUN requested: Jobscript written to {file}. STOP")
        return

    cmd = [submit_command, file]

    submit_err = None
    try:
        submit_output = sp.run(cmd, universal_newlines=True)
    except sp.CalledProcessError as err:
        submit_err = err.stderr

    if submit_output == "":
        submit_output = "empty (e.g. local computation)"

    try:
        timestr = time.strftime("%Y/%m/%d_%H:%M:%S")
        with open(submit_log, "a") as f:
            f.write(f"{timestr}: {submit_output}\n")
            if submit_err is not None:
                f.write(f"{timestr} [STDERR]: \n{submit_err}\n")
    except (IndexError, ValueError):
        print("Error during slurm submission: {:s}".format(submit_err))
