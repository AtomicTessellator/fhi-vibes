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

    if command:
        dct.update({"command": command})

    dct.update({"logfile": str(Path(log_folder) / dct["name"])})

    # write jobscribt to file
    generate_jobscript(dct, file=file)

    if dry:
        print(f"DRY RUN requested: Jobscript written to {file}. STOP")
        return

    Path(log_folder).mkdir(exist_ok=True)

    cmd = [submit_command, file]

    submit_output = sp.run(cmd, text=True, capture_output=True)

    if submit_output == "":
        submit_output = "empty (e.g. local computation)"

    print(submit_output.stdout)

    try:
        timestr = time.strftime("%Y/%m/%d_%H:%M:%S")
        with open(submit_log, "a") as f:
            f.write(f"{timestr}: {submit_output.stdout}\n")
            if submit_output.stderr:
                f.write(f"{timestr} [STDERR]: \n{submit_output.stderr}\n")
    except (IndexError, ValueError):
        print("Error during slurm submission: {:s}".format(submit_output))
