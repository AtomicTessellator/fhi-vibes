"""Functions that generate FWActions after performing Aims Calculations"""
from pathlib import Path

import numpy as np


def get_calc_times(workdir=None, calc_dirs=None):
    """Get the calculation for a set of Aims Calculations"""
    if calc_dirs is None:
        if workdir:
            backup_dir = Path(workdir) / "backups"
        else:
            backup_dir = Path("./backups/")

        calc_dirs = list(backup_dir.glob("backup.?????/*"))
        if not calc_dirs:
            calc_dirs = list(backup_dir.glob("backup.?????.*"))
        print(calc_dirs)

    calc_times = list()
    for direc in calc_dirs:
        aims_out = Path(direc) / "aims.out"
        if aims_out.exists():
            lines = np.array(open(str(aims_out)).readlines())
            time_text = "          Detailed time accounting                     :  max(cpu_time)    wall_clock(cpu1)\n"
            time_line = np.where(lines == time_text)[0][0] + 1
            time = float(lines[time_line].split(":")[1].split("s")[0])
        else:
            time = 0.0
        calc_times.append(time)
    return calc_times
