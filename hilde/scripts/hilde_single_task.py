""" create a configuration file and working directory """

import shutil
from argparse import ArgumentParser
from pathlib import Path
from hilde.settings import Settings, default_config_name
from hilde import supported_tasks


def main():
    """ main routine """
    parser = ArgumentParser(description="create a configuration file and workdir")
    parser.add_argument("config_files", nargs="+", help="hilde.cfg etc.")
    parser.add_argument("-g", "--geometry", help="geometry file to use")
    parser.add_argument("-wd", "--workdir", help="name of the working directory")
    parser.add_argument("--dry", action="store_true")
    args = parser.parse_args()

    settings = Settings(args.config_files, write=False)
    print("Summary of settings:")
    settings.print()

    for task in supported_tasks:
        if task in settings:
            break
        exit("Task could not be identified.")

    print(f"Task to be performed:  {task}")

    if "workdir" in settings[task]:
        workdir = settings[task].pop("workdir")

    if args.workdir:
        workdir = args.workdir

    print(f"Working directory:     {workdir}")

    if Path(workdir).exists():
        exit("**Error: Working directory exsits, chance of dataloss.")
    else:
        Path(workdir).mkdir()

    if args.dry:
        exit()

    outfile = Path(workdir) / default_config_name
    settings.write(filename=outfile)
    print(f"Settings written to:   {outfile}")

    if args.geometry:
        outfile = Path(workdir) / args.geometry
        shutil.copy(args.geometry, outfile)
        print(f"Geometry written to:   {outfile}")

    # copy run script
    run_script = Path(settings.common.home_dir) / f'hilde/scripts/run/{task}.py'
    shutil.copy(run_script, Path(workdir) / 'run.py')
    print(f"Run script written to: {workdir}/run.py")



if __name__ == "__main__":
    main()
