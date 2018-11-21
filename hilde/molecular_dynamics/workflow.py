""" run molecular dynamics simulations using the ASE classes """

import shutil
from pathlib import Path
from subprocess import run
from ase.calculators.socketio import SocketIOCalculator
from hilde.trajectory.md import step2file, metadata2file
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.paths import cwd
from hilde.helpers.compression import backup_folder


_calc_dirname = "calculations"


def run_md(
    atoms,
    calc,
    md,
    maxsteps=25000,
    trajectory=None,
    metadata="md_metadata.yaml",
    socketio_port=None,
    walltime=1800,
    workdir=".",
    logfile=None,
    restart_cmd=None,
    **kwargs,
):
    """ run and MD for a specific time

    Args:
        atoms (Atoms): initial structure
        calc (Calculator): calculator to calculate forces (and stresses)
        md (MolecularDynamics): MD algorithm
        trajectory (str/Path): use to save trajectory as yaml
        socketio_port (int): use to run through socketio
        walltime (int, optional): Defaults to 1800.
        workdir (str, optional): Defaults to '.'.
    """

    watchdog = Watchdog(walltime=walltime)

    if trajectory is None:
        trajectory = (Path(workdir) / 'trajectory.yaml').absolute()
    else:
        trajectory = Path(trajectory).absolute()
    calc_dir = Path(workdir) / _calc_dirname

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc
    atoms.calc = calc

    # backup
    additional_files = []
    if logfile is not None:
        additional_files = [logfile]

    if calc_dir.exists():
        backup_folder(
            calc_dir,
            target_folder=workdir,
            additional_files=additional_files,
        )

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        calc_dir, mkdir=True
    ):

        if socketio_port is not None:
            atoms.calc = iocalc

        # log very initial step and metadata
        if md.nsteps == 0:
            metadata2file(atoms, calc, md, file=trajectory)
            step2file(atoms, atoms.calc, md, trajectory)

        # store MD metadata locally
        metadata2file(atoms, calc, md, file=metadata)

        while not watchdog() and md.nsteps < maxsteps:
            md.run(1)
            step2file(atoms, atoms.calc, md, trajectory)

    # backup and cleanup
    backup_folder(
        calc_dir,
        target_folder=workdir,
        additional_files=additional_files,
    )

    shutil.rmtree(calc_dir)

    # restart
    if md.nsteps < maxsteps:
        if restart_cmd is not None:
            print(f"restart with {restart_cmd}")
            run(restart_cmd.split())
    else:
        return True
