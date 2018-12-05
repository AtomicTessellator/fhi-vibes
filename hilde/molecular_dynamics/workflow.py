""" run molecular dynamics simulations using the ASE classes """

import shutil
from pathlib import Path
from ase.calculators.socketio import SocketIOCalculator

from hilde.settings import Settings
from hilde.trajectory.md import step2file, metadata2file
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.paths import cwd
from hilde.helpers.socketio import get_port
from hilde.helpers.compression import backup_folder as backup
from .initialization import setup_md


_calc_dirname = "calculations"


def run_md(
    atoms,
    calc,
    md=None,
    restart=True,
    maxsteps=25000,
    trajectory="trajectory.yaml",
    metadata="md_metadata.yaml",
    walltime=1800,
    workdir=".",
    backup_folder="backups",
    logfile="md.log",
    wf_extrapolation=True,
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

    watchdog = Watchdog(walltime=walltime, **kwargs)

    workdir = Path(workdir)
    trajectory = (workdir / trajectory).absolute()
    logfile = workdir / logfile
    calc_dir = workdir / _calc_dirname
    backup_folder = workdir / backup_folder

    # take the literal settings for running the task
    settings = Settings()

    # make sure forces are computed (aims only)
    # use wavefunction extrapolation
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

        if wf_extrapolation:
            calc.parameters["wf_extrapolation"] = "quadratic"

    if md is None:
        md_settings = {
            **kwargs,
            "logfile": str(logfile),
            "trajectory": trajectory,
            "workdir": workdir,
        }
        atoms, md, _ = setup_md(atoms, **md_settings)

    if md is None:
        raise RuntimeError("ASE MD algorithm has to be given")

    socketio_port = get_port(calc)
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
        backup(calc_dir, target_folder=backup_folder, additional_files=additional_files)

    something_happened = False
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
        settings.write()

        while not watchdog() and md.nsteps < maxsteps:
            something_happened = True
            md.run(1)
            step2file(atoms, atoms.calc, md, trajectory)

    # backup and cleanup if something new happened
    if something_happened:
        backup(calc_dir, target_folder=backup_folder, additional_files=additional_files)
        shutil.rmtree(calc_dir)

    # restart
    if md.nsteps < maxsteps:
        return False
    return True
