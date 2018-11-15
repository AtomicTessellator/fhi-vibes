import numpy as np

from ase.constraints import UnitCellFilter
from ase.optimize import BFGS

from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.paths import cwd


def relax(
    atoms,
    calc,
    fmax=0.01,
    unit_cell=True,
    maxsteps=100,
    trajectory="trajectory.yaml",
    socketio_port=None,
    walltime=1800,
    workdir=".",
):
    """ run bfgs relaxation  """

    raise Exception('under construction')

    bfgs_settings = {"logfile": "relax.log", "maxstep": 0.2}

    watchdog = Watchdog(walltime=walltime)

    trajectory = Path(trajectory).absolute()

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc
    atoms.calc = calc

    if unit_cell:
        opt_atoms = UnitCellFilter(atoms)
    else:
        opt_atoms = atoms

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        workdir, mkdir=True
    ):

        if socketio_port is not None:
            opt_atoms.calc = iocalc

        opt = BFGS(opt_atoms, **bfgs_settings)

        # log very initial step and metadata
        if md.nsteps == 0:
            metadata2file(atoms, calc, md, file=trajectory)
            step2file(atoms, atoms.calc, md, trajectory)

        for ii, converged in enumerate(opt.irun(fmax=fmax, steps=maxsteps)):
            if watchdog():
                break

    return converged
