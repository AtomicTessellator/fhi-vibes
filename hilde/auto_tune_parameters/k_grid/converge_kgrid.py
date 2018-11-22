from pathlib import Path
import numpy as np

from ase.constraints import UnitCellFilter
from ase.calculators.socketio import SocketIOCalculator

from hilde.auto_tune_params.k_grid.kpointoptimizer import KPointOptimizer
from hilde.helpers.paths import cwd
from hilde.trajectory.relaxation import metadata2file, step2file
from hilde.watchdogs import WallTimeWatchdog as Watchdog

def converge_kgrid(
    atoms,
    calc,
    func=lambda x: x.calc.get_property('energy',x)/len(x),
    loss_func=lambda x: x,
    dfunc_min=1e-4,
    even=True,
    unit_cell=True,
    maxsteps=100,
    trajectory="kpt_trajectory.yaml",
    logfile="kpoint_conv.log",
    socketio_port=None,
    walltime=1800,
    workdir=".",
    kpts_density_init=1.0,
):
    watchdog = Watchdog(walltime=walltime)

    workdir = Path(workdir).absolute()
    trajectory = Path(trajectory).absolute()

    kpt_settings = {
        "func": func,
        "loss_func": loss_func,
        "dfunc_min": dfunc_min,
        "even": even,
        "logfile": str(workdir / logfile),
    }

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc

    atoms.calc = calc
    opt_atoms = atoms

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        workdir / "calculation", mkdir=True
    ):
        if socketio_port is not None:
            atoms.calc = iocalc

        opt = KPointOptimizer(opt_atoms, **kpt_settings)

        # log very initial step and metadata
        if opt.nsteps == 0:
            metadata2file(atoms, calc, opt, trajectory)

        for ii, converged in enumerate(opt.irun(steps=maxsteps)):
            step2file(atoms, atoms.calc, opt, trajectory, unit_cell=unit_cell)
            if watchdog():
                break

    return converged, opt.kpts_density, calc
