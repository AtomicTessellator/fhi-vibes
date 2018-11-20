from pathlib import Path
import numpy as np

from ase.constraints import UnitCellFilter
from ase.calculators.socketio import SocketIOCalculator

from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.paths import cwd
from hilde.trajectory.relaxation import metadata2file, step2file


def relax(
    atoms,
    calc,
    linesearch=True,
    fmax=0.01,
    maxstep=0.2,
    unit_cell=True,
    maxsteps=100,
    trajectory="bfgs_trajectory.yaml",
    logfile="relax.log",
    socketio_port=None,
    walltime=1800,
    workdir=".",
    output="geometry.in.relaxed",
):
    """ run a BFGS relaxation

    Args:
        atoms ([type]): [description]
        calc ([type]): [description]
        fmax (float, optional): Defaults to 0.01. [description]
        maxstep (float, optional): Defaults to 0.2. [description]
        unit_cell (bool, optional): Defaults to True. [description]
        maxsteps (int, optional): Defaults to 100. [description]
        trajectory (str, optional): Defaults to "bfgs_trajectory.yaml". [description]
        logfile (str, optional): Defaults to "relax.log". [description]
        socketio_port ([type], optional): Defaults to None. [description]
        walltime (int, optional): Defaults to 1800. [description]
        workdir (str, optional): Defaults to ".". [description]
        output (str, optional): Defaults to "geometry.in.relaxed". [description]

    Returns:
        [bool]: converged
    """

    if linesearch:
        from ase.optimize.bfgslinesearch import BFGSLineSearch as BFGS
    else:
        from ase.optimize.bfgs import BFGS

    watchdog = Watchdog(walltime=walltime)

    workdir = Path(workdir).absolute()
    trajectory = Path(trajectory).absolute()

    bfgs_settings = {"logfile": str(workdir / logfile), "maxstep": maxstep}

    if "compute_forces" in calc.parameters:
        calc.parameters["compute_forces"] = True

    if "compute_analytical_stress" in calc.parameters:
        if unit_cell:
            calc.parameters["compute_analytical_stress"] = True
        else:
            calc.parameters["compute_analytical_stress"] = False

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
        workdir / "calculation", mkdir=True
    ):
        if socketio_port is not None:
            atoms.calc = iocalc

        opt = BFGS(opt_atoms, **bfgs_settings)

        # log very initial step and metadata
        if opt.nsteps == 0:
            metadata2file(atoms, calc, opt, trajectory)

        for ii, converged in enumerate(opt.irun(fmax=fmax, steps=maxsteps)):
            step2file(atoms, atoms.calc, opt, trajectory, unit_cell=unit_cell)
            if watchdog():
                break

    with cwd(workdir):
        atoms.write(output, format="aims", scaled=True)

    return converged
