""" Provide a full highlevel phonopy workflow """

from pathlib import Path

from ase.calculators.socketio import SocketIOCalculator

from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd
import hilde.phonopy.wrapper as ph
from hilde.trajectory.phonopy import metadata2file, step2file
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from .postprocess import postprocess


def phonopy(
    atoms,
    calc,
    supercell_matrix,
    kpt_density=None,
    displacement=0.01,
    trajectory="phonopy_trajectory.yaml",
    socketio_port=None,
    walltime=1800,
    workdir=".",
    primitive_file="geometry.in.primitive",
    supercell_file="geometry.in.supercell",
    force_constants_file="force_constants.dat",
    bandstructure_file="bandstructure.pdf",
    pickle_file="phonon.pick",
    db_path=None,
    **kwargs,
    #    fingerprint_file="fingerprint.dat",
):
    """ perform a full phonopy calculation """

    trajectory = Path(trajectory).absolute()
    workdir = Path(workdir).absolute()

    watchdog = Watchdog(walltime=walltime, buffer=1)

    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)

    if "compute_forces" in calc.parameters:
        calc.parameters["compute_forces"] = True

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc

    # grab first geometry for computation
    calculation_atoms = scs[0].copy()

    if kpt_density is not None:
        update_k_grid(calculation_atoms, calc, kpt_density)

    # save input geometries
    with cwd(workdir, mkdir=True):
        atoms.write(primitive_file, format="aims", scaled=True)
        supercell.write(supercell_file, format="aims", scaled=True)

    # perform calculation
    with cwd(workdir / "calculations", mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            # save inputs and supercell
            metadata2file(atoms, calc, phonon, trajectory)

            if socketio_port is not None:
                calc = iocalc

            for ii, cell in enumerate(scs):
                calculation_atoms.calc = calc

                if cell is None:
                    continue

                # update calculation_atoms and compute force
                calculation_atoms.positions = cell.positions
                _ = calculation_atoms.get_forces()

                step2file(calculation_atoms, calculation_atoms.calc, ii, trajectory)

                if watchdog():
                    raise ResourceWarning("**Watchdog: running out of time!")

    postprocess(
        phonon,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file=force_constants_file,
        bandstructure_file=bandstructure_file,
        pickle_file=pickle_file,
        db_path=db_path,
        **kwargs,
        # fingerprint_file=fingerprint_file,
    )


def initialize_phonopy_attach_calc(atoms, calc, supercell_matrix, displacement=0.01):
    """ phonopy preprocess returning supercells with attached calculator for FW """
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)
    for sc in scs:
        sc.calc = calc
    return phonon, supercell, scs
