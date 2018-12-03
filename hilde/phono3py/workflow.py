""" Provide a full highlevel phonopy workflow """

from pathlib import Path

from ase.calculators.socketio import SocketIOCalculator

from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd, move_to_dir
import hilde.phono3py.wrapper as ph3
from hilde.trajectory.phono3py import metadata2file, step2file, last_from_yaml
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.compression import backup_folder as backup
from hilde.settings import DEFAULT_CONFIG_FILE
from .postprocess import postprocess


_calc_dirname = "calculations"


def last_calculation_id(trajectory):
    """ return the id of the last computed supercell """
    disp_id = -1

    try:
        dct = last_from_yaml(trajectory)
        disp_id = dct["phonopy"]["id"]
    except (FileNotFoundError, KeyError):
        pass

    return disp_id


def phono3py(
    atoms,
    calc,
    supercell_matrix_2,
    supercell_matrix_3,
    kpt_density=None,
    displacement=0.03,
    trajectory="trajectory.yaml",
    socketio_port=None,
    walltime=1800,
    workdir=".",
    primitive_file="geometry.in.primitive",
    supercell_2_file="geometry.in.supercell.fc2",
    supercell_3_file="geometry.in.supercell.fc3",
    force_constants_file_2="second_order_force_constants.dat",
    force_constants_file_3="third_order_force_constants.dat",
    pickle_file="phonon3.pick",
    backup_folder="backups",
    q_mesh=[11, 11, 11],
    cutoff_pair_distance=10.0,
    symprec=1e-5,
    log_level=2,
    db_path=None,
    **kwargs,
    #    fingerprint_file="fingerprint.dat",
):
    """
    Perform a full phono3py calculation
    Args:
        atoms (Atoms Object): Primitive cell to perform the phono3py calculation on
        calc (Calculator Object): ASE calculator used to get the forces of a system
        supercell_matrix_2 (list or np.ndarray): supercell matrix for the second order force constant calculations
        supercell_matrix_3 (list or np.ndarray): supercell matrix for the third order force constant calculations
        kpt_density (float): k-grid density to get a standardized k_grid for all calculations
        displacement (float): size of the displacements used in finite difference calculations
        trajectory (str): trajectory file for calculation
        socketio_port (int): socket number if using a socket io calculator
        walltime (int): Wall time of the calculation in seconds
        workdir (str): base work directory for the force calculations
        primitive_file (str): file name to store the primitive cell geometry
        supercell_2_file (str): file name to store the second order force constant super cell geometry
        supercell_3_file (str): file name to store the third order force constant super cell geometry
        force_constants_file_2 (str): Second order force constant output file name
        force_constants_file_3 (str): Third order force constant output file name
        pickle_file (str): pickle file filename
        backup_folder (str): directory to store back up information for a restart
        q_mesh (list of ints): size of the qpoint mesh for thermal conductivity calcs
        cutoff_pair_distance (float): cutoff distance for force interactions
        symprec (float): symmetry percison for phono3py
        log_level (int): how much logging should be done
        db_path (str): Path to database
    """

    workdir = Path(workdir)
    trajectory = (workdir / trajectory).absolute()
    backup_folder = workdir / backup_folder
    calc_dir = workdir / _calc_dirname

    watchdog = Watchdog(walltime=walltime, buffer=1)

    # backup configuration.cfg
    if workdir.absolute() != Path().cwd():
        move_to_dir(DEFAULT_CONFIG_FILE, workdir, exist_ok=True)

    # Phonopy preprocess
    phonon3, sc_2, sc_3, scs_2, scs_3 = ph3.preprocess(
        atoms,
        supercell_matrix_2,
        supercell_matrix_3,
        q_mesh,
        displacement,
        cutoff_pair_distance,
        symprec,
        log_level,
    )

    # make sure forces are computed (aims only)
    if "compute_forces" in calc.parameters:
        calc.parameters["compute_forces"] = True

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc

    # grab first geometry for computation
    calculation_atoms = scs_2[0].copy()

    if kpt_density is not None:
        update_k_grid(calculation_atoms, calc, kpt_density)

    # save input geometries and configuration
    with cwd(workdir, mkdir=True):
        atoms.write(primitive_file, format="aims", scaled=True)
        sc_2.write(supercell_2_file, format="aims", scaled=False)
        sc_3.write(supercell_3_file, format="aims", scaled=False)

    # perform calculation
    with cwd(calc_dir, mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            # save inputs and supercell
            if last_calculation_id(trajectory) < 0:
                metadata2file(atoms, calc, phonon3, trajectory)

            if socketio_port is not None:
                calc = iocalc

            for ii, cell in enumerate([*scs_2, *scs_3]):
                calculation_atoms.calc = calc
                if "forces" in calculation_atoms.calc.results:
                    del (calculation_atoms.calc.results["forces"])
                if cell is None:
                    continue

                # if precomputed:
                if ii <= last_calculation_id(trajectory):
                    continue

                # update calculation_atoms and compute force
                calculation_atoms.positions = cell.positions
                _ = calculation_atoms.get_forces()

                step2file(calculation_atoms, calculation_atoms.calc, ii, trajectory)

                if watchdog():
                    break

    # backup and restart if necessary
    if ii < len(scs_2) + len(scs_3) - 1:
        backup_folder(calc_dir, target_folder=backup_folder)
        return False

    postprocess(
        phonon3,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file_2=force_constants_file_2,
        force_constants_file_3=force_constants_file_3,
        pickle_file=pickle_file,
        db_path=db_path,
        **kwargs,
        # fingerprint_file=fingerprint_file,
    )

    return True


def initialize_phono3py_attach_calc(
    atoms,
    calc,
    supercell_matrix_2,
    supercell_matrix_3,
    displacement=0.03,
    q_mesh=[11, 11, 11],
    cutoff_pair_distance=10.0,
    symprec=1e-5,
    log_level=2,
):
    """
    Setup a full phono3py calculation
    Args:
        atoms (Atoms Object): Primitive cell to perform the phono3py calculation on
        calc (Calculator Object): ASE calculator used to get the forces of a system
        supercell_matrix_2 (list or np.ndarray): supercell matrix for the second order force constant calculations
        supercell_matrix_3 (list or np.ndarray): supercell matrix for the third order force constant calculations
        displacement (float): size of the displacements used in finite difference calculations
        q_mesh (list of ints): size of the qpoint mesh for thermal conductivity calcs
        cutoff_pair_distance (float): cutoff distance for force interactions
        symprec (float): symmetry percison for phono3py
        log_level (int): how much logging should be done
    """
    phonon3, sc_2, sc_3, scs_2, scs_3 = ph3.preprocess(
        atoms,
        supercell_matrix_2,
        supercell_matrix_3,
        q_mesh,
        displacement,
        cutoff_pair_distance,
        symprec,
        log_level,
    )
    for sc in [*scs_2, *scs_3]:
        sc.calc = calc
    return phonon3, sc_2, sc_3, scs_2, scs_3, calc
