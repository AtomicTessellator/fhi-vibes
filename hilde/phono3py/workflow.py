""" Provide a full highlevel phonopy workflow """

from pathlib import Path

from ase.calculators.socketio import SocketIOCalculator

from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd, move_to_dir
import hilde.phono3py.wrapper as ph3
from hilde.phonopy.workflow import calculate
from hilde.trajectory.phonons import metadata2dict, step2file, last_from_yaml
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

def run(
    atoms,
    calc,
    supercell_matrix,
    kpt_density=None,
    displacement=0.03,
    symprec=1e-5,
    q_mesh=[11,11,11],
    cutoff_pair_distance=10.0,
    log_level=2,
    trajectory="trajectory.yaml",
    walltime=1800,
    workdir=".",
    primitive_file="geometry.in.primitive",
    supercell_file="geometry.in.supercell",
    force_constants_file="force_constants.dat",
    pickle_file="phonon3.pick",
    backup_folder="backups",
    db_path=None,
    **kwargs,
):

    calculation_atoms, calc, supercell, scs, phonon3, metadata = preprocess(
        atoms, calc, supercell_matrix, kpt_density, displacement, symprec, q_mesh, cutoff_pair_distance, log_level
    )

    calculate(
        atoms=calculation_atoms,
        calculator=calc,
        primitive=atoms,
        supercell=supercell,
        supercells_with_displacements=scs,
        metadata=metadata,
        trajectory=trajectory,
        walltime=walltime,
        workdir=workdir,
        primitive_file=primitive_file,
        supercell_file=supercell_file,
    )

    postprocess(
        phonon3,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file=force_constants_file,
        pickle_file=pickle_file,
        db_path=db_path,
        **kwargs,
        # fingerprint_file=fingerprint_file,
    )

    return True

def preprocess(
    atoms, calc, supercell_matrix, kpt_density=None, displacement=0.03, symprec=1e-5, q_mesh=[11,11,11], cutoff_pair_distance=10.0, log_level=2, **kwargs
):
    phonon3, sc, scs = ph3.preprocess(
        atoms,
        supercell_matrix,
        q_mesh,
        displacement,
        cutoff_pair_distance,
        symprec,
        log_level,
    )

    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    calculation_atoms = scs[0].copy()

    if kpt_density is not None:
        update_k_grid(calculation_atoms, calc, kpt_density)

    metadata = metadata2dict(atoms, calc, phonon3)

    return calculation_atoms, calc, sc, scs, phonon3, metadata

def preprocess_fireworks(
    atoms,
    calc,
    supercell_matrix,
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
        supercell_matrix (list or np.ndarray): supercell matrix for the third order force constant calculations
        displacement (float): size of the displacements used in finite difference calculations
        q_mesh (list of ints): size of the qpoint mesh for thermal conductivity calcs
        cutoff_pair_distance (float): cutoff distance for force interactions
        symprec (float): symmetry percison for phono3py
        log_level (int): how much logging should be done
    """
    phonon3, sc, scs = ph3.preprocess(
        atoms,
        supercell_matrix,
        q_mesh,
        displacement,
        cutoff_pair_distance,
        symprec,
        log_level,
    )
    for sc in scs:
        sc.calc = calc
    return phonon3, sc, scs, calc
