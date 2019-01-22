""" Provide a highlevel phono3py workflow for computing 3rd order force constants """

from pathlib import Path

from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid, k2d
from hilde.helpers.paths import cwd
from hilde.phonopy import metadata2dict
from hilde.tasks import calculate_socket, calc_dirname
from .postprocess import postprocess
from .wrapper import preprocess as phono3py_preprocess, defaults


def preprocess(
    atoms,
    calc,
    supercell_matrix=None,
    n_atoms_in_sc=None,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    kpt_density=None,
    up_kpoint_from_pc=False,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
):
    """ wrap the Phono3py preprocess for workflow """
    print(cutoff_pair_distance)
    # Phonopy preprocess
    phonon3, _, supercell, _, scs = phono3py_preprocess(
        atoms=atoms,
        fc3_supercell_matrix=supercell_matrix,
        n_atoms_in_sc_3=n_atoms_in_sc,
        disp=displacement,
        cutoff_pair_distance=cutoff_pair_distance,
        symprec=symprec,
    )
    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    if kpt_density is not None:
        update_k_grid(supercell, calc, kpt_density)
    elif up_kpoint_from_pc:
        kpt_density = k2d(atoms, calc.parameters['k_grid'])
        update_k_grid(supercell, calc, kpt_density)

    metadata = metadata2dict(atoms, calc, phonon3)

    to_run_scs = []
    for sc in scs:
        if sc:
            sc.calc = calc
            to_run_scs.append(sc)
    return calc, supercell, to_run_scs, phonon3, metadata


def run(
    atoms,
    calc,
    supercell_matrix=None,
    natoms_in_sc=None,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    kpt_density=None,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
    walltime=1800,
    workdir=".",
    trajectory="trajectory.yaml",
    primitive_file="geometry.in.primitive",
    supercell_file="geometry.in.supercell_fc3",
    force_constants_file="force_constants_3.dat",
    pickle_file="phono3py.pick",
    db_kwargs=None,
    up_kpoint_from_pc=False,
    **kwargs,
):
    '''
    Runs a Phono3py calculation
    Args:
        atoms (ASE Atoms Object): base atoms to use in the calculation
        calc (ACE Calculator Object): Calculator used to calculate the forces
        supercell_matrix (np.ndarray): Supercell matrix for the third order phonon calculation
        cutoff_pair_distance (float): Sets all interatomic interactions to zero for atoms farther apart than this number
        kpt_density (list of floats or float): k-point density for the calculation
        displacement (float): Finite displacement value
        symprec (float): Tolerance factor detecting the space group
        walltime (int): Max wall clock time for the calculation in seconds
        workdir (str): Path to the working directory
        trajectory (str): File name for the trajectory file
        primitive_file (str): File name for the primitive cell geometry file
        supercell_file (str): File name for the super cell geometry file
        force_constants_file (str): File name for the force constants output file
        pickle_file (str): File name for the pickled Phono3py object
        db_kwargs (dict): kwargs to add the calculation to a database
    Returns (bool): True if all calculations are completed
    '''
    if "fireworks" in kwargs:
        fw = kwargs.pop("fireworks")
        if not kpt_density:
            kpt_density = k2d(atoms, calc.parameters["k_grid"])
    if "analysis_workdir" in kwargs:
        del(kwargs["analysis_workdir"])

    calc, supercell, scs, phonon3, metadata = preprocess(
        atoms,
        calc,
        supercell_matrix,
        natoms_in_sc,
        cutoff_pair_distance,
        kpt_density,
        up_kpoint_from_pc,
        displacement,
        symprec,
    )

    # save input geometries and settings
    settings = Settings()
    with cwd(workdir, mkdir=True):
        atoms.write(primitive_file, format="aims", scaled=True)
        supercell.write(supercell_file, format="aims", scaled=False)

    with cwd(Path(workdir) / calc_dirname, mkdir=True):
        settings.write()

    completed = calculate_socket(
        atoms_to_calculate=scs,
        calculator=calc,
        metadata=metadata,
        trajectory=trajectory,
        walltime=walltime,
        workdir=workdir,
        primitive_file=primitive_file,
        supercell_file=supercell_file,
        **kwargs
    )
    if completed:
        postprocess(
            # phonon3,
            trajectory=trajectory,
            workdir=workdir,
            force_constants_file=force_constants_file,
            pickle_file=pickle_file,
            db_kwargs=db_kwargs,
            **kwargs,
        )

    return completed
