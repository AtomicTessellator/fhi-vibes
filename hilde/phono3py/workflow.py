""" Provide a highlevel phono3py workflow for computing 3rd order force constants """

from pathlib import Path

from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid, k2d
from hilde.helpers.paths import cwd
from hilde.helpers.config import AttributeDict as adict
from hilde.trajectory.phonons import metadata2dict
from hilde.tasks import calculate_socket, calc_dirname
from .postprocess import postprocess
from .wrapper import preprocess as phono3py_preprocess, defaults


def preprocess(
    atoms,
    calc,
    supercell_matrix,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    kpt_density=None,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
):

    """ wrap the Phono3py preprocess for workflow """

    # Phonopy preprocess
    phonon3, _, supercell, _, scs = phono3py_preprocess(
        atoms=atoms,
        fc3_supercell_matrix=supercell_matrix,
        disp=displacement,
        cutoff_pair_distance=cutoff_pair_distance,
        symprec=symprec,
    )
    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    if kpt_density is not None:
        update_k_grid(supercell, calc, kpt_density)

    metadata = metadata2dict(atoms, calc, phonon3)

    for sc in scs:
        if sc is not None:
            sc.calc = calc

    return calc, supercell, scs, phonon3, metadata


def run(
    atoms,
    calc,
    supercell_matrix,
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
    **kwargs,
):
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
        cutoff_pair_distance,
        kpt_density,
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
