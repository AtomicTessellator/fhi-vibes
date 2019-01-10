""" Provide a full highlevel phonopy workflow """

from pathlib import Path

from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid, k2d
from hilde.helpers.paths import cwd
from hilde.trajectory.phonons import metadata2dict
from hilde.tasks import calculate_socket, calc_dirname
from .postprocess import postprocess
from .wrapper import preprocess as phonopy_preprocess, defaults

def preprocess(
    atoms,
    calc,
    supercell_matrix=None,
    natoms_in_sc=None,
    kpt_density=None,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
):

    # Phonopy preprocess
    phonon, supercell, scs = phonopy_preprocess(
        atoms, supercell_matrix, natoms_in_sc, displacement, symprec
    )

    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    if kpt_density is not None:
        update_k_grid(supercell, calc, kpt_density)

    metadata = metadata2dict(atoms, calc, phonon)
    scs_return = []
    for sc in scs:
        if sc:
            sc.calc = calc
            scs_return.append(sc)

    return calc, supercell, scs, phonon, metadata


def run(
    atoms,
    calc,
    supercell_matrix=None,
    natoms_in_sc=None,
    kpt_density=None,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
    walltime=1800,
    workdir=".",
    trajectory="trajectory.yaml",
    primitive_file="geometry.in.primitive",
    supercell_file="geometry.in.supercell",
    force_constants_file="force_constants.dat",
    pickle_file="phonon.pick",
    db_kwargs=None,
    **kwargs,
):
    if "fireworks" in kwargs:
        fw = kwargs.pop("fireworks")
        if not kpt_density:
            kpt_density = k2d(atoms, calc.parameters["k_grid"])
    if "analysis_workdir" in kwargs:
        del(kwargs["analysis_workdir"])

    calc, supercell, scs, phonon, metadata = preprocess(
        atoms, calc, supercell_matrix, natoms_in_sc, kpt_density, displacement, symprec
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
            # phonon,
            trajectory=trajectory,
            workdir=workdir,
            force_constants_file=force_constants_file,
            pickle_file=pickle_file,
            db_kwargs=db_kwargs,
            **kwargs,
        )
        return True
