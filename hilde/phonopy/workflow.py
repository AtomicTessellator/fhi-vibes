""" Provide a full highlevel phonopy workflow

    Input: geometry.in and settings.in
    Output: geometry.in.supercell and trajectory.yaml """

from pathlib import Path

from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid, k2d
from hilde.helpers.paths import cwd
from hilde.tasks import calculate_socket, calc_dirname

# from .postprocess import postprocess
from .wrapper import preprocess, defaults
from . import metadata2dict


def run(
    atoms,
    calc,
    supercell_matrix,
    kpt_density=None,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
    walltime=1800,
    workdir=".",
    trajectory="trajectory.yaml",
    primitive_file="geometry.in.primitive",
    supercell_file="geometry.in.supercell",
    **kwargs,
):

    # Phonopy preprocess
    phonon, supercell, scs = preprocess(atoms, supercell_matrix, displacement, symprec)

    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    # update kpt density
    if kpt_density is not None:
        update_k_grid(supercell, calc, kpt_density)

    # save metadata
    metadata = metadata2dict(atoms, calc, phonon)

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
    )

    return completed
