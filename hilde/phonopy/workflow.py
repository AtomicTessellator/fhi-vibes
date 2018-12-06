""" Provide a full highlevel phonopy workflow """

from pathlib import Path

from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd
from hilde.trajectory.phonons import metadata2dict
from hilde.tasks import calculate_socket, calc_dirname
from .postprocess import postprocess
from .wrapper import preprocess as phonopy_preprocess, defaults

def preprocess(
    atoms,
    calc,
    supercell_matrix,
    kpt_density=None,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
):

    # Phonopy preprocess
    phonon, supercell, scs = phonopy_preprocess(
        atoms, supercell_matrix, displacement, symprec
    )

    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    if kpt_density is not None:
        update_k_grid(supercell, calc, kpt_density)

    metadata = metadata2dict(atoms, calc, phonon)

    for sc in scs:
        sc.calc = calc

    return calc, supercell, scs, phonon, metadata


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
    force_constants_file="force_constants.dat",
    pickle_file="phonon.pick",
    db_path=None,
    **kwargs,
):

    calc, supercell, scs, phonon, metadata = preprocess(
        atoms, calc, supercell_matrix, kpt_density, displacement, symprec
    )

    # save input geometries and settings
    settings = Settings()
    with cwd(workdir, mkdir=True):
        atoms.write(primitive_file, format="aims", scaled=True)
        supercell.write(supercell_file, format="aims", scaled=False)

    with cwd(Path(workdir) / calc_dirname, mkdir=True):
        settings.write()

    calculate_socket(
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

    postprocess(
        # phonon,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file=force_constants_file,
        pickle_file=pickle_file,
        db_kwargs={"db_path": db_path},
        **kwargs,
    )

    return True

# def preprocess_fireworks(atoms, calc, supercell_matrix, displacement=0.01):
#     """ phonopy preprocess returning supercells with attached calculator for FW """
#     phonon, supercell, scs = phonopy_preprocess(atoms, supercell_matrix, displacement)
#     if calc.name == "aims":
#         calc.parameters["compute_forces"] = True
#     for sc in scs:
#         sc.calc = calc
#     return phonon, supercell, scs, calc
