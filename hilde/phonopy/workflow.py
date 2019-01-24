""" Provide a full highlevel phonopy workflow

    Input: geometry.in and settings.in
    Output: geometry.in.supercell and trajectory.yaml """

from pathlib import Path

from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd
from hilde.tasks import calculate_socket, calc_dirname
from hilde.helpers.warnings import warn
from hilde.helpers.restarts import restart

from .postprocess import postprocess
from .wrapper import preprocess, defaults
from . import metadata2dict


def run_phonopy(**kwargs):
    """ high level function to run phonopy workflow """

    args = bootstrap()
    args.update(kwargs)

    completed = calculate_socket(**args)

    if not completed:
        restart()
    else:
        print("Start postprocess.")
        postprocess(**args)
        print("done.")


def bootstrap():
    """ load settings, prepare atoms, calculator, and phonopy """

    settings = Settings()
    atoms = settings.get_atoms()

    if "phonopy" not in settings:
        warn("Settings do not contain phonopy instructions.", level=2)

    # Phonopy preprocess
    phonon, supercell, scs = preprocess(atoms, **settings.phonopy)

    calc = setup_aims({"compute_forces": True}, settings=settings, atoms=supercell)

    # save metadata
    metadata = metadata2dict(atoms, calc, phonon)

    return {
        "atoms_to_calculate": scs,
        "calculator": calc,
        "metadata": metadata,
        **settings.phonopy,
    }
