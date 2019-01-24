""" Provide a full highlevel phonopy workflow

    Input: geometry.in and settings.in
    Output: geometry.in.supercell and trajectory.yaml """

from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.tasks import calculate_socket
from hilde.helpers.warnings import warn
from hilde.helpers.restarts import restart

from .postprocess import postprocess
from . import metadata2dict


def run_phonopy(**kwargs):
    """ high level function to run phonopy workflow """

    args = bootstrap(**kwargs)

    completed = calculate_socket(**args)

    if not completed:
        restart()
    else:
        print("Start postprocess.")
        postprocess(**args)
        print("done.")


def bootstrap(name="phonopy", **kwargs):
    """ load settings, prepare atoms, calculator, and phonopy """

    if name.lower() == "phonopy":
        from hilde.phonopy.wrapper import preprocess
    elif name.lower() == "phono3py":
        from hilde.phono3py.wrapper import preprocess

    settings = Settings()
    atoms = settings.get_atoms()

    if name not in settings:
        warn(f"Settings do not contain name instructions.", level=2)

    # Phonopy preprocess
    phonopy_settings = {"atoms": atoms, **settings[name], **kwargs}

    phonon, supercell, scs = preprocess(**phonopy_settings)

    if "calculator" in kwargs:
        calc = kwargs["calculator"]
    else:
        calc = setup_aims({"compute_forces": True}, settings=settings, atoms=supercell)

    # save metadata
    metadata = metadata2dict(phonon, calc)

    return {
        "atoms_to_calculate": scs,
        "calculator": calc,
        "metadata": metadata,
        **phonopy_settings,
    }
