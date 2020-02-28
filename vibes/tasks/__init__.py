""" Provide highlevel access to calculator.calculate """
from vibes.helpers.utils import talk
from vibes.settings import Settings

from .calculate import calc_dirname, calculate_socket  # noqa: F401


def run():
    """ loader for vibes workflows:
            - phonopy
            - md """

    # load settings
    settings = Settings()

    if "phonopy" in settings:
        from vibes.phonopy import run_phonopy

        talk("launch phonoy workflow")

        run_phonopy()

    elif "md" in settings:
        from vibes.molecular_dynamics import run_md

        talk("launch MD workflow")

        run_md()
