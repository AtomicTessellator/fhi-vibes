""" Provide highlevel access to calculator.calculate """
from .calculate import calculate_socket, calc_dirname
from vibes.settings import Settings
from vibes.helpers.utils import talk


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
