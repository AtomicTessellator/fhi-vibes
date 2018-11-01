""" Provide an aims calculator without much ado """

from pathlib import Path
from ase.calculators.aims import Aims
from hilde.settings import Settings

def setup_aims(custom_settings={}, workdir=None):
    """Set up an aims calculator.

    Args:
        custom_settings (dict): Settings that replace the minimal defaults
        workdir (str): directory to work in

    Returns:
        Aims: ASE calculator object

    """
    try:
        settings = Settings('../../hilde.conf')
        command = settings.machine.aims_command
        species_dir = Path(settings.machine.basissetloc) / custom_settings['species_type']
        del(custom_settings['species_type'])
    except:
        command = custom_settings['aims_command']
        species_dir = custom_settings['species_dir']
    default_settings = {
        'aims_command': command,
        'outfilename' : "aims.out",
        'species_dir': str(species_dir),
        'output_level': 'MD_light',
        'relativistic': 'atomic_zora scalar',
        'xc': 'pw-lda',
        'k_grid': 3 * [2],
    }
    aims_settings = {**default_settings, **custom_settings}

    if workdir:
        return Aims(label=Path(workdir).absolute(),
                    **aims_settings)
    else:
        return Aims(**aims_settings)
