'''Set up and run a harmonic analysis'''
# from hilde.anharmonicity_score import get_anharmonicity_score
from hilde.helpers.converters import calc2dict
from hilde.helpers.warnings import warn
from hilde.helpers.restarts import restart
from hilde.settings import TaskSettings, Settings
from hilde.statistical_sampling.initialization import preprocess
from hilde.tasks.calculate import calculate_socket
from hilde.templates.aims import setup_aims


def run_statistical_sampling(**kwargs):
    """ high level function to run harmonic analysis workflow """

    args = bootstrap(**kwargs)

    try:
        get_anharmonicity_score(**args)
        exit("** Postprocess could be performed from previous calculations. Check!")
    except (FileNotFoundError, RuntimeError):
        completed = calculate_socket(**args)

    if not completed:
        restart()
    # else:
    #     print("Start postprocess.")
    #     get_anharmonicity_score(**args)
    #     print("done.")


def bootstrap(name="statistical_sampling", settings=None, **kwargs):
    """ load settings, prepare atoms, calculator, and phonopy """

    if settings is None:
        settings = TaskSettings(name=None, settings=Settings())

    if "atoms" not in kwargs:
        atoms = settings.get_atoms()
    else:
        atoms = kwargs["atoms"]

    stat_sample_settings = {"atoms": atoms}

    if name not in settings:
        warn(f"Settings do not contain {name} instructions.", level=1)
    else:
        stat_sample_settings.update(settings[name])

    # Phonopy preprocess
    stat_sample_settings.update(kwargs)
    td_cells, metadata = preprocess(**stat_sample_settings)

    calc = kwargs.get(
        "calculator",
        setup_aims(
            atoms=atoms, settings=settings, custom_settings={"compute_forces": True}
        ),
    )

    # save metadata
    metadata["calculator"] = calc2dict(calc)
    to_return = list()

    return {
        "atoms_to_calculate": td_cells,
        "calculator": calc,
        "metadata": metadata,
        "workdir": name,
        **stat_sample_settings,
    }
