"""set up ase calculators from settings"""
from ase.calculators.calculator import Calculator, get_calculator_class

from vibes import keys

from .aims import setup_aims  # noqa: F401


def legacy_update(settings: dict) -> dict:
    """replace legacy keynames in settings"""
    if "control" in settings:
        settings[f"{keys.calculator}.{keys.name}"] = "aims"
        settings[keys.calculator][keys.parameters] = settings.pop("control")

    if "control_kpt" in settings:
        settings[keys.calculator]["kpoints"] = settings.pop("control_kpt")

    if "basissets" in settings:
        settings[keys.calculator]["basissets"] = settings.pop("basissets")

    if "socketio" in settings:
        settings[keys.calculator]["socketio"] = settings.pop("socketio")


def from_settings(settings: dict = None) -> Calculator:
    """get calculator class and create the calculator from settings.parameters"""
    legacy_update(settings)
    calc_dict = settings[keys.calculator]
    calc_name = calc_dict.name

    cls = get_calculator_class(calc_name)

    return cls(**calc_dict.get(keys.parameters, {}))
