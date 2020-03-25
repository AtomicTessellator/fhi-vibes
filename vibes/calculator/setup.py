"""set up ase calculators from settings"""
from ase.calculators.calculator import get_calculator_class

from vibes import keys

from .aims import setup_aims  # noqa: F401


def from_settings(settings: dict = None):
    # replace legacy keynames
    if "control" in settings:
        settings[keys.calculator] = {keys.name: "aims"}
        settings[keys.calculator][keys.parameters] = settings.pop("control")

    if "control_kpt" in settings:
        dct = {"kpoints": settings.pop("control_kpt")}
        settings[keys.calculator].update(dct)

    if "basissets" in settings:
        dct = {"basissets": settings.pop("basissets")}
        settings[keys.calculator].update(dct)

    if "socketio" in settings:
        dct = {"socketio": settings.pop("socketio")}
        settings[keys.calculator].update(dct)

    # get calculator class and create the calculator
    calc_dict = settings[keys.calculator]
    calc_name = calc_dict.pop(keys.name)

    cls = get_calculator_class(calc_name)

    return cls(**calc_dict.get(keys.parameters, {}))
