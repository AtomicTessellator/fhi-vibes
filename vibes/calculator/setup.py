"""set up ase calculators from settings"""
from ase.calculators.calculator import Calculator, get_calculator_class

from vibes import keys


def from_settings(settings: dict = None) -> Calculator:
    """get calculator class and create the calculator from settings.parameters"""
    calc_dict = settings[keys.calculator]
    calc_name = calc_dict.name.lower()

    cls = get_calculator_class(calc_name)

    return cls(**calc_dict.get(keys.parameters, {}))
