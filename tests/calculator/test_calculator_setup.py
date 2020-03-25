from pathlib import Path

import pytest

from vibes import calculator
from vibes.settings import Settings

parent = Path(__file__).parent.absolute()

files = list(parent.glob("*.in"))


@pytest.mark.parametrize("file", files)
def test_create_calculator(file):
    settings = Settings(file, config_files=None)

    calculator.setup.from_settings(settings)


@pytest.mark.parametrize("file", files)
def test_create_calculator_context(file):
    settings = Settings(file, config_files=None)

    calculator.CalculatorContext(settings)


def test_create_aims_legacy():
    settings = Settings(parent / "aims.in", config_files=None)
    c1 = calculator.setup.from_settings(settings)

    settings = Settings(parent / "aims_old.in", config_files=None)
    c2 = calculator.setup.from_settings(settings)

    assert c1.parameters == c2.parameters, (c1, c2)
