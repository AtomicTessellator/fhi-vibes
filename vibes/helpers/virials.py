"""virials helpers"""


def get_virials(atoms):
    """return virials"""

    return atoms.calc.results["virials"]


def has_virials(atoms):
    """Check if we can obtain virials with get_virials"""

    return "virials" in atoms.calc.results
