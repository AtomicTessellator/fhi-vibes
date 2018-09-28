from ase.calculators.calculator import kptdensity2monkhorstpack


def d2k(*args, **kwargs):
    """ Wrapper for kptdensity2monkhorstpack to return list """
    return list(kptdensity2monkhorstpack(*args, **kwargs))