""" socket io helpers """
from hilde.helpers.warnings import warn


def get_port(calculator):
    """ return port of the calculator """

    port = None
    if "use_pimd_wrapper" in calculator.parameters:
        port = calculator.parameters["use_pimd_wrapper"][1]
    else:
        warn(f"{calculator.name} is not supported.")

    return port
