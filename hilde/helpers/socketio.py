""" socket io helpers """
from hilde.helpers.warnings import warn


def get_port(calculator):
    """ return port of the calculator """

    port = None
    if "use_pimd_wrapper" in calculator.parameters:
        port = calculator.parameters["use_pimd_wrapper"][1]
        warn(f"Use SocketIO with port {port}")
    else:
        warn(f"Socketio not used with calculator {calculator.name}")

    return port
