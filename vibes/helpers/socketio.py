""" socket io helpers """
import numpy as np
from ase import units

from vibes.helpers import talk
from vibes.helpers.warnings import warn


_prefix = "socketio"


def get_port(calculator, prefix=_prefix):
    """return port of the calculator

    Args:
        calculator: calculator to get the port of
        prefix: prefix for messages from this function
    Returns:
        port: the port for socketio
    """

    port = None

    if not hasattr(calculator, "parameters"):
        warn(f"{prefix} No parameters found in calculator {calculator.name}.", level=1)
        return port

    if "use_pimd_wrapper" in calculator.parameters:
        port = calculator.parameters["use_pimd_wrapper"][1]
        if "UNIX:" in calculator.parameters["use_pimd_wrapper"][0]:
            talk(f"Use SocketIO with unixsocket and port {port}", prefix=prefix)
        else:
            talk(f"Use SocketIO with port {port}", prefix=prefix)
    else:
        talk(f"Socketio not used with calculator {calculator.name}", prefix=prefix)

    return port


def get_unixsocket(calculator, prefix=_prefix):
    """return the unixsocket of the calculator

    Args:
        calculator: calculator to get the unixsocket of
        prefix: prefix for messages from this function
    Returns:
        unixsocket: the unixsocket used for socketio
    """

    unixsocket = None

    if not hasattr(calculator, "parameters"):
        warn(f"{prefix} No parameters found in calculator {calculator.name}.", level=1)
        return unixsocket

    if "use_pimd_wrapper" in calculator.parameters:
        host = calculator.parameters["use_pimd_wrapper"][0]
        if host[:5] == "UNIX:":
            unixsocket = host[5:]
            talk(f"Use SocketIO with unixsocket {unixsocket}", prefix=prefix)
        else:
            talk(f"Use SocketIO with localhost", prefix=prefix)
    else:
        talk(f"SocketIO not used with calculator {calculator.name}", prefix=prefix)

    return unixsocket


def get_stresses(atoms):
    """Use Socket to get atomic stresses

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        atoms of the calculation to get teh stress of

    Returns
    -------
    The atomic stress

    Raises
    ------
    AssertionError
        If STRESSREADY is not sent in response to GETSTRESSES message
    """
    if "socketio" not in atoms.calc.name.lower():
        return atoms.calc.get_stresses()
    atoms.calc.server.protocol.sendmsg("GETSTRESSES")
    msg = atoms.calc.server.protocol.recvmsg()
    assert msg == "STRESSREADY"
    natoms = atoms.calc.server.protocol.recv(1, np.int32)
    stresses = atoms.calc.server.protocol.recv((int(natoms), 3, 3), np.float64)
    return stresses * units.Hartree


def socket_stress_off(calc):
    """Turn stresses computation off via socket

    Parameters
    ----------
    calc: ase.calculators.calulator.Calculator
        calculator to turn off stress computation for
    """
    if "socketio" in calc.name.lower():
        calc.server.protocol.sendmsg("STRESSES_OFF")
    else:
        talk(f"Calculator {calc.name} is not a socket calculator.")
        calc.parameters["compute_heat_flux"] = False
        calc.parameters["compute_analytical_stress"] = True


def socket_stress_on(calc):
    """ Turn stresses computation on via socket

    Parameters
    ----------
    calc: ase.calculators.calulator.Calculator
        calculator to turn on stress computation for
    """
    if "socketio" in calc.name.lower():
        calc.server.protocol.sendmsg("STRESSES_ON")
    else:
        talk(f"Calculator {calc.name} is not a socket calculator.")
        calc.parameters["compute_heat_flux"] = True
        del calc.parameters["compute_analytical_stress"]
