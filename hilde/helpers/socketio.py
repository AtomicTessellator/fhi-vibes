""" socket io helpers """
import numpy as np
from ase import units
from hilde.helpers import talk
from hilde.helpers.warnings import warn


def get_port(calculator, prefix="socketio"):
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
        talk(f"Use SocketIO with port {port}", prefix=prefix)
    else:
        talk(f"Socketio not used with calculator {calculator.name}", prefix=prefix)

    return port


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
