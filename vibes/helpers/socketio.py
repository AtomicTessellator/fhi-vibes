""" socket io helpers """
import numpy as np
from ase import units

from vibes.helpers import talk
from vibes.helpers.warnings import warn

from . import stresses as stresses_helper

_prefix = "socketio"


def get_socket_info(calculator, prefix=_prefix):
    """return port of the calculator

    Args:
        calculator: calculator to get the port of
        prefix: prefix for messages from this function
    Returns:
        host: The host for socketio
        port: the port for socketio
        unixsocket: get_unixsocket
    """

    port = None
    unixsocket = None
    if not hasattr(calculator, "parameters"):
        warn(f"{prefix} No parameters found in calculator {calculator.name}.", level=1)
        return port

    if "use_pimd_wrapper" in calculator.parameters:
        port = calculator.parameters["use_pimd_wrapper"][1]
        host = calculator.parameters["use_pimd_wrapper"][0]
        if "UNIX:" in host:
            unixsocket = calculator.parameters["use_pimd_wrapper"][0]
            talk(f"Use SocketIO with unixsocket file {unixsocket}", prefix=prefix)
        else:
            talk(f"Use SocketIO with host {host} and port {port}", prefix=prefix)
    else:
        talk(f"Socketio not used with calculator {calculator.name}", prefix=prefix)

    return port, unixsocket


def get_stresses(atoms):
    """Use Socket to get intensive atomic stresses in eV/AA^3 in Nx3x3 shape.
    Raw stresses are supposed to be extensive and the volume is divided out.

    """
    if "socketio" not in atoms.calc.name.lower():
        return stresses_helper.get_stresses(atoms)
    # assume these are extensive stresses
    atoms.calc.server.protocol.sendmsg("GETSTRESSES")
    msg = atoms.calc.server.protocol.recvmsg()
    assert msg == "STRESSREADY"
    natoms = atoms.calc.server.protocol.recv(1, np.int32)
    stresses = atoms.calc.server.protocol.recv((int(natoms), 3, 3), np.float64)
    return stresses * units.Hartree / atoms.get_volume()


def socket_stress_off(calculator):
    """Turn stresses computation off via socket

    Parameters
    ----------
    calculator: ase.calculators.calulator.Calculator
        calculator to turn off stress computation for
    """
    if "socketio" in calculator.name.lower():
        calculator.server.protocol.sendmsg("STRESSES_OFF")
    else:
        talk(f"Calculator {calculator.name} is not a socket calculator.")
        calculator.parameters["compute_heat_flux"] = False
        calculator.parameters["compute_analytical_stress"] = True


def socket_stress_on(calculator):
    """ Turn stresses computation on via socket

    Parameters
    ----------
    calculator: ase.calculators.calulator.Calculator
        calculator to turn on stress computation for
    """
    if "socketio" in calculator.name.lower():
        calculator.server.protocol.sendmsg("STRESSES_ON")
    else:
        talk(f"Calculator {calculator.name} is not a socket calculator.")
        if "aims" in calculator.name.lower():
            talk(f"Switch on `compute_heat_flux` for {calculator.name}")
            calculator.parameters["compute_heat_flux"] = True
            del calculator.parameters["compute_analytical_stress"]
