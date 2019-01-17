""" socket io helpers """
import numpy as np
from hilde.helpers.warnings import warn
from hilde.konstanten.einheiten import atomic_units


def get_port(calculator):
    """ return port of the calculator """

    port = None
    if "use_pimd_wrapper" in calculator.parameters:
        port = calculator.parameters["use_pimd_wrapper"][1]
        warn(f"Use SocketIO with port {port}")
    else:
        warn(f"Socketio not used with calculator {calculator.name}")

    return port


def get_stresses(atoms):
    """ Use Socket to get atomic stresses """
    atoms.calc.server.protocol.sendmsg("GETSTRESSES")
    msg = atoms.calc.server.protocol.recvmsg()
    assert msg == "STRESSREADY"
    natoms = atoms.calc.server.protocol.recv(1, np.int32)
    stresses = atoms.calc.server.protocol.recv((int(natoms), 3, 3), np.float64)
    return stresses * atomic_units.eV
