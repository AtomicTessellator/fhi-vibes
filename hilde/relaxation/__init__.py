""" Relaxation.
 * Optimizers from ASE
 * SocketIO
 * Yaml Trajectory """

from hilde.trajectory import input2dict
from .bfgs import relax as bfgs_relax


def metadata2dict(atoms, calc, opt):
    """ convert metadata information to plain dict """
    opt_dict = opt.todict()

    return {"geometry_optimization": opt_dict, **input2dict(atoms, calc)}
