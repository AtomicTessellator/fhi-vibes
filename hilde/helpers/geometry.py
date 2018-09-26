import numpy as np
from hilde.konstanten.maths import perfect_fill, vol_sphere
# from collections import namedtuple

def inscribed_sphere_in_box(cell):

    # the normals of the faces of the box
    na = np.cross(cell[1, :], cell[2, :])
    nb = np.cross(cell[2, :], cell[0, :])
    nc = np.cross(cell[0, :], cell[1, :])
    na /= np.linalg.norm(na)
    nb /= np.linalg.norm(nb)
    nc /= np.linalg.norm(nc)
    # distances between opposing planes
    rr = 1.e10
    rr = min(rr, abs(na @ cell[0, :]))
    rr = min(rr, abs(nb @ cell[1, :]))
    rr = min(rr, abs(nc @ cell[2, :]))
    rr *= 0.5
    return rr

def get_cubicness(cell):
    """
    Purpose:
    Quantify 'how cubic' a given lattice or cell is by comparing the largest
    sphere that fits into the cell to a sphere that fits into a cubic cell
    of similar size
    Args:
        cell: Lattice of the cell
    Returns:
        ratio of radii of the two spheres:
          * 1 means perfectly cubic,
          * ratio**3 compared volumes of the two spheres
    """

    # perfect radius: 1/2 * width of the cube
    radius_perfect = np.linalg.det(cell)**(1/3) * .5
    radius_actual  = inscribed_sphere_in_box(cell)

    # volume = vol_sphere * inscribed_sphere_in_box(cell)**3 / np.linalg.det(cell)
    # Fill = namedtuple('Fill', ['volume', 'radius'])
    # fill = Fill(volume=volume, radius=radius)

    return radius_actual / radius_perfect

