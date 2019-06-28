""" Copied from hiphive:
    File: https://gitlab.com/materials-modeling/hiphive/blob/1382b48a6aae530c4cc4e3dd29a0342f3071659a/hiphive/md_tools/spectral_energy_density.py
    License: https://gitlab.com/materials-modeling/hiphive/blob/master/LICENSE
"""

import numpy as np


def compute_sed(traj, ideal, prim, k_points):
    """ Computes spectral energy density for a trajectory.

    Parameters
    ----------
    traj : list of ase.atoms.Atoms
        trajectory with atoms objects with velocities
    ideal : ase.atoms.Atoms
        ideal atoms object
    prim : ase.atoms.Atoms
        compatible primitive cell. Must be aligned correctly
    k_points : list
        list of k points in cart coord (2pi must be included)

    Returns
    -------
    density: np.ndarray(dtype=float, shape=(len(k_points), velocities.shape[2]))
        The spectral density of the trajectory
    """

    velocities = []
    for atom in traj:
        velocities.append(atom.get_velocities())
    velocities = np.array(velocities)
    velocities = velocities.transpose(1, 2, 0).copy()
    velocities = np.fft.fft(velocities, axis=2)

    masses = prim.get_masses()
    indices, offsets = _index_offset(ideal, prim)

    # pos = np.dot(k_points, np.dot(offsets, prim.cell).T)
    pos = np.dot(k_points, ideal.positions.T)

    exppos = np.exp(1.0j * pos)
    density = np.zeros((len(k_points), velocities.shape[2]))
    for alpha in range(3):
        for b in range(len(masses)):
            tmp = np.zeros(density.shape, dtype=np.complex)
            for i in range(len(indices)):
                index = indices[i]
                if index != b:
                    continue
                # sum similar atoms with the same phase
                tmp += np.outer(exppos[:, i], velocities[i, alpha])

            density += masses[b] * np.abs(tmp) ** 2

    return density


def _index_offset(atoms, prim, atol=1e-3, rtol=0.0):
    """Computes the index_offset

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The atoms object
    prim: ase.atoms.Atoms
        The primitive cell
    atol: float
        The absolute tolerance
    rtol: float
        The relative tolerance

    Returns
    -------
    index: np.ndarray(int)
        The index array
    offset: np.ndarray(int)
        The offset

    Raises
    ------
    ValueError
        If prim is not compatible with atoms
    AssertionError
        the type of offset is not int
    """
    index, offset = [], []
    for pos in atoms.positions:
        spos = np.linalg.solve(prim.cell.T, pos)
        for i, spos2 in enumerate(prim.get_scaled_positions()):
            off = spos - spos2
            off_round = np.round(off)
            if not np.allclose(off, off_round, atol=atol, rtol=rtol):
                continue
            index.append(i)
            off = off_round.astype(int)
            assert np.allclose(off, off_round)
            offset.append(off)
            break
        else:
            raise ValueError("prim not compatible with atoms")

    index, offset = np.array(index), np.array(offset)
    return index, offset
