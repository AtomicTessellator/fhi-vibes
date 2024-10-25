"""use spglib to symmetrize q points"""
import collections

import numpy as np
import spglib as spg
from ase.atoms import Atoms

from vibes.structure.convert import to_spglib_cell


def get_ir_reciprocal_mesh(
    mesh: np.ndarray,
    atoms: Atoms,
    monkhorst: bool = True,
    symprec: float = 1e-5,
    eps: float = 1e-9,
) -> collections.namedtuple:
    """
    wrapper for spglib.get_ir_reciprocal_mesh

    Args:
    ----
        mesh: array specifying number of points per axis
        atoms: structure that determines symmetry
        monkhorst: return monkhorst-pack-style grid (gamma incl. when odd grid number)

    Returns:
    -------
        (points, mapping): frac. grid points in interval [-0.5, 0.5), mapping to ir.

    """
    mesh = np.asarray(mesh)
    spg_cell = to_spglib_cell(atoms)

    if monkhorst:  # apply shift
        is_shift = -(np.asarray(mesh) % 2) + 1
    else:
        is_shift = np.zeros(3)

    mapping, grid = spg.get_ir_reciprocal_mesh(
        mesh, cell=spg_cell, is_shift=is_shift
    )

    points = grid.astype(float) / mesh  # to frac. coords
    points += 0.5 * is_shift / mesh  # restore shift

    # map to [-0.5, 0.5)
    points = ((points + 0.5 + eps) % 1 - 0.5 - eps).round(decimals=14)

    data = {"points": points, "mapping": mapping}

    IrReciprocalMesh = collections.namedtuple("ir_reciprocal_mesh", data.keys())

    return IrReciprocalMesh(**data)
