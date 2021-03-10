""" use spglib to symmetrize q points """

import collections

import numpy as np
from ase import Atoms

from vibes.helpers import Timer

from .wrapper import get_symmetry_dataset


def get_ir_reciprocal_mesh(
    q_mesh: np.ndarray,
    primitive: Atoms,
    is_time_reversal: bool = True,
    symprec: float = 1e-5,
    eps: float = 1e-9,
    debug: bool = False,
) -> collections.namedtuple:
    """Take q-mesh in fractional coords (w.r.t. primitve) and reduce by symmetry

    Args:
        q_mesh: list of q points in frac coords w.r.t. to primitive reciprocal cell
        primitive: reference structure to determine symmetry from
        is_time_reversal: If True time reversal symmetry is preserved
        eps: finite zero

    Returns:
        (map2ir, map2full): mapping to ir. grid, mapping back to orig. grid

    """
    timer = Timer(prefix="symmetry")
    # get all pure rotations:
    spg_dataset = get_symmetry_dataset(primitive, symprec=symprec)
    all_rotations = spg_dataset["rotations"]
    all_translations = spg_dataset["translations"]

    rotations = []  # list of pure rotations w/o translations
    for rot, transl in zip(all_rotations, all_translations):
        if np.linalg.norm(transl) < eps:
            rotations.append(rot)

    # prepare indices of the irreducible prototypes and map
    prototype_indeces = []
    map_to_prototypes = np.arange(len(q_mesh), dtype=int)
    rot_to_prototypes = np.zeros(len(q_mesh), dtype=int)

    # for each q-point, check if it can be mapped under rotations to a prototype
    # while it's norm is preserved
    for iq, q in enumerate(q_mesh):

        prototype_found = False
        for ip in prototype_indeces:
            p = q_mesh[ip]
            # try to match protoype to q-point
            for ir, rot in enumerate(rotations):
                dq = rot @ p - q
                diff = np.linalg.norm(dq)
                # check if vectors are identical and norm is preserved
                if diff < eps and np.allclose(np.linalg.norm(q), np.linalg.norm(p)):
                    # if we're here, ip matches to iq under rotation ir
                    # print(f"match {iq} to {ip} under rot {ir}")
                    map_to_prototypes[iq] = ip
                    rot_to_prototypes[iq] = ir
                    prototype_found = True
                    break
        if not prototype_found:
            # no map found, append prototype
            # print(f"append prototype {iq}")
            prototype_indeces.append(iq)

    full2ir, ir2full = np.unique(map_to_prototypes, return_inverse=True)

    timer(f"q-mesh reduced from {len(ir2full)} to {len(full2ir)} points.")

    return collections.namedtuple("ir_map", ("full2ir", "ir2full"))(full2ir, ir2full)
