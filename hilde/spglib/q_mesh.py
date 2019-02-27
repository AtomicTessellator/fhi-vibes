""" use spglib to symmetrize q points """

import numpy as np
import spglib as spg
from hilde.structure.convert import to_spglib_cell
from hilde.helpers import Timer


def get_ir_reciprocal_mesh(q_mesh, primitive, is_time_reversal=True, symprec=1e-5):
    r""" reduce the given q_mesh by symmetry
    Input:
        q_mesh: list of q points in reduced coordinates
        primitive: reference structure to determine symmetry from
    Output:
        mapping: which q points maps to which irreducible
        reduced_q_points: the reduced set of q-points

    Example:
        https://gitlab.com/flokno/hilde/blob/devel/examples/harmonic_analysis/irreducible_q_points/ir_qpoints.ipynb

    Remarks:
        Maybe it would be nice to return the respective rotations as well?

    """

    timer = Timer()

    my_mesh = q_mesh.copy()

    n_lp = len(my_mesh)

    # find the highest q-point in each reciprocal direction
    n_qmesh = my_mesh.max(axis=0) - my_mesh.min(axis=0) + 1

    # use spglib to create a big enough mesh
    cell = to_spglib_cell(primitive)
    sym_settings = {"is_time_reversal": is_time_reversal, "symprec": symprec}

    mapping, spg_mesh = spg.get_ir_reciprocal_mesh(n_qmesh, cell, **sym_settings)

    # list the unique grid points
    ir_grid = spg_mesh[np.unique(mapping)]

    # dict translates between full and reduced grid point indices
    dct = {}
    for nn, ii in enumerate(np.unique(mapping)):
        dct[ii] = nn

    # mapping to the unique grid points
    ir_mapping = np.array([dct[ii] for ii in mapping])

    # match the q points from my list with the ones from spg
    match_list = -np.ones(n_lp, dtype=int)
    my_mapping = -np.ones(n_lp, dtype=int)
    for i1, q1 in enumerate(my_mesh):
        for i2, q2 in enumerate(spg_mesh):
            if np.linalg.norm(q1 % n_qmesh - q2 % n_qmesh) < 1e-12:
                # exchange the q point in spg list with my one
                ir_grid[ir_mapping[i2]] = q1
                match_list[i1] = i2
                my_mapping[i1] = ir_mapping[i2]

    # sanity checks
    assert -1 not in match_list, match_list
    assert -1 not in my_mapping, my_mapping
    assert len(np.unique(match_list)) == n_lp, (len(np.unique(match_list)), match_list)

    # `my_mapping` now maps my grid points to the unique grid points from spg
    # now I want to list only the grid points relevant for me
    # -> create the list of irreducible q points for my list of q points:

    my_ir_grid = ir_grid[np.unique(my_mapping)]  # % n_qmesh

    # create dict that translates between my grid points and my unique grid points
    dct = {}
    for nn, ii in enumerate(np.unique(my_mapping)):
        dct[ii] = nn

    my_ir_mapping = np.array([dct[ii] for ii in my_mapping])

    # verify that my_ir_grid indeed containts the reduced q points
    for i, q in enumerate(my_mapping):
        assert (my_ir_grid[my_ir_mapping[i]] % n_qmesh == ir_grid[q] % n_qmesh).all()

    assert len(my_ir_grid) == len(np.unique(my_ir_grid, axis=0)), my_ir_grid

    timer(
        f"number of q points reduced from {len(my_ir_mapping)} to "
        f"{len(np.unique(my_ir_mapping))}"
    )

    return my_ir_mapping, my_ir_grid
