""" Utility functions for working with Brillouin zones """

import collections

import numpy as np
from ase import Atoms
from ase.dft import kpoints

from vibes.helpers import Timer
from vibes.helpers.latex import latexify_labels
from vibes.spglib import get_symmetry_dataset


def get_paths(atoms: np.ndarray) -> list:
    """nothing but atoms.get_bravais_lattice().special_path.split(',')"""
    return atoms.cell.get_bravais_lattice().special_path.split(",")


def get_bands(atoms: Atoms, paths: np.ndarray = None, npoints: int = 50) -> list:
    """Get the recommended BZ path(s) for atoms

    Args:
        atoms: The structure to get the recommended high-symmetry point path
        paths: Paths connecting high-symmetry points
        npoints: Number of points for each band

    Returns:
        bands: The recommended BZ path(s) for atoms

    """
    if paths is None:
        paths = get_paths(atoms)
    bands = []
    for path in paths:
        points = kpoints.parse_path_string(path)[0]  # [:-1]
        ps = [points.pop(0)]
        for _, p in enumerate(points):
            ps.append(p)
            bands.append(atoms.cell.bandpath("".join(ps), npoints=npoints).kpts)
            ps.pop(0)
    return bands


def get_labels(paths: list, latex: bool = False) -> list:
    """Get the labels for a given path for printing them with latex

    Args:
        paths: Paths connecting high-symmetry points
        latex: If True convert labels to Latex format

    Returns:
        labels: The labels for the high-symmetry path

    """
    if len(paths) == 1:
        labels = kpoints.parse_path_string(paths[0])[0]
        labels.append("|")
    else:
        labels = []
        for path in paths:
            points = kpoints.parse_path_string(path)[0]
            labels.extend(points)
            labels.append("|")

    # discard last |
    labels = labels[:-1]

    for ii, ll in enumerate(labels):
        if ll == "|":
            labels[ii] = f"{labels[ii-1]}|{labels[ii+1]}"
            labels[ii - 1], labels[ii + 1] = "", ""

    labels = [ll for ll in labels if ll]

    if latex:
        return latexify_labels(labels)

    return labels


def get_bands_and_labels(
    atoms: Atoms, paths: list = None, npoints: int = 50, latex: bool = False
) -> tuple:
    """Combine get_bands() and get_labels()

    Args:
        atoms: The structure to get the recommended high-symmetry point path
        paths: Paths connecting high-symmetry points
        npoints: Number of points for each band
        latex: If True convert labels to Latex format

    Returns:
        bands: The recommended BZ path(s) for atoms
        labels: The labels for the high-symmetry path

    """
    if paths is None:
        paths = get_paths(atoms)

    bands = get_bands(atoms, paths, npoints=npoints)
    labels = get_labels(paths, latex=latex)

    return collections.namedtuple("bands_and_labels", ("bands", "labels"))(
        bands, labels
    )


def get_ir_grid(
    q_points: np.ndarray,
    primitive: Atoms,
    is_time_reversal: bool = True,
    symprec: float = 1e-5,
    eps: float = 1e-9,
) -> collections.namedtuple:
    """Take q-points in fractional coords (w.r.t. primitve) and reduce by symmetry

    Args:
        q_points: list of q points in frac coords w.r.t. to primitive reciprocal cell
        primitive: reference structure to determine symmetry from
        is_time_reversal: If True time reversal symmetry is preserved
        eps: finite zero

    Returns:
        q_grid
          .points: q-points in the given grid (fractional w.r.t to primitive)
          .points_cartesian: q-points in cart. coords
          .ir_indices: indices of the irreducible prototypes
          .ir_points: the irreducible q-points
          .ir_weights: the corresponding weights
          .map2ir_points: map from points to correspnding ir. point (ir_points)
          .map2ir_indices: map from points to corresponding ir. index (ir_indices)
          .spg_data: spg_dataset including rotations in frac. and cart. coords
          .symop2ir: index of symmetry operation that transforms point to ir_point

    """
    timer = Timer(message="reduce q-grid", prefix="symmetry")
    # get all pure rotations:
    spg_dataset = get_symmetry_dataset(primitive, index_maps=True, symprec=symprec)
    rotations = spg_dataset.rotations  # _cartesian
    # translations = spg_dataset.translations
    dummy_indices = np.arange(len(rotations), dtype=int)

    # prepare indices of the irreducible prototypes and map
    ir_indices = []
    map2ir_indices = np.arange(len(q_points), dtype=int)
    symop2ir = np.zeros(len(q_points), dtype=int)

    # for each q-point, check if it can be mapped under rotations to a ir. prototype
    for iq, q in enumerate(q_points):

        prototype_found = False
        for ir in ir_indices:
            p = q_points[ir]  # pick a ir. prototype

            # check if norm is preserved. REM: this is only true in Cart. coords
            # if not abs(q.dot(q) - p.dot(p)) < eps:
            #     continue

            # try to match protoype to q-point
            # use vectorization
            # create list of all rotations applied to q
            # multiply from the right, we look for q . S = p
            # REM: This would only be equivalent to q = S . P in Cart. coords
            rotated_qs = (q[None, :, None] * rotations).sum(axis=-2)
            # compute the differences to q point
            diffs = np.square(rotated_qs - p[None, :]).sum(axis=-1)
            matched_rotations = dummy_indices[diffs < eps]
            if len(matched_rotations) > 0:  # match found
                ig = matched_rotations[0]  # index of symmetry op.
                symop2ir[iq] = ig
                map2ir_indices[iq] = ir
                prototype_found = True
                # print(f"{iq:2d} -> {ir:2d} under op. {ig:2d}")
                # print(q, rotated_qs[ig])
                break

        if not prototype_found:
            # no map found, append prototype
            # print(f"append prototype {iq}")
            ir_indices.append(iq)

    # get map to ir_points and weights
    ir_indices_check, map2ir_points, ir_weigths = np.unique(
        map2ir_indices, return_inverse=True, return_counts=True
    )
    # sanity check
    assert np.allclose(ir_indices, ir_indices_check)

    # prepare and return results including cartesian rotations
    q_points = q_points.copy()
    q_points_cart = primitive.cell.reciprocal().cartesian_positions(q_points)

    ir_points = q_points[ir_indices]

    timer(f"q-points reduced from {len(q_points)} to {len(ir_points)} points.")

    data = {
        "points": q_points,
        "points_cartesian": q_points_cart,
        "ir_indices": ir_indices,
        "ir_points": ir_points,
        "ir_weights": ir_weigths,
        "map2ir_points": map2ir_points,
        "map2ir_indices": map2ir_indices,
        "spg_data": spg_dataset,
        "symop2ir": symop2ir,
    }

    QGrid = collections.namedtuple("q_grid", data.keys())

    return QGrid(**data)
