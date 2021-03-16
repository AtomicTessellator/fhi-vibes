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
    return atoms.cell.get_bravais_lattice().special_path().split(",")


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

    return bands, labels
