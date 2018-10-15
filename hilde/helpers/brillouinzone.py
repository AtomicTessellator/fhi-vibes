"""
Utility functions for working with Brillouin zones
"""

from ase.dft.kpoints import get_cellinfo, special_paths, bandpath

def get_paths(atoms):
    """ Get recommended path connencting high symmetry points in the BZ.

    Args:
        atoms (Atoms): atoms object

    Returns:
        list: Recommended special points als list
        >>> ['GXL', 'KU']

    """
    cellinfo = get_cellinfo(atoms.cell)
    paths = special_paths[cellinfo.lattice].split(',')
    return paths


def get_bands(atoms, paths):
    """ Get the recommended BZ path(s) for atoms """
    bands = []
    for path in paths:
        for ii, _ in enumerate(path[:-1]):
            bands.append(
                bandpath(path[ii : ii+2], atoms.cell)[0])
    return bands


def get_labels(paths):
    """ Get the labels for a given path for printing them with latex """
    if len(paths) == 1:
        labels = [*paths]
    else:
        labels = [*'|'.join(paths)]
        for ii, l in enumerate(labels):
            if l == '|':
                labels[ii] = f"{labels[ii-1]}|{labels[ii+1]}"
                labels[ii-1], labels[ii+1] = '', ''
            if l == 'G':
                labels[ii] = '\\Gamma'
        labels = [l for l in labels if l]

    latexify = lambda sym: "$\\mathrm{\\mathsf{" + str(sym)  + "}}$"
    return [latexify(sym) for sym in labels]


def get_bands_and_labels(atoms, paths=None):
    """ Combine get_bands() and get_labels() """
    if paths is None:
        paths = get_paths(atoms)

    return get_bands(atoms, paths), get_labels(paths)