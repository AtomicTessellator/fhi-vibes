import numpy as np
from hilde.konstanten.numerics import eps as eps_default


def clean_matrix(matrix, eps=eps_default):
    """ clean from small values"""
    matrix = np.array(matrix)
    for ij in np.ndindex(matrix.shape):
        if abs(matrix[ij]) < eps:
            matrix[ij] = 0
    return matrix


def get_3x3_matrix(matrix, dtype=int):
    """ get a 3x3 matrix """
    if np.size(matrix) == 1:
        supercell_matrix = matrix * np.eye(3)
    elif np.size(matrix) == 3:
        supercell_matrix = np.diag(matrix)
    elif np.size(matrix) == 9:
        supercell_matrix = np.asarray(matrix).reshape((3, 3))
    else:
        raise Exception(
            "Supercell matrix must have 1, 3, 9 elements, has {}".format(
                np.size(matrix)
            )
        )

    return supercell_matrix.astype(dtype)

