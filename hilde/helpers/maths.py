import numpy as np
from hilde.konstanten.numerics import eps as eps_default

def clean_matrix(matrix, eps=eps_default):
    """ clean from small values"""
    matrix = np.array(matrix)
    for ij in np.ndindex(matrix.shape):
        if abs(matrix[ij]) < eps:
            matrix[ij] = 0
    return matrix
