from . import supercell as sc
import numpy as np

def find_cubic_cell(cell,
                    target_size=100,
                    deviation=0.2,
                    lower_limit=-2, upper_limit=2,
                    verbose=False):

    smatrix = sc.supercell.find_optimal_cell(cell, np.eye(3),
                                          target_size=target_size,
                                          deviation=deviation,
                                          lower_limit=lower_limit, upper_limit=upper_limit,
                                          verbose=verbose
                                          )
    return smatrix
