import numpy as np
from hilde.helpers.supercell import find_cubic_cell
from pathlib import Path

def get_smatrix(at, n_target=64):
    """
    Return the supercell matrix for atoms with target size
    Args:
        at: pAtoms
            The primitive cell
        n_target: int
            The target number of atoms in the supercell
    Returns:
        The supercell matrix that produces the cubic super cell
    """
    target_size = n_target / len(at)
    return find_cubic_cell(cell=at.cell, target_size=target_size)

def setup_workdir(at):
    """
    Set up a working directory for a single point calculation
    Args:
        at: pAtoms
            The primitive cell of the calculation
    Returns:
        The Path object to the work directory
    """
    vol = at.get_volume()
    wd = Path('./{}/{:.3f}_EMT'.format(at.sysname, vol)).absolute()
    wd.mkdir(parents=True, exist_ok=True)
    return wd

def reshape_fc(fc_arr):
	'''
	Reshapes force constant array from a linear 1D array to a 4D array
	Args:
		fc_arr: np.ndarray(shape=(3*3*nAtoms*nAtoms,))
			Linear force constant array
	Returns:
		fc_arr: np.ndarray(shape=(nAtoms, nAtoms, 3, 3))
			properly formatted force constant array
	'''
	return fc_arr.reshape(int(np.sqrt(len(fc_arr.flatten())) / 3),
		                  int(np.sqrt(len(fc_arr.flatten())) / 3),
		                  3,
		                  3)
