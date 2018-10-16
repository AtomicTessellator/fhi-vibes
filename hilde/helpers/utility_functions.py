import numpy as np

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
		                  3
		                 )