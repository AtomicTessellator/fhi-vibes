import numpy as np

def reshapeFC(forceConstantArr):
	'''
	Reshapes force constant array from a linear 1D array to a 4D array
	Args:
		forceConstantArr: np.ndarray(shape=(3*3*nAtoms*nAtoms,))
			Linear force constant array
	Returns:
		forceCosntantArr: np.ndarray(shape=(nAtoms, nAtoms, 3, 3))
			properly formatted force constant array
	'''
	return forceConstantArr.reshape(int(np.sqrt(len(self.force_constants))/3), int(np.sqrt(len(self.force_constants))/3),3,3)