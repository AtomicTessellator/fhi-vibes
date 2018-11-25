module_name = __name__

import numpy as np


def mutate_kgrid(atoms):
    new_atoms = dict(atoms)
    new_atoms["calculator_parameters"]["k_grid"] += np.ones(3)
    return new_atoms, new_atoms["calculator_parameters"]["k_grid"]


mutate_kgrid.name = f"{module_name}.{mutate_kgrid.__name__}"
