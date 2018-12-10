""" Provide highlevel access to functions that can be used with fireworks PyTasks """

from .utility_tasks import add_phonon_to_db
from .utility_tasks import add_result_to_spec
from .utility_tasks import check_convergence
from .utility_tasks import get_relaxed_structure
from .utility_tasks import mod_calc
from .utility_tasks import transfer_spec
from .utility_tasks import update_calc_in_db

from .general_py_task import atoms_calculate_task
from .general_py_task import general_function_task

from hilde.phonopy.workflow import preprocess as ph_preprocess
from hilde.phono3py.workflow import preprocess as ph3_preprocess

def preprocess_ph_ph3(atoms, calc, phonopy_settings=None, phono3py_settings=None):
    if phonopy_settings:
        phonopy_preprocess = ph_preprocess(atoms, calc, **phonopy_settings)
    else:
        phonopy_preprocess = None

    if phono3py_settings:
        phono3py_preprocess = ph3_preprocess(atoms, calc, **phono3py_settings)
    else:
        phono3py_preprocess = None
    return phonopy_preprocess, phono3py_preprocess

