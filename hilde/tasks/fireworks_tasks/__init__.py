""" Provide highlevel access to functions that can be used with fireworks PyTasks """
from .phonopy_tasks import analyze_phonopy, initialize_phonopy
from .elec_struct_tasks import calculate_sp, calculate_multiple_sp

analyze_phonopy_py_task = analyze_phonopy.__module__ + "." + analyze_phonopy.__name__
calculate_py_task = calculate_sp.__module__ + "." + calculate_sp.__name__
initialize_phonopy_py_task = initialize_phonopy.__module__ + "." + initialize_phonopy.__name__
calculate_mult_sp_py_task = calculate_multiple_sp.__module__ + "." + calculate_multiple_sp.__name__
