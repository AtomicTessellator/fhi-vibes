from ase.io import read
from hilde.helpers.k_grid import update_k_grid
from hilde.phonopy.workflow import phonopy
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.templates.aims import setup_aims

atoms = read("geometry.in")

workflow = Settings(["workflow.cfg"])

if "relaxation" in workflow:
    relax_settings = Settings(["hilde.cfg"])
    for key, val in  workflow.relaxation_control.items():
        relax_settings.control[key] = val
    calc = setup_aims(settings=relax_settings)
    update_k_grid(atoms, calc, workflow.relaxation_control_kpt.density)
    relax(
        atoms,
        calc,
        socketio_port=relax_settings.socketio.port,
        **workflow.relaxation
    )
print("to_phonons")
if "phonopy" in workflow:
    phonopy_settings = Settings(["hilde.cfg"])
    for key, val in workflow.phonopy_control.items():
        phonopy_settings.control[key] = val
    calc = setup_aims(settings=phonopy_settings)
    update_k_grid(atoms, calc, workflow.phonopy_control_kpt.density)
    phonopy(atoms, calc, socketio_port=phonopy_settings.socketio.port, **workflow.phonopy)
