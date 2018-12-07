from hilde import Settings, setup_aims
from hilde.phonopy.workflow import run as phonopy_workflow
from hilde.fireworks import Launchpad

settings = Settings()

atoms, calc = setup_aims(settings=settings)

launchpad = LaunchPad()
launchpad.add_wf(atoms=atoms, calc=calc, workflow=phonopy_workflow)
