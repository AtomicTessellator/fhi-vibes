from ase.io import read

from hilde.settings import Settings
from hilde.fireworks.workflows.workflow_generator import generate_workflow

settings = Settings()
print(settings)
atoms = read(settings.geometry.file)
steps = [settings]
fw_settings = {"to_launchpad": True}
fw_settings["name"] = (
    atoms.symbols.get_chemical_formula()
)
generate_workflow(steps, atoms=atoms, fw_settings=fw_settings)
