from ase.io import read

from hilde.settings import Settings
from hilde.workflows.workflow_generator import run_step

atoms = read("geometry.in")

workflow = Settings(["workflow.cfg"])
print(workflow)
for config_file in workflow.workflow.step_files.split(","):
    run_step(config_file, "hilde.cfg", atoms)
