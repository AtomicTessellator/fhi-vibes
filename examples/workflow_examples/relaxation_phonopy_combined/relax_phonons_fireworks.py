from ase.io import read

from fireworks import Workflow

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.settings import Settings
from hilde.workflows.workflow_generator import get_step_fw

atoms = read("geometry.in")

workflow = Settings(["workflow.cfg"])
print(workflow)
fw_list = []
for config_file in workflow.workflow.step_files.split(","):
    for fw in get_step_fw(config_file, "hilde.cfg", atoms):
        fw_list.append(fw)

fw_dep = {}
for i,fw in enumerate(fw_list[:-1]):
    fw_dep[fw] = fw_list[i+1]
wf = Workflow(fw_list, fw_dep)
lp = LaunchPadHilde(port=27019)
lp.add_wf(wf)
