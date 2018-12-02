from ase.io import read

from fireworks import Workflow

from hilde import DEFAULT_CONFIG_FILE
from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.settings import Settings
from hilde.workflows.workflow_generator import get_step_fw

settings = Settings()
atoms = read(settings.geometry.file)

fw_list = []
for fw in get_step_fw(settings, atoms=atoms):
    fw_list.append(fw)

fw_dep = {}
for i,fw in enumerate(fw_list[:-1]):
    fw_dep[fw] = fw_list[i+1]

wf = Workflow(fw_list, fw_dep)
launchpad = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
launchpad.add_wf(wf)
