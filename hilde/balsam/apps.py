#!/usr/bin/env python3
"""Register the hilde applications"""
import os
from balsam.core.models import ApplicationDefinition as App

balsam_path = os.path.abspath(__file__)
bp_list = balsam_path.split("/")
balsam_path = "/".join(bp_list[:-1])
print(balsam_path, bp_list, __file__)
def register():
    os.chmod(f'{balsam_path}/pre/aims.py', 0o755)
    os.chmod(f'{balsam_path}/pre/phonons.py', 0o755)
    os.chmod(f'{balsam_path}/post/aims.py', 0o755)
    os.chmod(f'{balsam_path}/post/phonons.py', 0o755)
    if not App.objects.filter(name="hilde-run-aims").exists():
        aims_run = App(
            name = 'hilde-run-aims',
            executable = 'hilde run aims',
            preprocess = f'{balsam_path}/pre/aims.py',
            postprocess = f'{balsam_path}/post/aims.py',
            envscript = f'{balsam_path}/envs/theta_env.sh',
        )
        aims_run.save()

    if not App.objects.filter(name="hilde-run-phonopy").exists():
        phonopy_run = App(
            name = 'hilde-run-phonopy',
            executable = 'hilde run phonopy',
            preprocess = f'{balsam_path}/pre/phonons.py',
            postprocess = f'{balsam_path}/post/phonons.py',
            envscript = f'{balsam_path}/envs/theta_env.sh',
        )
        phonopy_run.save()

    if not App.objects.filter(name="hilde-run-md").exists():
        md_run = App(
            name = 'hilde-run-md',
            executable = 'hilde run md',
            preprocess = f'{balsam_path}/pre/md.py',
            postprocess = f'{balsam_path}/post/md.py',
            envscript = f'{balsam_path}/envs/theta_env.sh',
        )
        md_run.save()
