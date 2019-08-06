#!/usr/bin/env python3
"""Register the hilde applications"""
import os
from balsam.core.models import ApplicationDefinition as App
from hilde.settings import Settings

settings = Settings()

if "aimsbin" in settings.machine:
    aimsbin = settings.machine.aimsbin
else:
    aimsbin = settings.machine.aims_command

balsam_path = os.path.abspath(__file__)
bp_list = balsam_path.split("/")
balsam_path = "/".join(bp_list[:-1])
print(balsam_path, bp_list, __file__)


def register():
    os.chmod(f"{balsam_path}/pre/aims.py", 0o755)
    os.chmod(f"{balsam_path}/post/aims.py", 0o755)
    os.chmod(f"{balsam_path}/post/phonopy_setup.py", 0o755)
    os.chmod(f"{balsam_path}/post/phonopy_post.py", 0o755)

    if not App.objects.filter(name="run-aims").exists():
        aims_run = App(
            name="run-aims",
            executable=f"{aimsbin} > aims.out",
            preprocess=f"{balsam_path}/pre/aims.py",
            postprocess=f"{balsam_path}/post/aims.py",
            envscript=f"{balsam_path}/envs/theta_env.sh",
        )
        aims_run.save()

    if not App.objects.filter(name="run-aims-setup-phonopy").exists():
        aims_run = App(
            name="run-aims-setup-phonopy",
            executable=f"{aimsbin} > aims.out",
            preprocess=f"{balsam_path}/pre/aims.py",
            postprocess=f"{balsam_path}/post/phonopy_setup.py",
            envscript=f"{balsam_path}/envs/theta_env.sh",
        )
        aims_run.save()

    if not App.objects.filter(name="run-aims-post-phonopy").exists():
        aims_run = App(
            name="run-aims-post-phonopy",
            executable=f"{aimsbin} > aims.out",
            preprocess=f"{balsam_path}/pre/aims.py",
            postprocess=f"{balsam_path}/post/phonopy_post.py",
            envscript=f"{balsam_path}/envs/theta_env.sh",
        )
        aims_run.save()
