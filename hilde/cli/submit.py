"""`hilde run` part of the CLI"""

import click
from .misc import AliasedGroup, complete_filenames

paths = complete_filenames
_prefix = "hilde.submit"
_command = lambda c: f"hilde run {c}"


def _start(settings_file, name):
    """check if settings contain [slurm] and submit"""
    from hilde.settings import Settings
    from hilde.slurm.submit import submit as _submit

    settings = Settings(settings_file=settings_file)
    if "slurm" not in settings:
        raise click.ClickException(f"[slurm] settings not found in {settings_file}")

    dct = settings["slurm"]
    dct["name"] = name

    _submit(dct, command=_command(name))


@click.command(cls=AliasedGroup)
def submit():
    """submit a hilde workflow to slurm"""


@submit.command("aims")
@click.argument("settings", default="aims.in", type=paths)
@click.pass_obj
def aims_submit(obj, settings):
    """submit one or several aims calculations from SETTINGS (default: aims.in)"""

    _start(settings, "aims")


@submit.command("phonopy")
@click.argument("settings", default="phonopy.in", type=paths)
@click.pass_obj
def phonopy_run(obj, settings):
    """submit a phonopy calculation for SETTINGS (default: phonopy.in)"""

    _start(settings, "phonopy")


@submit.command("md")
@click.argument("settings", default="md.in", type=paths)
@click.pass_obj
def md_run(obj, settings):
    """submit an MD simulation from SETTINS (default: md.in)"""

    _start(settings, "md")
