"""CLI for hilde with click"""

import shutil
from pathlib import Path
import importlib.resources as pkg_resources

import click
import click_completion

from hilde import __version__ as hilde_version
from hilde.settings import Configuration
from hilde.helpers.utils import bold
from hilde._defaults import DEFAULT_CONFIG_FILE
from .cli_tracker import CliTracker
from . import info, template, run, utils, output, aigk
from .misc import AliasedGroup, check_path

click_completion.init()


@click.command(cls=AliasedGroup)
@click.version_option(hilde_version, "-V", "--version")
@click.option("-v", "--verbose", is_flag=True, hidden=True)
@click.option("--silent", is_flag=True, help="set verbosity level to 0")
@click.pass_context
def cli(ctx, verbose, silent):
    """hilde: lattice dynamics with pyhon"""
    if verbose:
        verbosity = 2
    elif silent:
        verbosity = 0
    else:
        verbosity = 1

    ctx.obj = CliTracker(verbose=verbosity)
    ctx.help_option_names = ["-h", "--help"]

    if verbosity > 1:
        click.echo(f"Welcome to hilde!\n")


cli.add_command(info.info)
cli.add_command(template.template)
cli.add_command(run.run)
cli.add_command(utils.utils)
cli.add_command(output.output)
cli.add_command(aigk.aiGK)

try:
    from hilde.fireworks.cli import fireworks

    cli.add_command(fireworks.fireworks)
except ImportError:
    pass


@cli.command("status", hidden=True)
@click.option("--verbose", is_flag=True)
@click.pass_obj
def hilde_status(obj, verbose):
    """check if everything is set up"""

    configfile = DEFAULT_CONFIG_FILE

    check_path(configfile)

    config = Configuration(config_file=configfile)

    if verbose:
        click.echo(f"This is the configuration found in {configfile}")
        click.echo(config)

    try:
        assert "machine" in config
        check_path(config.machine.basissetloc)
        check_path(shutil.which(config.machine.aims_command))
    except (AssertionError, click.FileError):
        msg = f"Check your configuration file in {configfile}"
        click.ClickException(msg)

    click.secho("It seems we are good to go!", bold=True)
