"""helpers for click"""

from pathlib import Path
import click

complete_filenames = click.Path(exists=True)


class AliasedGroup(click.Group):
    """modified groupd for shorthand arguments

    see https://click.palletsprojects.com/en/7.x/advanced/#command-aliases

    """

    def get_command(self, ctx, cmd_name):
        """geturn command based on prefix"""

        rv = click.Group.get_command(self, ctx, cmd_name)
        if rv is not None:
            return rv

        matches = [x for x in self.list_commands(ctx) if x.startswith(cmd_name)]

        if not matches:
            return None
        if len(matches) == 1:
            return click.Group.get_command(self, ctx, matches[0])

        ctx.fail("Too many matches: %s" % ", ".join(sorted(matches)))


def check_path(filename):
    """check if path exists"""
    if not Path(filename).exists():
        msg = f"\n  Current workdir is {Path().cwd()}"
        raise click.FileError(str(filename), hint=msg)
