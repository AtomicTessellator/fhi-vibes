""" Provide highlevel access to calculator.calculate """

from subprocess import run
from hilde.helpers.compression import backup_folder


def backup_and_restart(
    source_dir, target_folder=".", additional_files=[], restart_command=None
):
    backup_folder(source_dir, target_folder, additional_files)
