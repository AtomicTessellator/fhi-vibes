#! /usr/bin/env python
"""clean all files that are not under git control"""

import os
import subprocess as sp
from pathlib import Path

import click

parent = Path(__file__).parent

files = parent.glob("**/*")
all_files = sorted((file for file in files if file != parent), reverse=True)
git_files = sp.check_output(f"git ls-files {parent}".split(), encoding="utf-8")
git_files = [Path(file) for file in git_files.split()]


@click.command()
@click.option("--ignore", default=".coverage", show_default=True)
def cleanup(all_files=all_files, git_files=git_files, ignore=None):
    for file in all_files:
        if file not in git_files and file.name != ignore:
            try:
                os.remove(file)
                print(f"{file} removed")
            except IsADirectoryError:
                continue

    # clean empty folders
    for file in all_files:
        if file not in git_files and ignore not in file.name:
            try:
                os.rmdir(file)
                print(f"{file} removed")
            except OSError:
                continue


if __name__ == "__main__":
    cleanup()
