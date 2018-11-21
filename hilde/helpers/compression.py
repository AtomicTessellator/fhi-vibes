""" tools for compression """

import tarfile
from pathlib import Path
from hilde.helpers import cwd


def backup_filename(workdir="."):
    """ generate a backup file name """

    counter = 0
    file = lambda counter: Path(workdir) / f"backup.{counter:05d}.tgz"

    while file(counter).exists():
        counter += 1

    return file(counter)


def backup_folder(source_dir, target_folder=".", additional_files=[]):
    """ backup a folder as .tgz """

    output_filename = backup_filename(target_folder)

    make_tarfile(
        output_filename,
        source_dir,
        additional_files=additional_files,
    )


def make_tarfile(output_filename, source_dir, additional_files=[], arcname=None):
    """ create a tgz directory """

    outfile = Path(output_filename)

    with tarfile.open(outfile, "w:gz") as tar:
        tar.add(source_dir, arcname=arcname)
        for file in additional_files:
            tar.add(file)
