""" tools for compression """

import shutil
import tarfile
from pathlib import Path

from hilde.helpers import talk
from hilde.helpers.aims import peek_aims_uuid


def backup_filename(workdir=".", zip=False):
    """ generate a backup file name """

    counter = 0
    if zip:
        suffix = ".tgz"
    else:
        suffix = ""

    file = lambda counter: Path(workdir) / f"backup.{counter:05d}{suffix}"

    while file(counter).exists():
        counter += 1

    return file(counter)


def backup_folder(
    source_dir, target_folder=".", additional_files=[], zip=False, verbose=True
):
    """ backup a folder as .tgz """

    output_filename = backup_filename(target_folder, zip=zip)

    if not Path(source_dir).exists():
        talk(f"{source_dir} does not exists, nothing to back up.")
        return False

    try:
        Path(source_dir).rmdir()
        talk(f"{source_dir} is empty, do not backup and remove.")
        return False
    except OSError:
        pass

    # peek aims output file
    info_str = ""
    aims_uuid = peek_aims_uuid(Path(source_dir) / "aims.out")
    if aims_uuid:
        info_str = f"aims_uuid was:      {aims_uuid}"

    if zip:
        make_tarfile(output_filename, source_dir, additional_files=additional_files)
    else:
        shutil.move(source_dir, output_filename)

    message = []
    message += [f"Folder:             {source_dir}"]
    message += [f"was backed up in:   {output_filename}."]
    message += [info_str]

    if verbose:
        talk(message)

    return True


def make_tarfile(output_filename, source_dir, additional_files=[], arcname=None):
    """ create a tgz directory """

    outfile = Path(output_filename)

    with tarfile.open(outfile, "w:gz") as tar:
        tar.add(source_dir, arcname=arcname)
        for file in additional_files:
            tar.add(file)
