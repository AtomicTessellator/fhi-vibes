""" tools for backup """

import shutil
import tarfile
from glob import glob
from pathlib import Path

from vibes.helpers import talk
from vibes.helpers.aims import peek_aims_uuid
from vibes.keys import default_backup_folder  # noqa: F401

_default_files = ("aims.out", "control.in", "geometry.in")


def backup_filename(workdir="."):
    """generate a backup file name for the current backup directory"""

    counter = 0

    file = lambda counter: Path(workdir) / f"backup.{counter:05d}"

    while glob(f"{file(counter)}*"):
        counter += 1

    return file(counter)


def backup_folder(
    source_dir, target_folder=".", additional_files=None, zip=True, verbose=True
):
    """backup a folder as .tgz

    Args:
        source_dir: path to the source direcotry
        target_folder: Path to where the backups should be stored
        additional_files: Additional files to backup
        zip: True if backup folder should be compressed
        verbose: If True perform more verbose logging

    Returns:
        True if source_dir exists and is not empty
    """

    output_filename = backup_filename(target_folder)

    if not Path(source_dir).exists():
        talk(f"{source_dir} does not exists, nothing to back up.", prefix="backup")
        return False

    try:
        Path(source_dir).rmdir()
        talk(f"{source_dir} is empty, do not backup and remove.", prefix="backup")
        return False
    except OSError:
        pass

    if not Path(target_folder).exists():
        Path(target_folder).mkdir()

    # peek aims output file
    info_str = ""
    aims_uuid = peek_aims_uuid(Path(source_dir) / "aims.out")
    if aims_uuid:
        info_str = f"aims_uuid was:      {aims_uuid}"
        output_filename = f"{output_filename}.{aims_uuid[:8]}"

    if zip:
        output_filename = f"{output_filename}.tgz"
        make_tarfile(
            output_filename,
            source_dir,
            additional_files=additional_files,
            arcname=Path(source_dir).stem,
        )
    else:
        shutil.move(source_dir, output_filename)

    message = []
    message += [f"Folder:             {source_dir}"]
    message += [f"was backed up in:   {output_filename}"]
    message += [info_str]

    if verbose:
        talk(message, prefix="backup")

    return True


def make_tarfile(
    output_filename, source_dir, additional_files=None, arcname=None, only_defaults=True
):
    """create a tgz directory

    Args:
        output_filename: Path to the output file
        source_dir: Path to the source directory
        additional_files: Additional files to include in the tar file
        arcname: Path to the archive file
        only_defaults: only use the default file names for backup
    """

    outfile = Path(output_filename)

    files = Path(source_dir).glob("*")

    with tarfile.open(outfile, "w:gz") as tar:
        for file in files:
            if only_defaults and file.name not in _default_files:
                continue
            tar.add(file, arcname=Path(arcname) / file.name)

        if additional_files is not None:
            for file in additional_files:
                tar.add(file)
