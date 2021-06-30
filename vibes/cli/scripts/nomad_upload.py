""" Summarize output from ASE.md class (in md.log) """

import json
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

import click

from vibes import Settings
from vibes.helpers import Timer, talk


_addr = "http://nomad-lab.eu/prod/rae/api/uploads/?token="
upload_folder_dry = "nomad_upload_dry"


def upload_command(
    files: list, token: str, name: str = None, nojson: bool = False,
) -> str:
    """Generate the NOMAD upload command

    Args:
        folder: The folder to upload
        token: The NOMAD token
        name: nomad upload name
        nojson: don't retrieve json summary

    Returns:
        str: The upload command
    """
    name_str = ""
    if name is not None:
        name_str = f"&name={name}"

    json_str = " -H 'Accept: application/json'"
    if nojson:
        json_str = ""

    file_str = " ".join((str(f) for f in sorted(files)))
    cmd = f"tar cf - {file_str} | "

    cmd += f'curl "{_addr}{token}{name_str}" {json_str} -T'
    cmd += " - | xargs -0 echo"

    return cmd


@click.command(context_settings={"show_default": True})
@click.argument("files", nargs=-1, type=Path)
@click.option("-t", "--token", type=str, help="NOMAD upload token")
@click.option("--dry", is_flag=True, help="dry run (no upload)")
@click.option("-n", "--name", type=str, help="custom name for the upload")
@click.option("--tmp_folder", default="nomad_upload", help="prefix for tmp folder")
@click.option("--summary_file", default="nomad.json", help="name for summary file")
def upload(files, token, dry, name, tmp_folder, summary_file):
    """upload folders with calculations to NOMAD

    Args:
        files: The files to upload
        token: The NOMAD token
        dry: only show upload command
        name: nomad upload name
        tmp_folder: name for (tmp) folder used to upload
        summary_file: if given, write nomad summary json file
    """
    timer = Timer("Perform Nomad upload")

    settings = Settings()

    if not token and "nomad" in settings:
        token = settings.nomad.token

    if token is None:
        exit("** Token is missing, chech your .vibesrc or provide manually")

    # from ASE
    if files is None:
        exit("No folders specified -- another job well done!")

    # copy files to tmpdir, unzip if necessary
    tmp_files = []
    if dry:
        tmp_dir = Path(tmp_folder)
        tmp_dir.mkdir(exist_ok=True)
    else:
        tmp_dir = Path(tempfile.mkdtemp(prefix=tmp_folder, dir="."))
    for file in files:
        path = tmp_dir / file
        if path.suffix == ".tgz":
            path = tmp_dir / Path(file).parent / path.stem

        if path.exists():
            talk(f'Clean up "{path}"')
            shutil.rmtree(path)

        talk(f'Exctract "{file}" into "{path}"')
        try:
            with tarfile.open(file) as f:
                f.extractall(path=path)
        except IsADirectoryError:
            shutil.copytree(file, path)
        except tarfile.ReadError:
            shutil.copy(file, path)
        tmp_files.append(path)

    # upload
    cmd = upload_command(tmp_files, token, name=name)
    print(f"Upload command:\n{cmd}")

    if dry:
        talk(f".. dry run requested, stop here.")
        return

    outp = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    shutil.rmtree(tmp_dir)
    timer(f"Nomad upload finished")

    summary = json.dumps(json.loads(outp.stdout), indent=2)
    if summary_file is not None:
        Path(summary_file).write_text(summary)
        talk(f"JSON summary written to {summary_file}")
    else:
        print("# json summary")
        print(summary)


if __name__ == "__main__":
    upload()
