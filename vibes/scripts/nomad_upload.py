""" Summarize output from ASE.md class (in md.log) """

import subprocess
from argparse import ArgumentParser

from vibes import Settings
from vibes.helpers import Timer

new_addr = "http://repository.nomad-coe.eu/app/api/uploads/?token="
old_addr = "http://nomad-repository.eu:8000"


def upload_command(files: list, token: str, legacy: bool = False, name: str = None):
    """Generate the NOMAD upload command

    Args:
        folder: The folder to upload
        token: The NOMAD token
        legacy: use old NOMAD
        name: nomad upload name

    Returns:
        str: The upload command
    """
    name_str = ""
    if name is not None:
        name_str = f"&name={name}"

    file_str = " ".join((str(f) for f in sorted(files)))
    cmd = f"tar czf - {file_str} | "

    if legacy:
        cmd = +(f"curl -XPUT -# -HX-Token:{token} " f"-N -F file=@- {old_addr}")
    else:
        cmd += f'curl "{new_addr}{token}{name_str}" -T'

    cmd += " - | xargs echo"

    return cmd


def nomad_upload(files=None, token=None, legacy=False, dry=False, name: str = None):
    """upload folders with calculations to NOMAD

    Args:
        files: The files to upload
        token: The NOMAD token
        legacy: use old NOMAD
        dry: only show upload command
        name: nomad upload name
    """
    timer = Timer()

    settings = Settings()

    if not token and "nomad" in settings:
        token = settings.nomad.token

    if token is None:
        exit("** Token is missing, chech your .vibesrc or provide manually")

    # from ASE
    if files is None:
        exit("No folders specified -- another job well done!")

    cmd = upload_command(files, token, legacy, name=name)

    print(f"Upload command:\n{cmd}")

    if not dry:
        subprocess.check_call(cmd, shell=True)
        timer(f"Nomad upload finished")


def main():
    """ main routine """
    parser = ArgumentParser(description="Upload folder to Nomad")
    parser.add_argument("folders", nargs="+", help="folder containing data to upload")
    parser.add_argument("--token", help="Nomad token for uploads")
    parser.add_argument("--dry", action="store_true", help="only show command")
    args = parser.parse_args()

    nomad_upload(args.folders, args.token, args.dry)


if __name__ == "__main__":
    main()
