""" Summarize output from ASE.md class (in md.log) """

import subprocess
from argparse import ArgumentParser
from hilde import Settings
from hilde.helpers import Timer

new_addr = "https://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/coe/api/uploads/?curl=True"
old_addr = "http://nomad-repository.eu:8000"


def upload_command(folder, token, legacy):
    """Generate the NOMAD upload command

    Args:
        folder: The folder to upload
        token: The NOMAD token
        legacy: use old NOMAD

    Returns:
        str: The upload command
    """
    if legacy:
        cmd = (
            f"tar cf - {folder} | curl -XPUT -# -HX-Token:{token} "
            f"-N -F file=@- {old_addr} | xargs echo"
        )
    else:
        cmd = (
            f"tar czf - {folder} | "
            f'curl -X PUT -H "X-Token: {token}" -F file=@- "{new_addr}"'
            "| xargs echo"
        )

    return cmd


def nomad_upload(folders, token=None, legacy=False, dry=False):
    """upload folders with calculations to NOMAD

    Args:
        folders: The folders to upload
        token: The NOMAD token
        legacy: use old NOMAD
        dry: only show upload command
    """
    timer = Timer()

    settings = Settings()

    if not token and "nomad" in settings:
        token = settings.nomad.token

    if token is None:
        exit("** Token is missing, chech your .hilderc or provide manually")

    # from ASE
    if not folders:
        exit("No folders specified -- another job well done!")

    for ii, folder in enumerate(folders):

        cmd = upload_command(folder, token, legacy)

        if dry:
            print(f"Upload command {ii+1}:\n{cmd}")
        else:
            print(f"Upload folder {folder:30} ({ii+1} of {len(folders)})")
            print(f"  with command  {cmd}")

            subprocess.check_call(cmd, shell=True)

    if not dry:
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
