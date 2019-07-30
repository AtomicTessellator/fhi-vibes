#! /usr/bin/env python
"""clean all files that are not under git control"""

import os
from pathlib import Path
import subprocess as sp

parent = Path(__file__).parent

all_files = sorted([file for file in parent.glob("**/*") if file != parent], reverse=True)

git_files = sp.check_output(f"git ls-files {parent}".split(), encoding="utf-8")

git_files = [Path(file) for file in git_files.split()]

for file in all_files:
    if file not in git_files:
        try:
            os.remove(file)
            print(f"{file} removed")
        except IsADirectoryError:
            continue

# clean empty folders
for file in all_files:
    if file not in git_files:
        try:
            os.rmdir(file)
            print(f"{file} removed")
        except OSError:
            continue
