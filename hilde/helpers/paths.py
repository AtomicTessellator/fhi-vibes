import contextlib
import os
import sys
import functools
from pathlib import Path
from warnings import warn


@contextlib.contextmanager
def cwd(path, mkdir=False, debug=False):
    """ Purpose: change cwd intermediately
        Usage:
            with cwd(some_path):
                do so some stuff in some_path

            do so some other stuff in old cwd
    """
    CWD = os.getcwd()

    if os.path.exists(path) == False and mkdir:
        os.makedirs(path)

    if debug:
        os.chdir(path)
        yield
        os.chdir(CWD)
        return

    os.chdir(path)
    yield
    os.chdir(CWD)


def move_to_dir(file, folder):
    file = Path(file)
    if file.exists():
        folder.mkdir()
        file.rename(folder / file)
    else:
        warn(f"** move_to_dir: {file} does not exist.")


# would be nice to have?
# def decor_cwd(path, mkdir=False, debug=False):
#     def decorator_cwd(func):
#         @functools.wraps(func)
#         def execute_func(path, *args, mkdir=False, debug=False, **kwargs):
#             with cwd(path, mkdir=mkdir, debug=debug):
#                 values = func(*args, **kwargs)
#             return values
#         return execute_func
#     return decorator_cwd
