import contextlib
import os, sys

@contextlib.contextmanager
def cwd(path, mkdir=False):
    """ Purpose: change cwd intermediately
        Usage:
            with cwd(some_path):
                do so some stuff in some_path

            do so some other stuff in old cwd
    """
    CWD = os.getcwd()

    if os.path.exists(path) == False and mkdir:
        os.makedirs(path)

    os.chdir(path)
    try:
        yield
    except Exception as inst:
        print(inst)
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        os.chdir(CWD)
