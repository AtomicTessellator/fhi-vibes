import contextlib
import os, sys

@contextlib.contextmanager
def cd(path):
    """ Purpose: change cwd intermediately
        Usage:
            with cd(some_path):
                do so some stuff in some_path

            do so some other stuff in old cwd
    """
    CWD = os.getcwd()

    os.chdir(path)
    try:
        yield
    except:
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        os.chdir(CWD)
