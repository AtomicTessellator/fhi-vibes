import contextlib
import os, sys

@contextlib.contextmanager
def cd(path):
    CWD = os.getcwd()

    os.chdir(path)
    try:
        yield
    except:
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        os.chdir(CWD)