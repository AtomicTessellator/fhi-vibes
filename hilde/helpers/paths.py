import contextlib
import os, sys
import functools

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
    try:
        yield
    except Exception as inst:
        print(inst)
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        os.chdir(CWD)

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
