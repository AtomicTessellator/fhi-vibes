""" A simple timer """

import sys
from time import time

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    tqdm = lambda x, *args: x


def progressbar(func):
    """ show progressbar when looping """
    return tqdm(func, file=sys.stdout)


class Timer:
    def __init__(self, message=None):
        self.time = time()

        if message:
            print(message)

    def __call__(self, info_str=""):
        """ print how much time elapsed """

        time_str = f"{time() - self.time:.3f}s"

        if info_str.strip():
            print(f".. {info_str} in {time_str}")
        else:
            print(f".. time elapsed: {time_str}")


# print in bold
def bold(text):
    """ print test in bold face """
    return "\033[1m" + text + "\033[0m"
