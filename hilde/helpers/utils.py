""" A simple timer """

import sys
from time import time
import inspect

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
        return float(time_str[:-1])


# print in bold
def bold(text):
    """ print test in bold face """
    return "\033[1m" + text + "\033[0m"


def talk(message):
    """ https://stackoverflow.com/a/2654130/5172579 """

    curframe = inspect.currentframe()
    frame = inspect.getouterframes(curframe, 2)[1]

    file = frame[1].split("hilde")[-1]

    print(f"[hilde]: {message}\n")
    print(f".. from file hilde{file}, line {frame[2]}, function {frame[3]}")
