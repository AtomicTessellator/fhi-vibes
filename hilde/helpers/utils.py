"""A simple timer"""

import sys
from time import time, strftime
import inspect
import click

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    tqdm = lambda x, *args, **kwargs: x


def progressbar(func):
    """ show progressbar when looping """
    return tqdm(func, file=sys.stdout)


# print in bold
def bold(text):
    """ print test in bold face """
    return "\033[1m" + text + "\033[0m"


def talk(message):
    """hilde message output. Use instead of print. Sensitive to CLI context

     https://stackoverflow.com/a/2654130/5172579

     """
    # see if we are in a CLI context
    verbose = 1
    try:
        ctx = click.get_current_context()
        verbose = ctx.obj.verbose
    except RuntimeError:
        pass

    if verbose == 1:
        print_msg(message)
    elif verbose > 1:
        curframe = inspect.currentframe()
        frame = inspect.getouterframes(curframe, 2)[1]

        file = frame[1].split("hilde")[-1][1:]

        timestr = strftime("%H:%M:%S %Y/%m/%d")

        print(f"[{timestr} from {file}, l. {frame[2]} in {frame[3]}()]", flush=True)
        print_msg(message, indent=2)
        print()


def print_msg(message, indent=0):
    """print for talk"""
    indent = indent * " "
    if isinstance(message, list):
        for msg in message:
            print(f"{indent}{msg}", flush=True)
    else:
        print(f"{indent}{message}", flush=True)


class Timer:
    """simple timer"""

    def __init__(self, message=None, use_talk=True):
        self.time = time()

        if use_talk:
            self.print = talk
        else:
            self.print = lambda msg: print(msg, flush=True)

        if message:
            self.print(message)

    def __call__(self, info_str=""):
        """ print how much time elapsed """

        time_str = f"{time() - self.time:.3f}s"

        if info_str.strip():
            self.print(f".. {info_str} in {time_str}")
        else:
            self.print(f".. time elapsed: {time_str}")
        return float(time_str[:-1])
