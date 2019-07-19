"""A simple timer"""

import signal
from time import time, strftime
import inspect
import click
from son.progressbar import progressbar
from hilde.helpers.warnings import warn


# print in bold
def bold(text):
    """ print text in bold face """
    return "\033[1m" + text + "\033[0m"


def talk(message, prefix=None, verbosity=1):
    """hilde message output. Use instead of print. Sensitive to CLI context

    https://stackoverflow.com/a/2654130/5172579

    Args:
        message (str): message to print
        prefix (str): prefix for the message
        verbosity (int): verbosity level (0, 1, 2)
    """
    # see if we are in a CLI context
    verbose = 1
    try:
        ctx = click.get_current_context()
        verbose = ctx.obj.verbose
    except (RuntimeError, AttributeError):
        pass

    if verbose == 1:
        print_msg(message, prefix=prefix)
    elif verbose > 1:
        curframe = inspect.currentframe()
        frame = inspect.getouterframes(curframe, 2)[1]

        file = frame[1].split("hilde")[-1][1:]

        timestr = strftime("%H:%M:%S %Y/%m/%d")

        print(f"[{timestr} from {file}, l. {frame[2]} in {frame[3]}()]", flush=True)
        print_msg(message, prefix=prefix, indent=2)
        print()


def print_msg(message, prefix=None, indent=0):
    """print for talk

    Args:
        message (str): message to print
        prefix (str): prefix for message
        indent (int): number of spaces to indent by
    """
    indent = indent * " "
    if not prefix:
        pref = "[hilde]"
    else:
        pref = f"[{prefix}] "
    if isinstance(message, list):
        for msg in message:
            print(f"{indent}{pref:12}{msg}", flush=True)
    else:
        print(f"{indent}{pref:12}{message}", flush=True)


class Timer:
    """simple timer with Timeout function"""

    def __init__(self, message=None, use_talk=True, timeout=None, verbose=True):
        """Initialize

        Args:
            message: Message to print at initialization
            use_talk: If true use talk and not print
            timeout: set a timeout after which a TimeoutError is raised
            verboes: be verbose

        Timeout inspired by
            https://www.jujens.eu/posts/en/2018/Jun/02/python-timeout-function/
        """
        self.time = time()
        self.verbose = verbose

        if use_talk and talk:
            self.print = talk
        else:
            self.print = lambda msg: print(msg, flush=True)

        self.message = message
        if message and verbose:
            self.print(message)

        self.timeout = timeout
        if timeout:
            signal.signal(signal.SIGALRM, self.raise_timeout)
            signal.alarm(timeout)

    def __call__(self, info_str=""):
        """print how much time elapsed, optionally print `info_str`"""
        time_str = f"{time() - self.time:.3f}s"

        if info_str.strip() and self.verbose:
            self.print(f".. {info_str} in {time_str}")
        elif self.verbose:
            self.print(f".. time elapsed: {time_str}")

        # stop signal alarm if it was initialized
        if self.timeout:
            signal.signal(signal.SIGALRM, signal.SIG_IGN)

        return float(time_str[:-1])

    def raise_timeout(self, signum, frame):
        """raise TimeoutError"""
        warn(f"Timeout of {self.timeout}s is approaching, raise TimeoutError", level=1)
        raise TimeoutError


def raise_timeout(signum, frame):
    """raise TimeoutError"""
    raise TimeoutError
