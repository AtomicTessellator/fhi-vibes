"""A simple timer"""
import sys
import time
import threading
import itertools

import signal
import inspect
import click
from son.progressbar import progressbar
from hilde.helpers.warnings import warn


# print in bold
def bold(text):
    """ print text in bold face """
    return "\033[1m" + text + "\033[0m"


def talk(message, prefix=None, verbose=True):
    """hilde message output. Use instead of print. Sensitive to CLI context

    https://stackoverflow.com/a/2654130/5172579

    Args:
        message (str): message to print
        prefix (str): prefix for the message
        verbosity (int): verbosity level (0, 1, 2)
    """
    # see if we are in a CLI context
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

        timestr = time.strftime("%H:%M:%S %Y/%m/%d")

        print(f"[{timestr} from {file}, l. {frame[2]} in {frame[3]}()]", flush=True)
        print_msg(message, prefix=prefix, indent=2)
        print()


def print_msg(message, prefix=None, indent=0, width=15):
    """print for talk

    Args:
        message (str): message to print
        prefix (str): prefix for message
        indent (int): number of spaces to indent by
        width (int): width of prefix
    """
    indent = indent * " "
    if not prefix:
        pref = "[hilde]"
    else:
        pref = f"[{prefix}] "
    if isinstance(message, list):
        for msg in message:
            print(f"{indent}{pref:{width}}{msg}", flush=True)
    else:
        print(f"{indent}{pref:{width}}{message}", flush=True)


class Timer:
    """simple timer with Timeout function"""

    prefix = None

    def __init__(self, message=None, timeout=None, prefix=None, verbose=True):
        """Initialize

        Args:
            message: Message to print at initialization
            timeout: set a timeout after which a TimeoutError is raised
            prefix: prefix for `talk'
            verboes: be verbose

        Timeout inspired by
            https://www.jujens.eu/posts/en/2018/Jun/02/python-timeout-function/
        """
        self.time = time.time()
        self.verbose = verbose

        self.print = talk

        if prefix:
            self.prefix = prefix

        self.message = message
        if message and verbose:
            self.print(message, prefix=self.prefix)

        self.timeout = timeout
        if timeout:
            signal.signal(signal.SIGALRM, self.raise_timeout)
            signal.alarm(timeout)

    def __call__(self, info_str=""):
        """print how much time elapsed, optionally print `info_str`"""
        time_str = f"{time.time() - self.time:.3f}s"

        if info_str.strip() and self.verbose:
            self.print(f".. {info_str} in {time_str}", prefix=self.prefix)
        elif self.verbose:
            self.print(f".. time elapsed: {time_str}", prefix=self.prefix)

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


class Spinner:
    """Spinner for command line feedback

    Inspired by:
        https://stackoverflow.com/a/39504463/5172579
    """

    busy = False
    delay = 0.1

    spinning_cursor = itertools.cycle(["-", "|", "\\", "/"])

    def __init__(self, prefix="working", delay=None, file=sys.stdout, verbose=True):
        self.prefix = prefix
        self.file = file
        self.verbose = verbose
        self.msg = ""
        self.isatty = True
        self.spinner_generator = self.spinning_cursor
        if delay and float(delay):
            self.delay = delay

        if not hasattr(file, "isatty") or not file.isatty():
            self.isatty = False

    def get_msg(self):
        msg = f"{self.prefix}: "
        if self.busy:
            if self.isatty:
                msg += f"{next(self.spinner_generator)}"
            else:
                msg += f"working"
        else:
            msg += "finished." + "\n"
        return msg

    def print(self, newline=False):
        if self.verbose:
            self.msg = self.get_msg()
            if newline:
                self.file.write("\n")
            self.file.write(self.msg)
            self.file.flush()

    def rewind(self):
        if self.verbose:
            self.file.write(len(self.msg) * "\b")
            self.file.flush()

    def spinner_task(self):
        while self.busy:
            self.print()
            time.sleep(self.delay)
            self.rewind()
        self.print()

    def __enter__(self):
        self.busy = True
        if self.isatty:
            threading.Thread(target=self.spinner_task).start()
        else:
            self.print()

    def __exit__(self, exception, value, tb):
        self.busy = False
        if self.isatty:
            time.sleep(self.delay)
        else:
            self.print(newline=True)

        if exception is not None:
            return exception
