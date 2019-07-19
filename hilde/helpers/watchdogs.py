""" A watchdog keeping an eye on the time """

import os
import subprocess as sp
from time import time, strftime
from pathlib import Path
from hilde.helpers import warn, talk

_prefix = "watchdog"


def str2time(string):
    """Convert string of the shape D-HH:MM:SS to seconds

    Parameters
    ----------
    string: str
        string representing time

    Returns
    -------
    float
        time in seconds
    """

    d, h, m, s = 0, 0, 0, 0
    # split days
    l1 = string.split("-")
    if len(l1) == 2:
        ds = l1[0]
        d = int(ds)
        s2 = l1[1]
    else:
        s2 = l1[0]

    l2 = s2.split(":")

    if len(l2) == 1:
        s = int(l2[0])
    elif len(l2) == 2:
        m = int(l2[0])
        s = int(l2[1])
    elif len(l2) == 3:
        h = int(l2[0])
        m = int(l2[1])
        s = int(l2[2])

    return s + 60 * m + 3600 * h + 86400 * d


def get_time(jobid):
    """get current job time

    Parameters
    ----------
    jobid: int
        Job ID for the job watchdog is acting on

    Returns
    -------
    float
        the remaining time
    """
    try:
        squeue = sp.check_output(["squeue", "-l", "-j", jobid]).decode("utf-8")
        line = squeue.split("\n")[2]
        str_time = line.split()[5]
        return str2time(str_time)
    except sp.CalledProcessError:
        warn("subprocess.CalledProcessError was raised, return 0 time left.", level=1)
        return 0


def get_timelimit(jobid):
    """get job time limit

    Parameters
    ----------
    jobid: int
        Job ID for the job watchdog is acting on

    Returns
    -------
    float
        the time limit in seconds
    """
    squeue = sp.check_output(["squeue", "-l", "-j", jobid]).decode("utf-8")
    line = squeue.split("\n")[2]
    timelimit = line.split()[6]
    return str2time(timelimit)


class WallTimeWatchdog:
    """ Watches the walltime """

    def __init__(
        self,
        walltime=None,
        history=10,
        buffer=2,
        log="watchdog.log",
        verbose=True,
        **kwargs,
    ):
        """ Watchdog that controls the walltime everytime it is called

        Parameters
        ----------
        walltime: int
            Walltime in seconds
        history: int
            Defaults to 5. How many steps should be used to project the runtime
        buffer: int
            Defaults to 2. How many steps of buffer before watchdog should alert.
        log: str
            Path to log file
        verbose:
            If True print more logging information
        """

        if walltime is None:
            if verbose:
                talk("walltime not set, disable watchdog", prefix=_prefix)
            self.walltime = None
        else:
            self.walltime = walltime + time()

        self.buffer = buffer
        self.start_time = time()
        self.history = [time()]
        self.n_calls = 0
        self.logfile = None
        self.max_depth = history
        self.verbose = verbose

        if log is not None:
            self.logfile = Path(log)

    def __call__(self):
        """ Call the watchdog

        Returns
        -------
        bool
            Are we approaching the walltime or is a 'stop' flag present?
        """

        if self.walltime is None:
            return False

        # update history
        self.history.append(time())

        stop_file = Path("stop")
        if stop_file.exists():
            import sys

            stop_file.unlink()

            with self.logfile.open("a") as f:
                f.write("*** stop file found")
            sys.exit("*** Watchdog: stop flag was found: remove it and exit.")

        # is sufficient time left?
        time_is_up = time() + self.buffer_time > self.walltime

        # delete last step from history
        if len(self.history) > self.max_depth:
            self.history = self.history[1:]

        # log the step
        self.log()
        self.n_calls += 1

        if time_is_up and self.verbose:
            warn("Watchdog: running out of time!")

        # return information if time is up
        return time_is_up

    @property
    def increment_per_step(self):
        """ compute increment per step based on history """
        hist = self.history

        if len(hist) < 2:
            return 0

        return (hist[-1] - hist[0]) / (len(hist) - 1)

    @property
    def time_left(self):
        """ how much time is left? """
        return self.walltime - time()

    @property
    def buffer_time(self):
        """ approximate additional time the number of buffer steps would need """
        return self.increment_per_step * self.buffer

    @property
    def elapsed(self):
        """ Return elapsed time since start """
        return time() - self.start_time

    def log(self, mode="a"):
        """Log some timings

        Parameters
        ----------
        mode: str
            "a" appends to log, "w" overwrites
        """

        if self.logfile is None:
            return

        info_str = ""
        if self.n_calls == 0:
            mode = "w"
            info_str = f"# Walltime Watchdog \n"
            info_str += f"#   walltime:     {self.time_left:.0f}s\n"
            info_str += f"#   buffer steps: {self.buffer}\n"
            info_str += f"# {'Time':17s} " + " ".join(
                f"{s:>10s}"
                for s in ("n_call", "increment", "buffer_time", "time_left", "elapsed")
            )
            info_str += "\n"

        timestr = strftime("%Y/%m/%d %H:%M:%S")

        info_str += f"{timestr} {self.n_calls:10d} " + " ".join(
            f"{s:10.1f}"
            for s in (
                self.increment_per_step,
                self.buffer_time,
                self.time_left,
                self.elapsed,
            )
        )
        info_str += "\n"

        with self.logfile.open(mode) as f:
            f.write(info_str)


class SlurmWatchdog(WallTimeWatchdog):
    """Watch the slurm walltime"""

    def __init__(self, buffer=2, history=10, log="watchdog.log", verbose=True):
        """ Watchdog that controls the walltime everytime it is called

        Parameters
        ----------
        walltime: int
            Walltime in seconds
        buffer: int
            Defaults to 2. How many steps of buffer before watchdog should alert.
        history: int
            Defaults to 5. How many steps should be used to project the runtime
        log: str
            Path to log file
        verbose:
            If True print more logging information
        """
        # check jobid
        try:
            jobid = os.environ["SLURM_JOB_ID"]
            self.jobid = jobid
            walltime = get_timelimit(jobid)
            super().__init__(walltime, history, buffer, log, verbose)
        except KeyError:
            if verbose:
                msg = "seems we are not on a cluster, nothing to do for watchdog"
                talk(msg, prefix=_prefix)
            super().__init__(None, history, buffer, log, verbose=False)

    @property
    def elapsed(self):
        """ Return elapsed time since start """
        return get_time(self.jobid)
