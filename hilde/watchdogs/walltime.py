""" A watchdog keeping an eye on the time """

from time import time, strftime
from pathlib import Path


class WallTimeWatchdog:
    def __init__(self, walltime, history=10, buffer=3, log="watchdog.log"):
        """ Watchdog that controls the walltime everytime it is called

        Args:
            walltime (int): Walltime in seconds
            history (int, optional):
                Defaults to 5. How many steps should be used to project the runtime
            buffer (int, optional):
                Defaults to 5. How many steps of buffer before watchdog should alert.
        """

        self.buffer = buffer
        self.start_time = time()
        self.walltime = walltime + time()
        self.history = [time()]
        self.n_calls = 0
        self.logfile = None
        self.max_depth = history

        if log is not None:
            self.logfile = Path(log)

    def __call__(self):
        """ Call the watchdog

        Returns:
            bool: Are we approaching the walltime?
        """

        # update history
        self.history.append(time())

        # is sufficient time left?
        time_is_up = time() + self.buffer_time > self.walltime

        # delete last step from history
        if len(self.history) > self.max_depth:
            self.history = self.history[1:]

        # log the step
        self.log()
        self.n_calls += 1

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
        return self.walltime - time()

    @property
    def buffer_time(self):
        """ approximate additional time the number of buffer steps would need """
        return self.increment_per_step * self.buffer

    @property
    def elapsed(self):
        """ Return elapsed time since start """
        return time() - self.start_time

    def log(self):
        """ Log some timings """

        if self.logfile is None:
            return

        info_str = ""
        if self.n_calls == 0:
            info_str = f"\n# Walltime Watchdog \n"
            info_str += f"#   walltime:     {self.time_left:.0f}s\n"
            info_str += f"#   buffer steps: {self.buffer}\n"
            info_str += "# {:17s} {:>7s} {:>10s} {:>10s} {:>10s}\n".format(
                "Time", "n_calls", "increment", "buffer_time", "time_left"
            )

        timestr = strftime("%Y/%m/%d %H:%M:%S")
        info = (self.n_calls, self.increment_per_step, self.buffer_time, self.time_left)
        info_str += "{} {:7d} {:10.1f} {:10.1f} {:10.1f}\n".format(timestr, *info)

        with self.logfile.open("a") as f:
            f.write(info_str)
