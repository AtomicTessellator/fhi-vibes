""" A watchdog keeping an eye on the time """

from time import time


class WallTimeWatchdog:
    def __init__(self, walltime, history=10, buffer=5):
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

        self.max_depth = history

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

