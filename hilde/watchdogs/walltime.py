""" A watchdog keeping an eye on the time """

from time import time


class WallTimeWatchdog:
    def __init__(self, walltime, history=5, buffer=5):
        """ Watchdog that controls the walltime everytime it is called

        Args:
            walltime (int): Walltime in seconds
            history (int, optional):
                Defaults to 5. How many steps should be used to project the runtime
            buffer (int, optional):
                Defaults to 5. How many steps of buffer before watchdog should alert.
        """

        self.buffer = buffer
        self.walltime = walltime + time()
        self.history = [time()]

        self.max_depth = history

    def __call__(self):
        """ Call the watchdog

        Returns:
            bool: Are we approaching the walltime?
        """

        hist = self.history
        self.history.append(time())

        time_increment_per_step = (hist[-1] - hist[0]) / (len(hist) - 1)

        # delete last step
        if len(self.history) > self.max_depth:
            self.history = self.history[1:]

        time_left = self.walltime - time()

        projected_time = time() + time_increment_per_step * self.buffer

        return projected_time > self.walltime
