""" A simple timer """

from time import time


class Timer:
    def __init__(self):
        self.time = time()

    def __call__(self, info_str=""):
        """ print how much time elapsed """

        time_str = f"{time() - self.time:.3f}s"

        if info_str.strip():
            print(f".. {info_str} finished in {time_str}")
        else:
            print(f".. time elapsed: {time_str}")
