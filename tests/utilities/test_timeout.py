import time

from vibes.helpers.utils import Timeout


def test_timeout():
    timeout = Timeout(1, kill=False)

    time.sleep(0.5)

    timeout()

    print("made it to here")

    try:
        time.sleep(2)
        print("shouldn't make it to here")
    except TimeoutError:
        pass
