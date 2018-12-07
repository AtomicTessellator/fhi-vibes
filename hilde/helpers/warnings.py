import sys
import inspect


def warn(message, level=1):
    " https://stackoverflow.com/a/2654130/5172579 "

    curframe = inspect.currentframe()
    frame = inspect.getouterframes(curframe, 2)[1]

    if level == 1:
        typ = "Message"
    elif level == 2:
        typ = "Warning"
    elif level == 3:
        typ == "Error"

    stars = "*" + "*" * level

    print(f"{stars} {typ} from file {frame[1]}, line {frame[2]}, function {frame[3]}:")
    print(f"--> {message}")

    if level == 3:
        sys.exit("stop")
