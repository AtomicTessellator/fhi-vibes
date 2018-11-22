import warnings

def _warning(message, category=UserWarning, filename="", lineno=-1, *args):
    print(message)


def warn(message):
    warnings.showwarning = _warning
    warnings.warn(message)
