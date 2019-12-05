"""Green-Kubo stuff"""
from vibes.helpers import talk as _talk, Timer

_prefix = "GreenKubo"

Timer.prefix = _prefix


def talk(msg, **kw):
    """wrapper for `utils.talk` with prefix"""
    return _talk(msg, prefix=_prefix, **kw)
