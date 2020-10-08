"""Green-Kubo stuff"""
from vibes.helpers import Timer, talk, warn  # noqa: F401


_prefix = "GreenKubo"

Timer.prefix = _prefix


def _talk(msg, **kw):
    """wrapper for `utils.talk` with prefix"""
    return talk(msg, prefix=_prefix, **kw)
