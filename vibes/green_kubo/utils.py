"""Green-Kubo stuff"""
import pandas as pd

from vibes import defaults, keys
from vibes.correlation import get_correlation_time_estimate
from vibes.helpers import Timer, talk, warn

_prefix = "GreenKubo"

Timer.prefix = _prefix


def _talk(msg, **kw):
    """wrapper for `utils.talk` with prefix"""
    return talk(msg, prefix=_prefix, **kw)


def get_avalanche_data(
    series: pd.Series, Fmax: float = defaults.Fmax, verbose: bool = True, **kwargs
) -> (pd.Series, float, int):
    """Return (avalanche_function, avalanche_time in fs, avalanche_index) as dict

    Args:
        series: time series
        Fmax: max value of avalanche function to determine tmax
        verbose: be verbose
        kwargs: for F_avalanche

    Returns:
        dict
    """
    F = F_avalanche(series, verbose=verbose, **kwargs)
    t = t_avalanche(F, Fmax=Fmax, verbose=verbose)
    n = len(F[:t])

    return {keys.avalanche_function: F, keys.time_avalanche: t, keys.avalanche_index: n}


def F_avalanche(
    series, delta="auto", min_delta: int = 100, ps: bool = False, verbose: bool = True
):
    """Compute Avalanche Function (windowed noise/signal ratio)

    as defined in J. Chen et al. / Phys. Lett. A 374 (2010) 2392
    See also: Parzen, Modern Probability Theory and it's Applications, Chp. 8.6, p. 378f

    Args:
        series (pandas.Series): some time resolved data series
        delta (int): no. of time steps for windowing, or `auto`
        min_delta (int): minimal window size
        ps (bool): series.index given in ps (default: fs)
        verbose (bool): be verbose

    Returns:
        F(t, delta) = abs( sigma(series) / E(series)),
        where sigma is the standard deviation of the time series in an interval
        delta around t, and E is the expectation value around t.

    When `delta='auto'`, estimate the correlation time and use that size for binning
    the time steps for estimating std/E
    """
    series.copy().dropna()

    if ps:
        series.index *= 1000

    if delta == "auto":
        # estimate correlation time
        tau, y0 = get_correlation_time_estimate(series, verbose=True)

        delta = len(series[series.index < tau])

        if verbose:
            _talk(f"estimated correlation time: {tau/1000:7.3f} ps")
            _talk(f"-> choose delta of size:    {delta:5d} data points")

    delta = max(delta, min_delta)

    sigma = series.rolling(window=delta, min_periods=0).std()
    E = series.rolling(window=delta, min_periods=0).mean()

    F = (sigma / E).abs().dropna()

    F.name = keys.avalanche_function

    if ps:
        F.index /= 1000

    return F


def t_avalanche(series, Fmax=defaults.Fmax, verbose=True):
    """get avalanche time for series from F_avalanche

    Args:
        Fmax (float): max. allowed value for avalanche function
    """
    try:
        tmax = series[series > Fmax].index[0]
    except IndexError:
        warn("Avalanche time could not be determinded, return max. index")
        tmax = series.index.max()

    if verbose:
        _talk(f"-> avalanche time with max. F of {Fmax:2d}: {tmax:17.2f} fs")

    return tmax
