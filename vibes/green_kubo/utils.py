"""Green-Kubo stuff"""
import numpy as np
from scipy.optimize import curve_fit

from vibes.helpers import Timer
from vibes.helpers import talk as _talk

_prefix = "GreenKubo"

Timer.prefix = _prefix


def talk(msg, **kw):
    """wrapper for `utils.talk` with prefix"""
    return _talk(msg, prefix=_prefix, **kw)


def F_avalanche(series, delta="auto", verbose=True):
    """Compute Avalanche Function (windowed noise/signal ratio)

    as defined in J. Chen et al. / Phys. Lett. A 374 (2010) 2392
    See also: Parzen, Modern Probability Theory and it's Applications, Chp. 8.6, p. 378f

    Args:
        series (pandas.Series): some time resolved data series
        delta (int): no. of time steps for windowing, or `auto`
        verbose (bool): be verbose

    Returns:
        F(t, delta) = abs( sigma(series) / E(series)),
        where sigma is the standard deviation of the time series in an interval
        delta around t, and E is the expectation value around t.

    When `delta='auto'`, estimate the correlation time and use that size for binning
    the time steps for estimating std/E
    """
    if delta == "auto":
        # estimate correlation time

        exp = lambda x, y0, tau: y0 * np.exp(-x / tau)

        # where is jc drops below 0 the first time
        idx = series[series < series.max() / np.e].index[0]
        x = series[:idx]
        y0 = x.iloc[0]
        tau = idx
        bounds = ([0.75 * y0, 1.5 * y0], [0.5 * idx, 2 * idx])
        (y0, tau), _ = curve_fit(exp, x.index, x, bounds=bounds)

        delta = len(series[series.index < tau])

        if verbose:
            talk(f"Pre-est. correlation time (drop below 1/e): {idx:10.2f} fs")
            talk(f".. use {len(x)} data points to fit exponential")
            talk(f".. estimated correlation time:              {tau:10.2f} fs")
            talk(f"-> choose delta of size: {delta:20d} data points")

        assert 0.5 < idx / tau < 2, (idx, tau)

    sigma = series.rolling(window=delta, min_periods=0).std()
    E = series.rolling(window=delta, min_periods=0).mean()

    F = (sigma / E).abs().dropna()

    return F


def t_avalanche(series, Fmax=1, verbose=True, **kwargs):
    """get avalanche time for series from F_avalanche

    Args:
        Fmax (float): max. allowed value for avalanche function
    """
    F = F_avalanche(series, verbose=verbose, **kwargs)

    tmax = F[F > Fmax].index[0]

    if verbose:
        talk(f"-> avalanche time with max. F of {Fmax:2d}: {tmax:17.2f} fs")

    return tmax
