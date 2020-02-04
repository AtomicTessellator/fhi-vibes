import numpy as np
import scipy.signal as sl

from vibes.helpers.warnings import warn


def get_correlation(f1, f2, normalize=True, window=True):
    """Compute correlation function for signal f1 and signal f2

    Reference:
        https://gitlab.com/flokno/python_recipes/-/blob/master/mathematics/
        correlation_function/autocorrelation.ipynb

    Args:
        f1: signal 1
        f2: signal 2
        normalize: normalize correlation function
        window: apply Hann window
    Returns:
        the correlation function
    """
    Nt = min(len(f1), len(f2))

    if Nt != max(len(f1), len(f2)):
        msg = "The two signals are not equally long: "
        msg += f"len(f1), len(f2) = {len(f1)}, {len(f2)}"
        warn(msg)

    corr = sl.correlate(f1[:Nt], f2[:Nt])[Nt - 1 :]

    if normalize:
        corr /= np.arange(Nt, 0, -1)

    if window:
        corr *= sl.windows.hann(2 * Nt)[Nt:]

    return corr


def get_autocorrelation(f, **kwargs):
    """compute autocorrelation function for signal f

    Frontend to `vibes.helpers.correlation.correlation`
    """

    return get_correlation(f, f, **kwargs)
