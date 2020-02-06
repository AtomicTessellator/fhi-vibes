""" test green kubo cumulative kappa"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate as si
from scipy import signal as sl

from vibes.correlation import get_autocorrelation
from vibes.integrate import get_cumtrapz
from vibes.settings import ConfigDict

parent = Path(__file__).parent
folder = parent / "lammps"

# read info
info = ConfigDict(folder / "info.cfg").info

# read reference
flux_df = pd.read_csv(folder / "flux.csv.bz2", index_col=info.index_col)
kappa_df = pd.read_csv(folder / "kappa.csv.bz2", index_col=info.index_col)


def _ref_corr(series, step=10, nmax=100000, window=True):
    """return correlation function as pd.Series"""
    J = series.iloc[:nmax:step]
    time = J.index

    Nt = len(J)

    corr = sl.correlate(J, J)[Nt - 1 :]

    # Normalize
    corr /= np.arange(Nt, 0, -1)
    if window:  # window
        corr *= sl.windows.hann(2 * Nt)[Nt:]

    return pd.Series(corr, index=time)


def test_j_corr(flux=flux_df.flux):
    c1 = _ref_corr(flux)
    c2 = get_autocorrelation(flux)

    assert (c1 - c2).std() < 1e-5


def test_kappa(flux=flux_df.flux, ref_kappa=kappa_df.kappa):
    time = np.arange(2000, 100000, 2000)

    step = info.step
    tmax = info.tmax
    scale = info.scale

    k = []
    for nmax in time:
        # corr = get_corr(flux_df.flux, step=step, nmax=nmax, window=False)[:]
        J = flux.iloc[:nmax:step]
        corr = get_autocorrelation(J, window=False, verbose=False)
        k.append(si.trapz(corr[:tmax]) * scale * step)

    kappa = pd.Series(k, index=time)

    assert (kappa - ref_kappa).std() < 0.1

    c = get_autocorrelation(flux, verbose=False) * scale
    k2 = get_cumtrapz(c)

    assert (ref_kappa.iloc[-1] - k2[tmax]) < 0.001


if __name__ == "__main__":
    test_j_corr()
    test_kappa()
