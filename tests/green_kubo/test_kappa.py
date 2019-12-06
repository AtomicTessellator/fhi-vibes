""" test green kubo cumulative kappa"""
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import vibes.green_kubo.heat_flux as hf

parent = Path(__file__).parent


def test_kappa(datafile=parent / "heat_flux.nc"):
    # read the heat flux dataset
    DS = xr.load_dataset(datafile)

    # compute cumulative kappa
    kappa = hf.get_kappa(DS).kappa

    k = np.trace(kappa, axis1=1, axis2=2) / 3
    k = pd.Series(k, index=DS.time)

    assert abs(k.iloc[10] - 0.3401563148255578) < 1e-5, k.iloc[10]


def test_j_corr(datafile=parent / "heat_flux.nc"):
    # read the heat flux dataset
    DS = xr.load_dataset(datafile)

    # Compute heat flux autoccorrelation function
    jcorr = hf.get_kappa(DS).Jcorr

    # x component
    jx = jcorr[:, 0, 0].to_series() * 1000

    assert abs(jx[0] - 10.560280496761488) < 1e-5, jx[0]


if __name__ == "__main__":
    test_kappa()
    test_j_corr()
