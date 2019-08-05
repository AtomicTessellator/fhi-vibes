""" test green kubo cumulative kappa"""
from pathlib import Path
import xarray as xr
import hilde.green_kubo.heat_flux as hf

parent = Path(__file__).parent


def test_kappa(datafile=parent / "heat_flux.nc"):
    # read the heat flux dataset
    DS = xr.load_dataset(datafile)

    # compute cumulative kappa
    kappa = hf.get_cumulative_kappa(DS)

    k = kappa.sum(axis=(1, 2)).to_series() / 3

    assert abs(k.iloc[-1] - 0.4980131763627837) < 1e-5


def test_j_corr(datafile=parent / "heat_flux.nc"):
    # read the heat flux dataset
    DS = xr.load_dataset(datafile)

    # Compute heat flux autoccorrelation function
    jcorr = hf.get_heat_flux_aurocorrelation(DS)

    j = jcorr.sum(axis=(1, 2)).to_series()

    assert abs(j[0] - 27.89026177401778) < 1e-5


if __name__ == "__main__":
    test_kappa()
