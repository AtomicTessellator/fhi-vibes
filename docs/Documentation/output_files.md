Output Files
===

## Trajectories

### `trajectory.son`

This file contains metadata and calculation results for a set of related calculations that, e.g., come from a geometry optimization, a phonopy  calculation, or a molecular dynamics run. The file format is [`son`](https://flokno.github.io/son/), a slight extension to `json` allowing to add data sequentially.

??? example "Example `trajectory.son`"
    ```
    {"MD": {
      "type": "molecular-dynamics",
      ...},
    "calculator": {
      "calculator": "Aims",
      "calculator_parameters": {
        "xc": "pbesol",
        "k_grid": [2, 2, 2],
        ...}},
    "atoms": {
      "cell":
        [[ 8.33141234000000e+00, -8.33141234000000e+00,  0.00000000000000e+00],
         [ 8.33141234000000e+00,  8.33141234000000e+00,  0.00000000000000e+00],
         [ 0.00000000000000e+00,  0.00000000000000e+00,  1.24971185100000e+01]],
      "positions":
        [[ 0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00],
         [ 4.16570617000000e+00,  0.00000000000000e+00,  0.00000000000000e+00],
       …}
    “primitive”: { …
    …}
    ===
    {"atoms": {
      "info": {
        "nsteps": 0,
        "dt":  4.91134739423203e-01,
        "aims_uuid": "D985353A-F8FD-4635-A939-E129A7E6E146"},
      "positions":
        [[ 0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00],
         [ 4.16570617000000e+00,  0.00000000000000e+00,  0.00000000000000e+00], …}
    ---
    ...
    ```

### `trajectory.nc`

A [`NetCDF`](https://www.unidata.ucar.edu/software/netcdf/) file containing and [`xarray.Dataset`](http://xarray.pydata.org/en/stable/io.html?highlight=netcdf#netcdf) with post-processed data

??? example "Load `trajectory.nc`"
    ```python
    >>> import xarray as xr
    >>>
    >>> ds = xr.open_dataset('trajectory.nc')
    >>>
    >>> print(ds)
    <xarray.Dataset>
    Dimensions:              (I: 180, a: 3, b: 3, time: 10001)
    Coordinates:
      * time                 (time) float64 0.0 2.0 4.0 6.0 ... 2e+04 2e+04 2e+04
    Dimensions without coordinates: I, a, b
    Data variables:
        positions            (time, I, a) float64 ...
        displacements        (time, I, a) float64 ...
        velocities           (time, I, a) float64 ...
        momenta              (time, I, a) float64 ...
        forces               (time, I, a) float64 ...
        energy_kinetic       (time) float64 ...
        energy_potential     (time) float64 ...
        stress               (time, a, b) float64 ...
        stress_kinetic       (time, a, b) float64 ...
        stress_potential     (time, a, b) float64 ...
        temperature          (time) float64 ...
        cell                 (time, a, b) float64 ...
        positions_reference  (I, a) float64 ...
        lattice_reference    (a, b) float64 ...
        pressure             (time) float64 ...
        pressure_kinetic     (time) float64 ...
        pressure_potential   (time) float64 ...
    Attributes:
        name:             trajectory
        system_name:      O3Sb2
        natoms:           180
        time_unit:        fs
        timestep:         1.9999999999999978
        nsteps:           10000
        symbols:          ['O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',...
        masses:           [ 15.999  15.999  15.999  15.999  15.999  15.999  15.99...
        atoms_reference:  {"pbc": [true, true, true],\n"cell": \n[[ 1.50410061436...
        atoms_primitive:  {"pbc": [true, true, true],\n"cell": \n[[ 5.01366871456...
        atoms_supercell:  {"pbc": [true, true, true],\n"cell": \n[[ 1.50410061436...
        volume:           3051.2387320953862
        raw_metadata:     {"MD": {\n  "type": "molecular-dynamics",\n  "md-type":...
        hash:             097714462c68b9f8cbdf08e6a29c0bfda7922c01
    ```



## Numerical Data

### `.dat` files

Files with `.dat` suffix are 1D or 2D arrays that can be read with [`numpy.loadtxt`](https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html):

??? example "Example `frequencies.dat`"
    ```
    0.000000000000000000e+00
    0.000000000000000000e+00
    0.000000000000000000e+00
    9.328537621518436795e-01
    1.319254442145902706e+00
    1.615750112078766065e+00
    1.865707524303687359e+00
    2.085924425237942970e+00
    2.285015721907639907e+00
    2.468099064244776208e+00
    2.638508884291805412e+00
    2.798561286455531594e+00
    ```

??? example "Load `frequencies.dat`"
    ```python
    >>> import numpy as np
    >>>
    >>> data = np.loadtxt('frequencies.dat')
    >>>
    >>> print(data)
    [0.         0.         0.         0.93285376 1.31925444 1.61575011
     1.86570752 2.08592443 2.28501572 2.46809906 2.63850888 2.79856129]
    ```

### `.csv` files
Files `.csv` suffix are standard [`comma-separated values`](https://en.wikipedia.org/wiki/Comma-separated_values) that can be parsed, e.g., with `pandas`.

??? example "Example `trajectory.csv`"
    ```
    time,temperature,energy_kinetic,energy_potential,pressure_kinetic,pressure_potential,pressure
    0.0,286.2001210048815,6.658958660176626,-13092419.8403197,0.0014549191863471355,-0.0001421961672687167,0.0013127230190784188
    1.9999999999999978,289.5694460043815,6.7373520438187935,-13092419.9182694,0.0014720473956910023,-0.00012147676349475005,0.0013505706321962523
    3.9999999999999956,292.46310420107915,6.80467818694143,-13092419.9839958,0.0014867575181546967,,
    5.999999999999994,294.9143609337456,6.86171106725956,-13092420.0383063,0.0014992186605137458,,
    ```

??? example  "Load `trajectory.csv`"
    ```python
    >>> import pandas as pd
    >>>
    >>> df = pd.read_csv('trajectory.csv')
    >>>
    >>> print(df)
              time  temperature  energy_kinetic  energy_potential  pressure_kinetic  pressure_potential  pressure
    0          0.0   286.200121        6.658959     -1.309242e+07          0.001455           -0.000142  0.001313
    1          2.0   289.569446        6.737352     -1.309242e+07          0.001472           -0.000121  0.001351
    2          4.0   292.463104        6.804678     -1.309242e+07          0.001487                 NaN       NaN

    [3 rows x 7 columns]
    ```

### `.json` files
These are plain [`JSON`](https://www.json.org/) files that  can be parsed with the [python builtin `json` module](https://docs.python.org/3/library/json.html)

??? example "Example `md_describe.json`"
    ```
    {
     "time": {
      "count": 10001.0,
      "mean": 9999.999999999987,
      "std": 5774.368710084239,
      "min": 0.0,
      "25%": 4999.999999999994,
      "50%": 9999.999999999987,
      "75%": 14999.999999999984,
      "max": 19999.999999999975
     },
     "temperature": {
      "count": 10001.0,
      "mean": 309.14540993328325,
      "std": 14.54050396748833,
      "min": 257.7711666913883,
      "25%": 299.2299218510854,
      "50%": 309.2638007517673,
      "75%": 319.0663317618751,
      "max": 355.05955146050223
     },
     "energy_kinetic": {
      "count": 10001.0,
      "mean": 7.19282192299974,
      "std": 0.3383108800851539,
      "min": 5.997508095931853,
      "25%": 6.962120325100335,
      "50%": 7.1955764975375205,
      "75%": 7.4236499467457895,
      "max": 8.261096699662168
     },
     "energy_potential": {
      "count": 10001.0,
      "mean": -13092416.82839676,
      "std": 2.2319201670206765,
      "min": -13092420.901731301,
      "25%": -13092419.756083699,
      "50%": -13092415.5125696,
      "75%": -13092415.1902583,
      "max": -13092414.1603322
     },
     "pressure_kinetic": {
      "count": 10001.0,
      "mean": 0.0015715632359058732,
      "std": 7.391771228878895e-05,
      "min": 0.00131039852390555,
      "25%": 0.001521157129150454,
      "50%": 0.0015721650842653186,
      "75%": 0.001621996965507344,
      "max": 0.0018049711226602926
     },
     "pressure_potential": {
      "count": 1032.0,
      "mean": 0.00027321549417593297,
      "std": 0.002184381533766492,
      "min": -0.006481330360441526,
      "25%": -0.001147865102920749,
      "50%": 0.00023375660603530354,
      "75%": 0.001799169124854924,
      "max": 0.00778187950662052
     },
     "pressure": {
      "count": 1032.0,
      "mean": 0.001844596382543767,
      "std": 0.0021801900840728618,
      "min": -0.004976051568943327,
      "25%": 0.00041292059984573813,
      "50%": 0.0018072058785148773,
      "75%": 0.0033777412146081182,
      "max": 0.009323241222005936
     }
    }
    ```

??? example "Load `md_describe.json`"
    ```python
    >>> import json
    >>>
    >>> data = json.load(open('md_describe.json'))
    >>>
    >>> pprint(data)
    {'energy_kinetic': {'25%': 6.962120325100335,
                        '50%': 7.1955764975375205,
                        '75%': 7.4236499467457895,
                        'count': 10001.0,
                        'max': 8.261096699662168,
                        'mean': 7.19282192299974,
                        'min': 5.997508095931853,
                        'std': 0.3383108800851539},
     'energy_potential': {'25%': -13092419.756083699,
                          '50%': -13092415.5125696,
                          '75%': -13092415.1902583,
                          'count': 10001.0,
                          'max': -13092414.1603322,
                          'mean': -13092416.82839676,
                          'min': -13092420.901731301,
                          'std': 2.2319201670206765},
     'pressure': {'25%': 0.00041292059984573813,
                  '50%': 0.0018072058785148773,
                  '75%': 0.0033777412146081182,
                  'count': 1032.0,
                  'max': 0.009323241222005936,
                  'mean': 0.001844596382543767,
                  'min': -0.004976051568943327,
                  'std': 0.0021801900840728618},
     'pressure_kinetic': {'25%': 0.001521157129150454,
                          '50%': 0.0015721650842653186,
                          '75%': 0.001621996965507344,
                          'count': 10001.0,
                          'max': 0.0018049711226602926,
                          'mean': 0.0015715632359058732,
                          'min': 0.00131039852390555,
                          'std': 7.391771228878895e-05},
     'pressure_potential': {'25%': -0.001147865102920749,
                            '50%': 0.00023375660603530354,
                            '75%': 0.001799169124854924,
                            'count': 1032.0,
                            'max': 0.00778187950662052,
                            'mean': 0.00027321549417593297,
                            'min': -0.006481330360441526,
                            'std': 0.002184381533766492},
     'temperature': {'25%': 299.2299218510854,
                     '50%': 309.2638007517673,
                     '75%': 319.0663317618751,
                     'count': 10001.0,
                     'max': 355.05955146050223,
                     'mean': 309.14540993328325,
                     'min': 257.7711666913883,
                     'std': 14.54050396748833},
     'time': {'25%': 4999.999999999994,
              '50%': 9999.999999999987,
              '75%': 14999.999999999984,
              'count': 10001.0,
              'max': 19999.999999999975,
              'mean': 9999.999999999987,
              'min': 0.0,
              'std': 5774.368710084239}}
    ```



## Force Constants

### `FORCE_CONSTANTS`

These are force constants in the [`phonopy` format](https://phonopy.github.io/phonopy/input-files.html?highlight=force_const#force-constants-and-force-constants-hdf5) in the compact form `(n_primitive, n_supercell, 3, 3)`. They can be parsed with `phonopy.file_IO.parse_FORCE_CONSTANTS`.

### `FORCE_CONSTANTS_remapped`

These are force constants mapped to `(3 * n_supercell, 3 * n_supercell)` shape. They can be parsed with `numpy.loadtxt` similar to [`.dat.` files](#dat-files).


## Green-Kubo results


### `greenkubo.nc`

A [`NetCDF`](https://www.unidata.ucar.edu/software/netcdf/) file containing and [`xarray.Dataset`](http://xarray.pydata.org/en/stable/io.html?highlight=netcdf#netcdf) with post-processed data

??? example "Load `trajectory.nc`"
    ```python
    >>> import xarray as xr
    >>>
    >>> ds = xr.open_dataset('greenkubo.nc')
    >>>
    >>> print(ds)
    <xarray.Dataset> Size: 34MB
    Dimensions:                                        (time: 3001, a: 3, b: 3,
                                                        s: 6, q: 108, ia: 6,
                                                        q_ir: 24, nq: 9,
                                                        q_int: 8000, i: 2, J: 216)
    Coordinates:
      * time                                           (time) float64 24kB 0.0 .....
      * nq                                             (nq) int64 72B 4 6 ... 18 20
    Dimensions without coordinates: a, b, s, q, ia, q_ir, q_int, i, J
    Data variables: (12/42)
        heat_flux                                      (time, a) float64 72kB -1....
        heat_flux_autocorrelation                      (time, a, b) float64 216kB ...
        heat_flux_autocorrelation_filtered             (time, a, b) float64 216kB ...
        heat_flux_autocorrelation_cumtrapz             (time, a, b) float64 216kB ...
        heat_flux_autocorrelation_cumtrapz_filtered    (time, a, b) float64 216kB ...
        kappa                                          (a, b) float64 72B 1.423 ....
        ...                                             ...
        interpolation_w_sq                             (s, q_int) float64 384kB 0...
        interpolation_tau_sq                           (s, q_int) float64 384kB 0...
        force_constants                                (i, J, a, b) float64 31kB ...
        thermal_conductivity_corrected                 (a, b) float64 72B 2.324 ....
        volume                                         (time) float64 24kB 5.873e...
        temperature                                    (time) float64 24kB 307.5 ...
    Attributes: (12/19)
        name:               trajectory
        system_name:        CuI
        natoms:             216
        time_unit:          fs
        timestep:           5.0
        nsteps:             12000
        ...                 ...
        sigma:              0.6311824798505474
        st_size:            749635900
        hash_raw:           4967e1958c2f78b579e3669a6582ab830b075df8
        gk_window_fs:       1195.5434782608695
        gk_prefactor:       1134529729.719537
        filter_prominence:  0.2
