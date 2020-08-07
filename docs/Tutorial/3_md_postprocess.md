# Postprocessing

In this tutorial, you learn how to access the data you obtained from the MD calculation. It further showcases the interplay of `xarray`, `pandas`, `numpy` , and `matplotlib` for convenient exploration of data. We begin by showing how to convert the `xarray.Dataset` to plain `numpy` arrays if you prefer to work that way.

!!! warning
	We assume you successfully ran the simulation of LDA-Silicon at $300\,{\rm K}$ from the [previous chapter](3_md_ab_initio.md) and created a `trajectory.nc` dataset. You can can also take the reference data from [here](https://gitlab.com/vibes-developers/vibes-tutorial-files/-/tree/master/3_molecular_dynamics/ab_initio/si_8).

## Load and Inspect Dataset

The data in `trajectory.nc` is stored in an [`xarray.Dataset`](http://xarray.pydata.org/en/stable/data-structures.html#dataset), which is a collection of labelled `numpy` arrays. Let's open the file and inspect the content:

```python
import xarray as xr

# load the trajectory dataset
dataset = xr.load_dataset("trajectory.nc")

print(dataset)
```

??? info "`print(dataset)`"

    ```
    <xarray.Dataset>
    Dimensions:                    (I: 8, Ia: 24, Jb: 24, a: 3, b: 3, time: 2501)
    Coordinates:
      * time                       (time) float64 0.0 4.0 8.0 ... 9.996e+03 1e+04
    Dimensions without coordinates: I, Ia, Jb, a, b
    Data variables:
        positions                  (time, I, a) float64 0.0 0.0 0.0 ... 4.073 1.298
        displacements              (time, I, a) float64 -2.98e-11 ... -0.05708
        velocities                 (time, I, a) float64 0.03397 0.01561 ... 0.04307
        momenta                    (time, I, a) float64 0.9542 0.4384 ... 1.21
        forces                     (time, I, a) float64 0.08174 0.08174 ... 0.237
        energy_kinetic             (time) float64 0.3102 0.2771 ... 0.3916 0.3169
        energy_potential           (time) float64 -6.299e+04 ... -6.299e+04
        stress                     (time, a, b) float64 -0.001528 ... -0.005702
        stress_kinetic             (time, a, b) float64 -0.001839 ... -0.002492
        stress_potential           (time, a, b) float64 0.0003113 ... -0.003209
        temperature                (time) float64 300.0 268.0 215.7 ... 378.7 306.4
        cell                       (time, a, b) float64 5.419 2.4e-11 ... 5.419
        positions_reference        (I, a) float64 2.98e-11 1.084e-11 ... 4.064 1.355
        lattice_reference          (a, b) float64 5.419 -2.4e-11 ... 4e-12 5.419
        force_constants_remapped   (Ia, Jb) float64 13.7 0.0 0.0 ... 9.099e-11 13.7
        forces_harmonic            (time, I, a) float64 0.0818 0.08181 ... 0.217
        energy_potential_harmonic  (time) float64 0.002015 0.02751 ... 0.3128 0.3856
        sigma_per_sample           (time) float64 0.006452 0.03701 ... 0.1663 0.1836
        pressure                   (time) float64 0.0009887 0.001096 ... nan 0.00426
        pressure_kinetic           (time) float64 0.0013 0.001161 ... 0.001328
        pressure_potential         (time) float64 -0.0003113 -6.477e-05 ... 0.002932
        aims_uuid                  (time) object '57870CD4-D00E-43EB-A6F6-D2C87EEDBD0C' ... 'FF5DB1AC-F1C5-433F-AE30-913BF062C6B6'
    Attributes:
        name:             trajectory
        system_name:      Si
        natoms:           8
        time_unit:        fs
        timestep:         4.000000000000006
        nsteps:           2500
        symbols:          ['Si', 'Si', 'Si', 'Si', 'Si', 'Si', 'Si', 'Si']
        masses:           [28.085 28.085 28.085 28.085 28.085 28.085 28.085 28.085]
        atoms_reference:  {"pbc": [true, true, true],\n"cell": \n[[ 5.41850551468...
        atoms_primitive:  {"pbc": [true, true, true],\n"cell": \n[[-1.00000000000...
        atoms_supercell:  {"pbc": [true, true, true],\n"cell": \n[[ 5.41850551468...
        volume:           159.08841208433154
        raw_metadata:     {"MD": {\n  "type": "molecular-dynamics",\n  "md-type":...
        hash:             2d33a63b08cd6018441fa85ece2ba97d357eefc8
        sigma:            0.1561093848741265
    ```

Each of the `Data variables` can be accessed as an attribute and converted to a plain `numpy` array by calling `.data`. To access e.g. the positions, we do

```python
positions_as_ndarray = dataset.positions.data
```

The arrays can be extracted and written to file e.g. via [`numpy.savetxt`](https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html) or by other means.

## Example: Analyze Pressure

As an example on how to perform postprocess directly from the `xarray` dataset, we will now evaluate the potential pressure observed during the simulation, as [introduced earlier](3_md_intro.md#example-pressure). We will use [xarray](http://xarray.pydata.org/) and [pandas](https://pandas.pydata.org/) for the analysis. For interactive data exploration, we recommend to run the code in a [jupyter notebook](https://jupyter.org/) and play around with the suggested parameters like windows sizes etc.

### Load trajectory dataset

We first load the trajectory dataset and visualize the temperature:

```python
import xarray as xr

# load the trajectory dataset "trajectory.nc" from disk into the xarray.Dataset
dataset  = xr.load_dataset("trajectory.nc")

# extract temperature and potential pressure from all the data and convert to pandas.DataFrame
df_temperature_pressure = all_data[["temperature", "pressure_potential"]].to_dataframe()

# attach a moving average (width=200 timesteps) of the temperature
df_temperature_pressure["temperature_mean"] = df_temperature_pressure.temperature.rolling(window=200).mean()

# plot temperature and temperature_mean as function of time
ax = df_temperature_pressure[["temperature", "temperature_mean"]].plot()

ax.set_xlabel("Time (fs)")
ax.set_ylabel("Temperature (K)")
```

??? info "`df_temperature_pressure.plot`"
	![image](assets/md_temperature.png)

Since the calculation starts with all atoms located at their equilibrium positions, the initial potential energy is zero and the kinetic energy corresponds to ~300K, since we have setup the velocities
using the Maxwell-Boltzmann distribution. In the first 500 steps of the trajectort, the kinetic energy is partially converted to potential energy at. In turn, the temperature drops from $300\,{\rm K}$ to about $150\,{\rm K}$. 
The missing thermal energy to obtaing a temperature of 300K is then gradually provided by the thermostat, bringing the nuclear temperature back to $\sim 300\,{\rm K}$ after a few $\rm ps$.

### Discard thermalization period
We can remove the thermalization period from the simulation data, e.g., by [shifting the dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.shift.html):

```python
# discard 500 steps (2ps) of thermalization
shift = 500

df_temperature_pressure = df_temperature_pressure.shift(-shift).dropna()

df_temperature_pressure[['temperature', 'temperature_mean']].plot()
```

??? info "`df_temperature_pressure.plot` after removing 2500 simulation steps ($5\,{\rm ps}$)"
	![image](assets/md_temperature_thermalized.png)

### Inspect the pressure
We are now ready to inspect the pressure observed in the simulation and plot it including its cumulative average:

```python
from ase.units import GPa

p = df_temperature_pressure.pressure_potential / GPa

ax = p.plot(alpha=0.75)

p.expanding().mean().plot(ax=ax, color="k")
```

??? info "Plot pressure"
	![image](assets/md_pressure.png)
	

### Expectation value and convergence estimation

[As discussed earlier](3_md_intro.md), the expectation value of the pressure is given by the mean of the observed pressures,

$$
\begin{align}
\left\langle p_{\rm Pot} \right\rangle
	= \lim_{N_{\rm t} \rightarrow \infty} \frac{1}{N_{\rm t}}
	\sum_n^{N_{\rm t}} 	
	p_{\rm Pot}({\bf R} (t_n))~.
\label{eq:<pPot>}
\end{align}
$$

In our finite simulation, $N_{\rm t} = 2000 < \infty$, so that

$$
\begin{align}
\left\langle p_{\rm Pot} \right\rangle
= \left\langle p_{\rm Pot} \right\rangle_{N_t = 2000} + \Delta~,
\label{eq:p_final}
\end{align}
$$

where $\left\langle p_{\rm Pot} \right\rangle_{N_t = 2000} = -0.076\,{\rm GPa}$ is the mean pressure observed during the finite simulation, and $\Delta$ is the (unknown) difference to the fully converged expectation value. 
Although, full converge would require an infinite trajetcory length and is thus formally never reachable, one can get arbitrarily close in practice and estimate the magnitude of the error $\Delta$.

We estimate this error by computing $\sigma_{\langle p \rangle}$, the [_standard error of the mean_](https://en.wikipedia.org/wiki/Standard_error):

$$
\begin{align}
\Delta \approx \sigma_{\langle p \rangle} = \frac{\sigma_p}{\sqrt{\tilde N_t}}~,
\label{eq:sigma_O}
\end{align}
$$

where $\sigma_p$ is the standard deviation of the pressure distribution observed during the simulation, and $\tilde N_t$ is an estimate of the number of _uncorrelated_ samples provided by the simulation. To this end, we estimate

$$
\begin{align}
\tilde N_t = N_t / \tau~,
\label{eq:N}
\end{align}
$$

where $\tau$ is the correlation time for the pressure.
The most straightforward way to compute $\tau$ is to evaluate the [autocorrelation function](https://en.wikipedia.org/wiki/Autocorrelation) and estimate its decay time:

```python
# estimate correlation time
import pandas as pd
from scipy import signal as si

# substract the mean pressure
pp = p - p.mean()

# get the autocorrelation function from
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html
corr = si.correlate(pp, pp)[len(pp) - 1 :]

# normalize to C(0) = 1
corr /= corr[0]

# create as pandas.Series for plotting
s = pd.Series(corr).rolling(min_periods=0, window=10).mean()
ax = s.plot()

# estimate correlation time from the drop below 0.1
tau = s.index.where(s < 0.1).min()
ax.axvline(tau)

ax.set_xlim(0, 100)
ax.set_title(f"$\\tau$ is {int(tau)} steps")
```

??? info "Plot pressure autocorrelation function"
	![image](assets/md_autocorr.png)

In the  present example, the observable decorrelates after about 10 time steps ($\equiv 40\,{\rm fs}$). We therefore estimate the number of uncorrelated samples to be 

$$
\begin{align*}
	\tilde N_t = N_t / 10 \approx 200
\end{align*}
$$

The standard deviation of the pressure distribution is

$$
\begin{align}
	\sigma_p = 0.239\,{\rm GPa}~,
\end{align}
$$

so that according to Eq. $\eqref{eq:sigma_O}$,

$$
\sigma_{\langle p \rangle} = \frac{0.239}{\sqrt{200}}\,{\rm GPa} \approx 0.053\,{\rm GPa}~.
$$


The final result for the pressure according to Eq. $\eqref{eq:p_final}$ is

$$
\begin{align*}
	\langle p_{\rm Pot} (300\,{\rm K}) \rangle = (-0.076 \pm 0.053)\,{\rm GPa}~,
\end{align*}
$$

which means that our result is converged within an estimated precision of $70\,\%$. **Remark:** This does _not_ mean that the true expectation lies within the given range. The estimated error is to be understood in the sense of a [confidence interval](https://en.wikipedia.org/wiki/Confidence_interval#Practical_example). The size of the error signals that the calculation is not fully converged and more sampling time would be necessary to report the observed pressure with confidence. You find reference for a total simulation time of $40\,{\rm ps}$ [here](https://gitlab.com/vibes-developers/vibes-tutorial-files/-/tree/master/3_molecular_dynamics/ab_initio/si_8_longer). How did the value and the error change?

Physics question: The observed potential pressure is negative. _Why?_ Do you expect a positive or negative lattice expansion at $300\,{\rm K}$?

??? info "Code snippet to compute the mean and the error estimator"

    ```python
    mean = p.mean()
    std = p.std()
    err = std / (len(p) / tau) ** 0.5
    
    print(f"Mean:  {mean:.5f} GPa")
    print(f"Std.:  {std:.5f} GPa")
    print(f"Error: {err:.5f} GPa ({abs(err / mean) * 100:.2f} %)")
    ```

### More examples

For more examples on how to directly work with the trajectory dataset in `trajectory.nc`, please have a look at  the [ASE Workshop Tutorial](https://gitlab.com/flokno/ase_workshop_tutorial_19) which analyzes _ab initio_ MD data for a perovskite.
