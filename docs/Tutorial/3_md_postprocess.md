# Postprocessing

!!! warning
	We assume you successfully ran the simulation of Lennard-Jones Argon at $20\,{\rm K}$ from the [previous chapter](3_md_canonical_sampling.md) and created a `trajectory.nc` dataset.



### Example: Pressure

We will now evaluate the potential pressure observed during the simulation, as [introduced earlier](3_md_intro.md#example-pressure). We will use [xarray](http://xarray.pydata.org/) and [pandas](https://pandas.pydata.org/) for the analysis. For interactive data exploration, we recommend to run the code in a [jupyter notebook](https://jupyter.org/) and play around with the suggested parameters like windows sizes etc.

#### Load trajectory dataset

We first load the trajectory dataset and visualize the temperature:

```python
import xarray as xr


# load the trajectory dataset
ds = xr.load_dataset("trajectory.nc")

# extract temperature and potential pressure and convert to pandas.DataFrame
df = ds[["temperature", "pressure_potential"]].to_dataframe()

# attach a moving average of the temperature
df["temperature_mean"] = df.temperature.rolling(window=200).mean()

# plot the dataframe
df[["temperature", "temperature_mean"]].plot()
```

??? info "`df.plot`"
	![image](../assets/md_temperature.png)
	
Since the calculation starts with all atoms located at their equilibrium positions, the initial potential energy is zero and the kinetic energy given to the system is converted to potential energy at early simulation times. In turn, the temperature drops from $20\,{\rm K}$ to about $10\,{\rm K}$. The missing energy is gradually provided by the thermostat, bringing the nuclear temperature back to $\sim 20\,{\rm K}$ after a few $\rm ps$.

#### Discard thermalization period
We can remove the thermalization period from the simulation data, e.g., by [shifting the dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.shift.html):

```python
# discard 2500 steps (5ps) of thermalization
shift = 2500

df = df.shift(-shift).dropna()

df[['temperature', 'temperature_mean']].plot()
```

??? info "`df.plot` after removing 2500 simulation steps ($5\,{\rm ps}$)"
	![image](../assets/md_temperature_thermalized.png)

#### Inspect the pressure
We are now ready to inspect the pressure observed in the simulation and plot it with a cumulative average:

```python
from ase.units import GPa

p = df.pressure_potential / GPa

ax = p.plot(alpha=0.75)

p.expanding().mean().plot(ax=ax, color="k")
```

??? info "Plot pressure"
	![image](../assets/md_pressure.png)
	
#### Convergence estimation
According to the [central-limit theorem](https://en.wikipedia.org/wiki/Central_limit_theorem), the standard deviation of the expectation value $\sigma_{\langle O \rangle}$ is given by the standard deviation of the sampled points $\sigma_O$ divided by the square root of the number of _uncorrelated_ samples, $\tilde N_t$:

$$
\begin{align}
\sigma_{\langle O \rangle} = \frac{\sigma_O}{\sqrt{\tilde N_t}}~.
\label{eq:sigma_O}
\end{align}
$$

Obviously, consecutive samples during an MD simulation are correlated. It is thereforce _not correct_ to compute the mean of the observable for all datapoints and give an error bar by computing the standard deviation divided by $\sqrt{N_t}$. We also need to estimate a correlation time $\tau$ for the observable.

The most straightforward way is to compute the correlation time by evaluating the [autocorrelation function](https://en.wikipedia.org/wiki/Autocorrelation) and estimate its decay time:

```python
# estimate correlation time
import pandas as pd
from scipy import signal as si

# substract the mean pressure
pp = p - p.mean()

# get the autocorrelation function from
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html
corr = si.correlate(pp, pp)[len(s) - 1 :]

# normalize to C(0) = 1
corr /= corr[0]

# create as pandas.Series for plotting
s = pd.Series(corr)
ax = s.plot()

# estimate correlation time from the drop below 0.1
tau = s.index.where(s < 0.1).min()
ax.axvline(tau)

ax.set_xlim(0, 500)
ax.set_title(f"$\\tau$ is {int(tau)} steps")
```

??? info "Plot pressure"
	![image](../assets/md_autocorr.png)

In the  present example, the observable decorrelates after about 192 time steps ($\equiv 768\,{\rm fs}$). We therefore estimate the number of uncorrelated sample to be 

$$
\begin{align*}
	\tilde N_t = N_t / 192
\end{align*}
$$

and

```python
mean = p.mean()
err = p.std() / (len(p) / tau) ** 0.5

print(f"Mean:  {mean:.5f} +/- {err:.5f} GPa")
print(f"Error: {err / mean * 100:.2f} %")
```

gives a final results of 

$$
\begin{align*}
	\langle p_{\rm Pot} \rangle = (0.1064 \pm 0.0009)\,{\rm GPa}~.
\end{align*}
$$


### More examples

For more examples on how to directly work with the trajectory dataset in `trajectory.nc`, please have a look at  the [ASE Workshop Tutorial](https://gitlab.com/flokno/ase_workshop_tutorial_19) which analyzes _ab initio_ MD data for a perovskite.