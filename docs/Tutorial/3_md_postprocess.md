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
	
#### Expectation value and convergence estimation

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

In our finite simulation, $N_{\rm t} = 5000 < \infty$, so that

$$
\begin{align}
\left\langle p_{\rm Pot} \right\rangle
= \left\langle p_{\rm Pot} \right\rangle_{N_t = 5000} + \Delta~,
\label{eq:p_final}
\end{align}
$$

where $\left\langle p_{\rm Pot} \right\rangle_{N_t = 5000} = 0.1064\,{\rm GPa}$ is the mean pressure observed during the finite simulation, and $\Delta$ is the (unknown) difference to the fully converged expectation value. While we are never able to fully converge a calculation, we can nevertheless estimate the magnitude of the error $\Delta$.

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

In the  present example, the observable decorrelates after about 192 time steps ($\equiv 768\,{\rm fs}$). We therefore estimate the number of uncorrelated samples to be 

$$
\begin{align*}
	\tilde N_t = N_t / 192 \approx 26
\end{align*}
$$

The standard deviation of the pressure distribution is

$$
\begin{align}
	\sigma_p = 0.0046\,{\rm GPa}~,
\end{align}
$$

so that according to Eq. $\eqref{eq:sigma_O}$,

$$
\sigma_{\langle p \rangle} = \frac{0.0046}{\sqrt{26}}\,{\rm GPa} \approx 0.0009\,{\rm GPa}~.
$$


The final result for the pressure according to Eq. $\eqref{eq:p_final}$ is

$$
\begin{align*}
	\langle p_{\rm Pot} (20\,{\rm K}) \rangle = (0.1064 \pm 0.0009)\,{\rm GPa}~,
\end{align*}
$$

which means that our result is converged within an estimated precision of $1\,\%$. **Remark:** This does _not_ mean that the true expectation lies within the given range. The estimated error is to be understood in the sense of a [confidence interval](https://en.wikipedia.org/wiki/Confidence_interval#Practical_example).

??? info "Code snippet to compute the mean and the error estimator"

    ```python
    mean = p.mean()
    std = p.std()
    err = std / (len(p) / tau) ** 0.5

    print(f"Mean:  {mean:.5f} GPa")
    print(f"Std.:  {std:.5f} GPa")
    print(f"Error: {err:.5f} GPa ({err / mean * 100:.2f} %)")
    ```

### More examples

For more examples on how to directly work with the trajectory dataset in `trajectory.nc`, please have a look at  the [ASE Workshop Tutorial](https://gitlab.com/flokno/ase_workshop_tutorial_19) which analyzes _ab initio_ MD data for a perovskite.