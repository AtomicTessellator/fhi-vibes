# Anharmonicity Quantification

!!! warning
	The tutorial assumes you are familiar with performing [phonon calculations](2_phonopy.md) and [molecular dynamics simulations](3_md_ab_initio.md).

## Background

As detailed [in our paper](https://arxiv.org/abs/2006.14672), we define the anharmonic contribution to the potential energy $\mathcal V ({\bf R})$ as

$$
\begin{align}
	\mathcal{V}^{\rm A}(\mathbf{R}) \equiv \mathcal{V}(\mathbf{R})-\mathcal{V}^{(2)}(\mathbf{R})~,
	\label{eq:VA}
\end{align}
$$

where $\mathcal{V}^{(2)}(\mathbf{R})$ at a given atomic configuration $\bf R$ is given by

$$
\begin{align}
	\mathcal{V}^{(2)}\left(\mathbf{R}=\mathbf{R}^{0}+\Delta \mathbf{R}\right)
	=\frac{1}{2} \sum_{I, J} \Phi_{\alpha \beta}^{I, J} \Delta R_{I}^{\alpha} \Delta R_{J}^{\beta}~,
\label{eq:V2}
\end{align}
$$

with the [harmonic force constants $\Phi^{IJ}$](2_phonopy_intro.md) obtained at the equilibrium configuration ${\bf R}^0$ as

$$
\begin{align}
	\Phi_{\alpha, \beta}^{I, J}
	=\left.\frac{\partial^{2} \mathcal{V}}{\partial R_{I}^{\alpha} \partial R_{J}^{\beta}}\right|_{\mathbf{R}^{0}}~.
	\label{eq:Phi}
\end{align}
$$

Likewise, we define the anharmonic contribution to the force components $F_{I, \alpha} ({\bf R})$ as

$$
\begin{align}
	F_{I, \alpha}^{\mathrm{A}}(\mathbf{R})
	=F_{I, \alpha}(\mathbf{R})-F_{t, \alpha}^{(2)}(\mathbf{R})~.
	\label{eq:FA}
\end{align}
$$

This is a depiction of Eq. $\eqref{eq:VA}$ and $\eqref{eq:FA}$ for a one-dimensional toy potential:

![image](assets/PES_sketch.png)

In order to estimate the strength of anharmonic effects in a material, we define the _anharmonicity measure_

$$
\begin{align}
\sigma^{\mathrm{A}}(T) \equiv \frac{\sigma\left[F^{\mathrm{A}}\right]_{T}}{\sigma[F]_{T}}=\sqrt{\frac{\sum_{I, \alpha}\left\langle\left(F_{I, \alpha}^{\mathrm{A}}\right)^{2}\right\rangle_{T}}{\sum_{I, \alpha}\left\langle\left(F_{I, \alpha}\right)^{2}\right\rangle_{T}}}~,
\label{eq:sigmaA}
\end{align}
$$

where $\langle \cdot \rangle_T$ denotes an [expectation value at a given temperature](3_md_postprocess.md#expectation-value-and-convergence-estimation), 

$$
\begin{align}
	\left\langle O \right\rangle
	= \lim _{N_{\mathrm{t}} \rightarrow \infty} 
	\frac{1}{N_{\mathrm{t}}} \sum_{n}^{N_{\mathrm{t}}} \left(t_{n}\right)~.
	\label{eq:meanO}
\end{align}
$$

$F_{I, \alpha} (t) \equiv F_{I, \alpha} [{\bf R} (t)]$ is the force component $\alpha$ on atom $I$ at time $t$, and $F^{\rm A}_{I, \alpha}$ is given by Eq. $\eqref{eq:FA}$. $\sigma^{\rm A} (T)$ therefore quantifies the _average strength of anharmonic force components $F_{I, \alpha}^{\rm A}$, normalized by the average strength of forces $F_{I, \alpha}$, observed at temperature $T$_.



## Evaluating anharmonicity with `FHI-vibes`

The necessary ingredient to evaluate Eq. $\eqref{eq:sigmaA}$ are:

- Atomic forces $F_{I, \alpha}$,
- harmonic force constants $\Phi^{IJ}$ to compute $F_{I, \alpha}^{(2)}$ according to Eq. $\eqref{eq:F2}$ for evaluating Eq. $\eqref{eq:FA}$ , and
- thermodynamic expectation values according to Eq. $\eqref{eq:meanO}$.

These ingredients can be obtained with `FHI-vibes` with the following workflow:

- Take the materials of interest and generate a reference structure, i.e., a primitive cell and a supercell.
- Obtain force constants for the supercell as introduced in the [phonons tutorial](2_phonopy.md).
- Run an MD simulation for the supercell as introduced in the [MD tutorial](3_md_canonical_sampling.md).

### Example: LDA-Silicon at room temperature

Assuming that  you performed the previous tutorials for LDA-Silicon in an 8-atoms supercell, we already have all the necessary ingredients available to evaluate $\sigma^{\rm A}$ for this system!

In a new working directory, copy your `trajectory.nc` dataset from the the [MD tutorial](3_md_ab_initio.md) and your force constants from the [phonopy tutorial](2_phonopy.md), i.e., the file `phonopy/output/FORCE_CONSTANTS`. You can attach the force constants to the trajectory dataset with the CLI tool `utils trajectory update`:

```
vibes utils trajectory update trajectory.nc -fc FORCE_CONSTANTS 
```

This will attach read the force constants from `FORCE_CONSTANTS` and attach them to the trajectory dataset. 

To evaluate Eq. $\eqref{eq:sigmaA}$, you can use the CLI tool `utils anharmonicity sigma`:

```
vibes utils anharmonicity sigma trajectory.nc
```

which will give you the total $\sigma^{\rm A}$ value (`sigma`), as well as an individual value for each atom species. The output should be 

```
DataFrame:
       sigma  sigma [Si]  sigma_atom_mean  sigma_mode
Si  0.156109    0.156109         0.156109    0.156026
```

This tells you that the average magnitude of anharmonic contributions to the forces, $F^{\rm A}$, in LDA-Silicon at $300\,{\rm K}$ is about $16\,\%$.

## Mode resolved anharmonicity

To obtain a mode-resolved $\sigma^{\rm A}_s$ similar to the analysis of Fig. 8 in [our paper](https://arxiv.org/pdf/2006.14672.pdf), you can run

```
vibes utils anharmonicity mode trajectory.nc 
```

which will produce a `.csv` file containing mode frequencies $\omega_s$ in THz and the respective mode-resolved anharmonicity $\sigma^{\rm A}_s$.

You can plot the file e.g. via

```python
import pandas as pd

s = pd.read_csv("sigmaA_mode_Si.csv", index_col=0)

ax = s.plot(marker=".", lw=0)

ax.set_xlim(0, 20)
ax.set_ylim(0, 0.5)

ax.set_xlabel("$\omega_s$ (THz)")
ax.set_ylabel(r"$\sigma^{\rm A}_s$")
```

??? info "Plot of $\sigma^{\rm A}_s$"
	![image](assets/sigma_mode_Si.png)
	
The plot won't look too impressive because we're using a small supercell and the anharmonicity in silicon is overall quite weak. But you should be good to go to investigate the anharmonicity of your material of choice by now -- Happy Computing 💪