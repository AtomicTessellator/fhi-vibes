# Molecular Dynamics

!!! info
	We assume that you are familiar with the basics of molecular dynamics simulations and you are here to learn how to perform them with `FHI-vibes`.

## Recap

Molecular Dynamics (MD) simulations aim at exploring the dynamical properties of a system defined by the Hamiltonian $\mathcal H$

$$
\begin{align}
\mathcal{H}(\mathbf{R}, \mathbf{P})=\sum_{I} \frac{\mathbf{P}_{I}^{2}}{2 M_{I}}+\mathcal{V}(\mathbf{R})~,
\label{eq:H}
\end{align}
$$

where ${\bf R} = \{ {\bf R}_I \}$ denotes the atomic positions, ${\bf P} = \{ {\bf P}_I \}$ the atomic momenta, $M_I$ the atomic mass, and $\mathcal V ({\bf R})$ is a many-body potential, e.g., given by an empirical force-field, or an _ab initio_ energy functional. The dynamical evolution of the  system is described by the equations of motion,

$$
\begin{align}
M_{I} \ddot{\mathbf{R}}_{I}(t)=\mathbf{F}_{I}(t)=-\nabla_{I} \mathcal{V}(\mathbf{R}(t))~,
\label{eq:Newton}
\end{align}
$$

where the acceleration $\ddot{\mathbf{R}}_{I}(t)$ of an atom at time $t$ is given by the atomic force ${\bf F}_I (t)$, i.e., the inverse gradient of the many-body potential $\mathcal V ({\bf R})$.

Equation $\eqref{eq:Newton}$ is solved numerically from a given initial condition $\{{\bf R} (t_0), {\bf P} (t_0) \}$, for example with the [Velocity Verlet integrator](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet).

The dynamical information encoded in the phase-space trajectory $\Gamma (t) = \{{\bf R} (t), {\bf P} (t) \}$ can then be used to evaluate expectation values of observables, which are typically functions of phase-space points:

$$
\begin{align}
	\langle O\rangle
	%&= \frac{1}{\mathcal{Z}} 
	%\int \mathrm{d} \mathbf{R} \mathrm{d} \mathbf{P} ~ 
	%	\mathrm{e}^{-\beta \mathcal{H}(\mathbf{R}, {\bf P})} 
	%	O(\mathbf{R}, {\bf P}) 
	%	\label{eq:O1} \\
	&=
	\lim_{T \rightarrow \infty} \frac{1}{T}
	\int_0^T {\rm d} t ~
	O(\mathbf{R} (t), {\bf P} (t))
	\label{eq:O2}~,
	\end{align}
$$

where Eq. $\eqref{eq:O2}$ holds for [ergodic systems](https://en.wikipedia.org/wiki/Ergodicity).

### Thermostats
A simulation as described above models a [_microcanonical ensemble_](https://en.wikipedia.org/wiki/Microcanonical_ensemble), where particle number $N$, volume $V$, and energy $E$ are conserved. 

In order to simulate a [_canonical ensemble_](https://en.wikipedia.org/wiki/Canonical_ensemble), where instead of the energy $E$ the temperature $T$ is the thermodynamic variable, one typically models the system including an interaction with a fictitious heat bath which allows to control the target temperature of the  system.

#### Langevin thermostat
`FHI-vibes` uses a Langevin thermostat for canonical sampling. In Langevin dynamics, modified equations of motion are used with

$$
\begin{align}
	\dot{\bf P}_I(t)
		= {\bf F}_I 
		- \gamma {\bf P}_I(t) 
		+ \sqrt{2 M_I \gamma T} \xi(t)
	\label{eq:Langevin}~,
\end{align}
$$

where $\gamma$ is a friction parameter and $\xi (t)$ is a white-noise term obeying $\langle\xi(t) \xi(0)\rangle=2 k_{\rm B} T \gamma \delta(t)$.