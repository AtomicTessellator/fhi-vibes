# Forceconstant Calculator

The force constant calculator `FCCalculator` implements properties derived from a harmonic model for the potential-energy surface,
$$
H ({\bf U}) = \frac{1}{2} \sum_{IJ, \alpha \beta} \Phi^{\alpha \beta}_{IJ} \, U_I^\alpha U_J^\beta~,\label{eq:H}
$$
where $\Phi^{\alpha \beta}_{IJ}$ denote the force constant matrix elements 
$$
\Phi^{\alpha \beta}_{IJ} 
= \left. \frac{\partial^2 E}{\partial R_I^\alpha \partial R_J^\beta} \right\vert_{{\bf R}^0}~,\label{eq:Phi}
$$


i.e., the second-order derivative of the potential energy $E$ w.r.t. to the atomic positions ${\bf R} = \{R_I^\alpha\}$ at a given configuration ${\bf R}^0$. Here, $I,J$ are the atom labels and $\alpha,\beta$ denote Cartesian components. ${\bf U} = \{U_I^\alpha\}$ denotes the atomic displacements w.r.t. to the reference position ${\bf R}^0$, i.e., 
$$
{\bf U} = {\bf R} - {\bf R}^0~.\label{eq:U}
$$

## Properties derived from the harmonic model

### Forces

Forces are given as a derivative of the potential energy and can therefore be obtained analytically from Eq.$\,\eqref{eq:H}$:
$$
{F_I^\alpha} = - \frac{\partial H}{\partial R_I^\alpha} = -\sum_{J, \beta} \Phi^{\alpha \beta}_{IJ} \, U_J^\beta~.\label{eq:F_I}
$$
This equations can be written as a matrix-vector product
$$
{\bf F} = - P {\bf U}~,\label{eq:F}
$$
where ${\bf F}$ and $\bf U$ are the $3N\times 1$ vectors representing atomic forces and displacements, and $P$ is the $3N \times 3N$ forceconstant matrix, where $N$ is the number of atoms in the simulation. In python, this is achieved by

```python
forces = -(fc @ displacements.flatten()).reshape(displacements.shape)
```

### Potential energy

The energy is given by computing the displacements $\bf U$ according to Eq.$\,\eqref{eq:U}$, and performing the matrix product given by Eq.$\,\eqref{eq:H}$, i.e., by scalar multiplication of the displacements $\bf U$ with the forces $\bf F$,
$$
E^{\rm pot} = - \frac{1}{2} {\bf F} \cdot {\bf U}~.
$$
In python:

```python
# energies: [N, 3] * [N, 3] -> [N]
energies = -(displacements * forces).sum(axis=1) / 2
# energy: [N] -> [1]
energy = energies.sum()
```

As a byproduct, atomic energy contributions (`energies`) are computed, given by
$$
E^{\rm pot}_I 
	= \frac{1}{2} \sum_{J, \alpha \beta} \Phi^{\alpha \beta}_{IJ} \, U_I^\alpha U_J^\beta~.
$$

### Stress

Stress is given in terms of a strain derivative of the potential energy
$$
\sigma^{\alpha \beta} = \frac{1}{V} \frac{\partial H}{\partial \epsilon^{\alpha\beta}}~,\label{eq:sigma}
$$
where $\epsilon^{\alpha \beta}$ describes a rotation-free, homogeneous straining of the atomic positions. The strain derivatives of the atomic displacements read [Knuth2015]
$$
\frac{\partial U_I^\gamma}{\partial \epsilon^{\alpha \beta}} = \delta^{\alpha \gamma} U_I^\beta~.\label{eq:dUde}
$$
In turn, the harmonic stress tensor reads
$$
\sigma^{\alpha \beta} 
	= \frac{1}{V} \frac{\partial H}{\partial \epsilon^{\alpha \beta}}
	= \underset{\sigma^{\alpha \beta}_{\rm HA}}{\underbrace{\frac{1}{V} \sum_{IJ, \gamma} \Phi^{\alpha \gamma}_{IJ} \, U_I^\beta U_J^\gamma}}
	+ \sigma_{\rm QHA}^{\alpha \beta}~,
$$
where the purely harmonic contribution is given by
$$
\sigma_{\rm HA}^{\alpha \beta} = \frac{1}{V} \sum_{IJ, \gamma} \Phi^{\alpha \gamma}_{IJ} \, U_I^\beta U_J^\gamma ~,
% \equiv \frac{2}{V} \, H~
\label{eq:sigma_ha}
$$
and the the _quasiharmonic_ contribution $\sigma_{\rm QHA}^{\alpha \beta}$ stems from the strain derivative of the force constants. The quasiharmonic contribution can be approximated by
$$
\frac{1}{V} \frac{\partial H}{\partial \epsilon^{\alpha \beta}} \approx 3 \left( \frac{\partial H}{\partial V} \right)^{\alpha \beta}~,
$$
which is strictly true only in cubic systems. By denoting the quasiharmonic force constants, i.e., the volume derivative $\partial \Phi / \partial V$ as $\Phi'$, we have
$$
\sigma_{\rm QHA}^{\alpha \beta} = \frac{3}{2} \sum_{IJ, \gamma} \Phi'^{\alpha \gamma}_{IJ} \, U_I^\beta U_J^\gamma~.
$$

### Virial stress

The virial stress in periodic systems is given by [Louwerse2006]
$$
\sigma_{\rm vir}^{\alpha \beta} 
	= - \frac{1}{2V} \sum_{IJ} (R_I^\alpha -R_J^\alpha) F_{IJ}^\beta~,
\label{eq:stress_vir}
$$
where $F_{IJ}^\beta$ denotes the pairwise force between atom $I$ and $J$ fulfilling Newton's 3rd law,
$$
F_{IJ}^\alpha = - F_{JI}^\beta~,
\label{eq:Newton3}
$$
and sums to the atomic force $F_I^\beta$,
$$
F_I^\beta = \sum_J F_{IJ}^\beta~.
\label{eq:F_I_sum}
$$
This pairwise force is given in the harmonic approximation as 
$$
F_{IJ}^\beta 
	= - \Phi_{IJ}^{\beta \gamma} U_J^\gamma + \Phi_{JI}^{\beta \gamma} U_I^\gamma~,
\label{eq:F_IJ}
$$
where summing over repeated indices is implied. This definition complies with both Eq.$\,\eqref{eq:Newton3}$ and $\eqref{eq:F_I_sum}$, by using the acoustic sum rule $\sum_J \Phi_{IJ}^{\alpha \beta} = 0$. Using this pairwise force in the definition of virial stress, Eq.$\,\eqref{eq:stress_vir}$, we have
$$
\begin{align}
\sigma_{\rm vir}^{\alpha \beta}
	=& - \frac{1}{V} \sum_{IJ} (R_I^\alpha - R_J^\alpha) \, \Phi_{IJ}^{\beta \gamma} U_J^\gamma \label{eq:stress_vir_ha} \\
	=& - \frac{1}{V} \sum_{IJ} (R_I^{0 \alpha} - R_J^{0 \alpha}) \, \Phi_{IJ}^{\beta \gamma} U_J^\gamma \nonumber \\
	& - \frac{1}{V} \sum_{IJ} (U_I^{\alpha} - U_J^{\alpha}) \, \Phi_{IJ}^{\beta \gamma} U_J^\gamma~,
\end{align}
$$
where in the second row ${\bf R} = {\bf R}^0 + {\bf U}$ was written out. However, the first contribution can be written as
$$
V \sigma_{\rm vir, 0}^{\alpha \beta}
	= - \sum_{IJ} (R_I^{0 \alpha} - R_J^{0 \alpha}) \, \Phi_{IJ}^{\beta \gamma} U_J^\gamma
	= \underset{\text{acoustic sum rule}}{\underbrace{\cancel{\sum_{IJ} R_J^{0 \alpha} \, \Phi_{IJ}^{\beta \gamma} U_J^\gamma}}}
	- \sum_{IJ} R_I^{0 \alpha} \, \Phi_{IJ}^{\beta \gamma} U_J^\gamma
	\equiv \sum_I R_I^{0 \alpha} F_I^\beta~,
$$
which does not contribute to thermal averages because $\langle F \rangle = 0$.

#### Remark on finite supercells

In finite supercells, care must be taken when computing the pair distance ${\bf R}_I - {\bf R}_J$. They should obey a minimum image convention (MIC). Since MIC is a non-continuous map as it involves wrapping of atoms by finite displacements, continuity of the map must be ensured, e.g., by averaging over MIC-equivalent pairs of atoms.

### Atomic virials

The harmonic virial stress can be decomposed into atomic contributions by writing
$$
\sigma_{\rm vir}^{\alpha \beta} = \sum_I \sigma_{{\rm vir}, I}^{\alpha \beta}~,
$$
with
$$
\sigma_{{\rm vir}, I}^{\alpha \beta}
	= - \frac{1}{V} \sum_{IJ} (R_I^\alpha - R_J^\alpha) \, \Phi_{IJ}^{\beta \gamma} U_J^\gamma \label{eq:stress_vir_I}~.
$$
This quantity is related to the _virial heat flux_ ${\bf J}_{\rm vir}$by multiplying with the atomic velocities $\dot{\bf U}$,
$$
{\bf J}_{\rm vir}^\alpha = \sum_{I, \gamma} \sigma_{{\rm vir}, I}^{\alpha \gamma} \dot{U}_I^\gamma~.
$$




## Literature

- `Knuth2015`: F. Knuth, C. Carbogno, V. Atalla, V. Blum, and M. Scheffler, Comput. Phys. Commun. **190**, 33 (2015).
- `Louwerse2006`: M. J. Louwerse and E. J. Baerends, Chem. Phys. Lett. **421**, 138 (2006).

