# Phonons

## Outline

During the course of this workshop, we have used periodic boundary conditions 
in first-principles calculations to efficiently describe macroscopic, crystalline materials. It is important
to realize that the application of periodic boundary conditions relies on the 
assumption that the nuclei constitute an immobile grid with fixed periodicity. 
However, thermodynamic fluctuations constantly lead to displacements from this
perfectly periodic grid of equilibrium positions --  even at zero temperature 
due to the quantum mechanical zero point motion. Accounting for this dynamics 
is essential to understand the physics of many fundamental material properties
such as the specific heat, the thermal expansion, as well as charge and heat
transport.

To introduce you to these effects, this tutorial consists of two parts:

1. **Phonons: Harmonic Vibrations in Solids**
   - Using `phonopy` via `vibes`
   - Supercell size convergence
   - Lattice Expansion: The Quasi-Harmonic Approximation
2. **Electron-Phonon Coupling: Band Gap Renormalization**
   - The Role of the Lattice Expansion
   - The Role of the Atomic Motion

In part I, we will compute the vibrational properties
of a solid using the _harmonic approximation_. In particular, we will discuss and 
investigate the convergence with respect to the supercell size used in the calculations.
Furthermore, we will learn how the harmonic approximation can be extended in a straightforward
fashion to approximatively account for a certain degree of anharmonic effects~(quasi-harmonic
approximation) and how this technique can be used to compute the thermal lattice expansion.

In [part 1](#Phonons), we will then go back to electronic structure theory and investigate
how the fact that the nuclei are not immobile affects the electronic band structure. Both the role
of the lattice expansion and of the atomic motion will be discussed and analyzed. 

## <a name="Phonons"></a> Phonons: Harmonic vibrations in solids

To determine the vibrations in a solid, we approximate the potential energy surface
for the nuclei by performing a Taylor expansion of the total energy $E$ around the equilibrium positions:

$$
\require{cancel}
\begin{aligned}
\def\vec#1{{\mathbf{#1}}}
\def\t#1{\text{#1}}
E \left(\{\vec{R}^0 + \Delta \vec{R}\}\right) 
& \approx 
E\left(\{\vec{R}^0\}\right) \\
& + \cancel{ \sum\limits_{I} \left.\frac{\partial E}{\partial \vec{R}_I}\right\vert_{\vec{R}^0} \Delta\vec{R}_{I} } \\
& + \frac{1}{2} \sum\limits_{I,J} \left.\frac{\partial^2 E}{\partial \vec{R}_I\partial \vec{R}_J}\right\vert_{\vec{R}^0} \Delta\vec{R}_{I}\Delta\vec{R}_{J} \\
& + \mathcal{O}(\Delta\vec{R}^3)
\end{aligned}
$$

The linear term vanishes, since no forces $\vec{F} = - \nabla E$ are acting on the system in equilibrium $\vec{R}^0$. 
Assessing the Hessian $\Phi_{IJ} = \frac{\partial^2 E}{\partial \vec{R}_I\partial \vec{R}_J}$ involves some additional
complications: In contrast to the forces $\vec{F}$, which only depend on the density, the Hessian $\Phi_{IJ}$ also depends
on its derivative with respect to the nuclear coordinates, i.e., on its _response_ to nuclear displacements. One can either 
use _Density Functional Perturbation Theory (DFPT)_ [[Baroni2001](references.md#baroni2001)] to compute the response 
or one can circumvent this problem by performing the second order derivative _numerically by finite differences_

$$
\begin{align}
\Phi_{IJ} 
= \left.\frac{\partial^2 E}{\partial \vec{R}_I\partial \vec{R}_J}\right\vert_{\vec{R}^0} 
= - \left.\frac{\partial }{\partial \vec{R}_I} \vec{F}_J\right\vert_{\vec{R}^0}
\approx - \frac{ \vec{F}_J(\vec{R}_I^0 + \varepsilon \,\vec{d}_I)}{\varepsilon}~,
\label{eq:FinDiff}
\end{align}
$$

as we will do in this tutorial.
The definition in Eq.$~\eqref{eq:FinDiff}$ is helpful to realize that the Hessian describes a coupling between different atoms, i.e., how the force acting on an atom $\vec{R}_J$ changes
if we displace atom $\vec{R}_I$, as you have already learned in tutorial 1. However, an additional complexity arises in the case of _periodic boundary conditions_,
since beside the atoms in the unit cell $\vec{R}_J$ we also need to account for the periodic images $\vec{R}_{J'}$. Accordingly, the Hessian is in principle a matrix 
of infinite size. In non-ionic crystals, however, the interaction between two atoms$~I$ and $J$ quickly decays with their distance$~\vec{R}_{IJ}$, so that we can compute the Hessian from
finite supercells, the size convergence of which must be accurately inspected (cf. Exercise \hyperref[ex2]{2}).

Once the real-space representation of the Hessian is computed, we can determine the _dynamical matrix_ by adding up the contributions
from all periodic images$~J'$ in the mass-scaled Fourier transform of the Hessian:

$$
\begin{align}
D_{IJ}({\vec{q}}) = \sum\limits_{J'} 
\frac{\t e^{\t{im} \left({\vec{q}}\cdot{\vec{R}_{JJ'}}\right)}}{\sqrt{M_I M_J}} 
\;\Phi_{IJ'} 
\quad .
\label{DynMat}
\end{align}
$$

In reciprocal space [[AshcroftMermin](references.md#AshcroftMermin)],
this _dynamical matrix_ determines the equation of motion for such a periodic array of harmonic 
atoms for each reciprocal vector$~\vec{q}$:

$$
\begin{align}
D(\vec{q}) \, \vec e_s (\vec{q}) = \omega_s^2(\vec{q}) \, \vec e_s (\vec{q})
\; .
\label{eq:eigenproblem}
\end{align}
$$

The dynamical matrix has dimension $3N_\t{A} \times 3N_\t{A}$, where $N_\t{A}$ 
is the number of 
atoms in the *primitive* unit cell. Equation$~\eqref{eq:eigenproblem}$ thus 
constitutes 
an eigenvalue problem with $3N_\t{A}$ solutions at each $\vec q$ point. The 
solutions 
are labelled by $s$ and are denoted as *phonon branches*. The lowest three 
branches are commonly called the *acoustic branches*, whereas in solids 
with more than one atom in the primitive unit cell, the remaining $(3N_\t{A} - 
3)$ branches are denoted as *optical branches*.

The eigenvalues$~\omega_s^2(\vec{q})$ and eigenvectors$~\vec e_s(\vec{q})$ of the 
dynamical matrix$~D(\vec{q})$
completely describe the dynamics of the system (in the harmonic approximation), which is nothing else than a superposition
of harmonic oscillators, one for each mode, i.e., for each 
eigenvalue$~\omega_s (\vec{q})$. 

From the set of eigenvalues $\{ \omega_s (\vec{q}) \}$, also denoted as 
spectrum, the *density of states* (DOS) can be obtained by 
calculating the number of states in an infinitesimal energy window 
$[\omega, \omega + \t d \omega]$:

$$
\begin{align}
g(\omega) = \sum_s\int\frac{\t d \vec{q}}{(2\pi)^3}\delta(\omega - 
\omega(\vec{q})) = \sum_s\int\limits_{\omega(\vec{q}) = 
\omega}\frac{\t{d}S}{(2\pi)^3}\frac{1}{\vert\nabla\omega(\vec{q})\vert}~.
\label{DOS}
\end{align}
$$

The DOS is a very useful quantity, since it allows to determine any integrals 
(the integrand of which only depends on $\omega$) by a 
simple integration over a one-dimensional variable$~\omega$ rather than a three-dimensional variable$~\vec{q}$. This is much 
easier to handle both in numerical and in analytical models. For instance, we can compute the associated thermodynamic 
potential[^footnote1], i.e., the (harmonic) *Helmholtz free 
energy*[^footnote2]

$$
\begin{equation}
F^{\mathrm{ha}}(T,V)  = \int \t d \omega\; g(\omega) \left ( 
\frac{\hbar\omega}{2} + k_\t{B}\, T~ 
\ln\left(1-\t{e}^{-\frac{\hbar\omega}{k_\t{B}\,T}}\right) 
\right)\;.
\label{HFVT}
\end{equation}
$$

In turn, this allows to calculate the heat capacity at constant 
volume

$$
\begin{equation}
C_V = - T \left(\frac{\partial^2 F^{\mathrm{ha}}(T,V)}{\partial T^2} \right)_V \;.
\label{CV}
\end{equation}
$$

For a comprehensive introduction to the field of lattice dynamics and its 
foundations, we refer to [[BornHuang](references.md#BornHuang)].

To compute the quantities introduced above, we will use 
`FHI-vibes`, which uses the package _phonopy_ [[Togo2015](references.md#Togo2015)] as a backend to compute vibrational properties via the finite-displacements method as outlined above. Please note that 
_phonopy_ makes extensive use of symmetry 
analysis [[Parlinski1997](references.md#Parlinski1997)], which allows to reduce numerical noise and to speed up the calculations considerably. Out system of choice will be fcc-diamond Silicon.

!!! warning
    In the following exercises, the computational settings, in particular the reciprocal space grid (tag `k_grid`), the basisset and supercell sizes, have been chosen to allow a rapid computation of the exercises. In a _real_ production calculation, the reciprocal space grid, the basis set, and the supercells would all have to be converged with much more care,  although the qualitative trends hold already with the present settings.

[^footnote1]: Given that the _Bose-Einstein distribution_ is used for
the derivation of the harmonic free energy in this case, we get the correct quantum-mechanical result including zero-point effects by this means.

[^footnote2]: A derivation of Eq.$\,\eqref{HFVT}$ will be presented in Tutorial 7 on 
Molecular Dynamics Simulations.

