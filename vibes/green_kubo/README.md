Green Kubo
===

## Thermal Conductivity

$$
\kappa^{\alpha \beta} = \frac{V}{k_\text{B} T^2} \lim_{t \to \infty} \int\limits_0^t \left\langle J^\alpha (\tau) J^\beta \right\rangle ~ \text{d} \tau
$$

## Interpolation

### Notation

- Mode index $s = (b, {\bf q})$, $-s = (b, -{\bf q})$, where $b$ denotes band index and $\bf q$ the wave vector.
- $u_s = \sum_I {\bf e}_{Is} \cdot u_I$: Mode displacement
- $p_s = \sum_I {\bf e}_{Is} \cdot p_I$: Mode momentum

- Complex amplitudes [BornHuang (38.38)]
    $$
    \begin{align}
    \def\D{\dagger}
    u_s &= a^\D_{-s} + a_s \\
    p_s &= {\rm i}\omega_s \left( a^\D_{-s} - a_s \right)~,
    \end{align}
    $$
    and
    $$
    \begin{align}
    a_s &= ... \\
    a^\D_{-s} &= ...
    \end{align}
    $$
    
- $n_s = a^\D_s a_s$: Occupation of mode $s$

### Procedure

1. compute harmonic flux in q-space ${\bf J}^{\rm ha-{\bf q}}$:
    $$
    {\bf J}^{\rm ha-{\bf q}} = \frac{1}{V} \sum_s \omega_s^2 n_s (t) {\bf v}_s
    $$

2. compute lifetimes $\tau_s = c_s^{-1} \int {\rm d}t \langle n_s(t) n_s \rangle$ with $c_s = \langle n_s^2 \rangle$ for commensurate q-points by fitting $g_s (t) = {\rm e}^{- t / \tau_s}$

3. define $\lambda_s = \omega^2_s \tau_s$ which accounts for $\tau_s \propto \omega_s^{-2}$ scaling of lifetimes [Herring], and interpolate $\lambda_s$ linearly in q on dense mesh in the BZ via linear barycentric interpolation [https://dahtah.wordpress.com/2013/03/06/barycentric-interpolation-fast-interpolation-on-arbitrary-grids/], [August Ferdinand Möbius: *Der barycentrische Calcul*, Verlag von Johann Ambrosius Barth, Leipzig, 1827.]

    1. $$
        {\bf q} = \sum_i \alpha_i {\bf q}_i^{\rm nn}~,
        $$

        with $\alpha_i > 0$ and $\sum_i \alpha_i = 1$. Then
        
    2. $$
    \lambda_b ({\bf q}) \approx \sum_i \alpha_i \lambda_b ({\bf q}^{\rm nn}_i)~,
       $$
    
        where ${\bf q}^{\rm nn}_i$ are the three nearest-neighboring commensurate q-points
    
    3. also get values at $\Gamma$ via this interpolation
    
4. compute thermal conductivity
    $$
    \kappa = c_V \sum_{b, {\bf q}} \omega^4_b ({\bf q}) \lambda_b ({\bf q}) v^2_b ({\bf q})~,
    $$
    where $c_V$ is fixed by requiring that $\kappa = \sum_{b, {\bf q}} c_b ({\bf q}) \, \omega^4_b ({\bf q}) \lambda_b ({\bf q}) v^2_b ({\bf q})$ at the commensurate q-points with $c_s = \langle n^2_s \rangle$. REM: $c_V$ is related to heat capacity, but it's not quite it. Choose different symbol? Mode heat capacity would be $\langle E^2_s \rangle$.

5. compute for q-meshes of increasing density $\rho_{\bf q}$,

6. extrapolate $\kappa_\rho \to_{\rho \to \infty} \kappa^{\rm interpolated}$ _assuming linear convergence with $1/\rho$ [Chui1970,PolyaSzegö1925] [proofs](https://www.whitman.edu/Documents/Academics/Mathematics/2014/owensla.pdf), REM: endpoint rule: evaluate at grid points, midpoint rule: evaluate between grid points -> naive sum is endpoint rule