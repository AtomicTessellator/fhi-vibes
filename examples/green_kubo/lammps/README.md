lammps Green Kubo Ar
===

[reference heat flux](https://lammps.sandia.gov/doc/compute_heat_flux.html)

## Background

$\begin{aligned} \mathbf{J} &=\frac{1}{V}\left[\sum_{i} e_{i} \mathbf{v}_{i}-\sum_{i} \mathbf{S}_{i} \mathbf{v}_{i}\right] \\ &=\frac{1}{V}\left[\sum_{i} e_{i} \mathbf{v}_{i}+\sum_{i<j}\left(\mathbf{f}_{i j} \cdot \mathbf{v}_{j}\right) \mathbf{x}_{i j}\right] \\ &=\frac{1}{V}\left[\sum_{i} e_{i} \mathbf{v}_{i}+\frac{1}{2} \sum_{i<j}\left(\mathbf{f}_{i j} \cdot\left(\mathbf{v}_{i}+\mathbf{v}_{j}\right)\right) \mathbf{x}_{i j}\right] \end{aligned}$

$\kappa=\frac{V}{k_{B} T^{2}} \int_{0}^{\infty}\left\langle J_{x}(0) J_{x}(t)\right\rangle d t=\frac{V}{3 k_{B} T^{2}} \int_{0}^{\infty}\langle\mathbf{J}(0) \cdot \mathbf{J}(t)\rangle d t$

## Flags

`compute ID group-ID heat/flux ke-ID pe-ID stress-ID`

This compute calculates 

- 6 quantities and stores them in a 6-component vector. The 

- first 3 components are the x, y, z components of the full heat flux vector, i.e. (Jx, Jy, Jz). 

- The next 3 components are the x, y, z components of just the convective portion of the flux, i.e. the first term in the equation for J above.

The heat flux can be output every so many timesteps (e.g. via the *thermo_style* custom command). Then as a post-processing operation, an **auto-correlation** can be performed, its **integral estimated**, and the **Green-Kubo formula above evaluated**.

The _fix ave/correlate command_ can calculate the auto-correlation. The *trap() function* in the variable command can calculate the integral.

## Output info

This compute calculates a 

- global vector of length 6 
  
  - (total heat flux vector, followed by convective heat flux vector), which 

- can be accessed by indices 1-6. These values

- can be used by any command that uses global vector values from a compute as input. See the Howto output doc page for an overview of LAMMPS output options.

The vector values calculated by this compute are 

- “extensive”, meaning they scale with the number of atoms in the simulation. They 

- can be divided by the appropriate volume to get a flux, which would then be an “intensive” value, meaning independent of the number of atoms in the simulation. Note that 

- if the compute is “all”, then the *appropriate volume to divide by is the simulation box volume*. However, 

- if a sub-group is used, it should be the volume containing those atoms.

#### Thermo command

[reference thermo](https://lammps.sandia.gov/doc/thermo.html)

```bash
thermo N
```

Compute and print thermodynamic info every _N_ timesteps

```bash
thermo_style style args
```

### fix ave/correlate command

[reference ave/correlate](https://lammps.sandia.gov/doc/fix_ave_correlate.html)

```bash
fix ID group-ID ave/correlate Nevery Nrepeat Nfreq value1 ...
```

The Nevery, Nrepeat, and Nfreq arguments specify on what timesteps the input values will be used to calculate correlation data.

| argument  | meaning                                                  |
| --------- | -------------------------------------------------------- |
| `Nevery`  | use input values every this many timesteps               |
| `Nrepeat` | # of correlation time windows to accumulate              |
| `Nfreq`   | calculate time window averages every this many timesteps |

#### Example

For example, if **Nevery=10**, **Nrepeat=200**, and **Nfreq=2000**, then 

- values on timesteps 0,10,20,…,2000 will be used to compute the final averages on timestep 2000 
- 200 averages will be computed: Cij(0), Cij(10), Cij(20), ..., and Cij(2000). 
- Cij(30) on timestep 2000 will be the average of 199 samples, namely 
  - Vi(0)*Vj(30), Vi(10)*Vj(40), …, Vi(1980)*Vj(1990), Vi(1990)*Vj(2000)
- Cij(30) on timestep 6000 will be the average of 599 samples, namely 
  - Vi(0)*Vj(30), Vi(10)*Vj(40), …, Vi(5980)*Vj(5990), Vi(5990)*Vj(6000)
- and so on and so on

**ave**

If the *ave* setting is running, then the accumulation is never zeroed. Thus the output of correlation data at any timestep is the 

- average over samples accumulated every *Nevery* steps since the fix was defined.

it can only be restarted by deleting the fix via the unfix command, or by re-defining the fix by re-specifying it.

## Parameter

| Parameter/Variable | Value                                                       |
| ------------------ | ----------------------------------------------------------- |
| `vol`              | _a_ = 5.476, 4x4x4 Box -> 9943.92 Å**3                      |
| `scale`            | `28.7289125255705` [=`convert` * `dt`  / (`V * kB * T**2`)] |
| `convert`          | `4.83166430676946e-16`                                      |
| `kB`               | `1.3806504e-23`                                             |
| `dt`               | `4 * 10` (because step `s = 10`)                            |

**remarks**
volume: _a_ = 5.476, 4x4x4 Box -> 9943.92 Å**3

### Sample input script

```bash
# Sample LAMMPS input script for thermal conductivity of solid Ar

units       real
variable    T equal 70
variable    V equal vol
variable    dt equal 4.0
variable    s equal 10      # sample interval     = Nevery
variable    p equal 200     # correlation length  = Nrepeat
variable    d equal $p*$s   # dump interval       = Nfreq

# convert from LAMMPS real units to SI

variable    kB equal 1.3806504e-23    # [J/K] Boltzmann
variable    kCal2J equal 4186.0/6.02214e23
variable    A2m equal 1.0e-10
variable    fs2s equal 1.0e-15
variable    convert equal ${kCal2J}*${kCal2J}/${fs2s}/${A2m}

# setup problem

dimension    3
boundary     p p p
lattice      fcc 5.376 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region       box block 0 4 0 4 0 4
create_box   1 box
create_atoms 1 box
mass         1 39.948
pair_style   lj/cut 13.0
pair_coeff   * * 0.2381 3.405
timestep     ${dt}
thermo       $d

# equilibration and thermalization

velocity     all create $T 102486 mom yes rot yes dist gaussian
fix          NVT all nvt temp $T $T 10 drag 0.2
run          8000

# thermal conductivity calculation, switch to NVE if desired

#unfix       NVT
#fix         NVE all nve

reset_timestep 0
compute      myKE all ke/atom
compute      myPE all pe/atom
compute      myStress all stress/atom NULL virial
compute      flux all heat/flux myKE myPE myStress
variable     Jx equal c_flux[1]/vol
variable     Jy equal c_flux[2]/vol
variable     Jz equal c_flux[3]/vol
fix          JJ all ave/correlate $s $p $d &
             c_flux[1] c_flux[2] c_flux[3] type auto file J0Jt.dat ave running
variable     scale equal ${convert}/${kB}/$T/$T/$V*$s*${dt}
variable     k11 equal trap(f_JJ[3])*${scale}
variable     k22 equal trap(f_JJ[4])*${scale}
variable     k33 equal trap(f_JJ[5])*${scale}
thermo_style custom step temp v_Jx v_Jy v_Jz v_k11 v_k22 v_k33
run          100000
variable     k equal (v_k11+v_k22+v_k33)/3.0
variable     ndens equal count(all)/vol
print        "average conductivity: $k[W/mK] @ $T K, ${ndens}
```
