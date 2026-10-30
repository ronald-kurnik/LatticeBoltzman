# Constricted Channel Flow with the Lattice Boltzmann Method

2-D simulation of air flowing through a duct that contains either a smooth (Venturi-like) or an abrupt constriction.  
The code uses the D2Q9 lattice Boltzmann method, is dimensionally scaled to a 10 cm duct, and reports static-pressure recovery.

## Physical setup

| Quantity              | Value                  |
|-----------------------|------------------------|
| Length \(L\)          | 10 cm                  |
| Height \(H\)          | 4 cm                   |
| Fluid                 | air                    |
| Density \(\rho\)      | \(1.225\,\mathrm{kg\,m^{-3}}\) |
| Viscosity \(\mu\)     | \(1.81\times10^{-5}\,\mathrm{Pa\,s}\) |
| Inlet velocity \(u_\mathrm{in}\) | \(0.05\,\mathrm{m\,s^{-1}}\) |
| Reynolds number \(\mathrm{Re}_L\) | \(\approx 338\)     |
| Outlet                | atmospheric pressure   |

## Lattice Boltzmann Method (D2Q9) – brief introduction

The lattice Boltzmann method (LBM) solves a discrete-velocity Boltzmann equation instead of the Navier–Stokes equations directly.  
On a square lattice the nine discrete velocities (D2Q9) are

\[
\mathbf{c}_i =
\begin{cases}
(0,0) & i=0 \\
(\pm 1,0),\ (0,\pm 1) & i=1\dots 4 \\
(\pm 1,\pm 1) & i=5\dots 8
\end{cases}
\]

with weights

\[
w_0=\frac{4}{9},\quad
w_{1-4}=\frac{1}{9},\quad
w_{5-8}=\frac{1}{36}.
\]

The central quantity is the particle distribution function \(f_i(\mathbf{x},t)\).  
Macroscopic density and velocity are recovered by moments:

\[
\rho = \sum_{i=0}^{8} f_i,
\qquad
\rho\mathbf{u} = \sum_{i=0}^{8} \mathbf{c}_i f_i.
\]

### Collision (BGK)

\[
f_i^* = f_i - \frac{1}{\tau}\bigl(f_i - f_i^\mathrm{eq}\bigr)
\]

where the equilibrium distribution is

\[
f_i^\mathrm{eq}
= w_i\rho\Bigl[
  1 + 3\mathbf{c}_i\cdot\mathbf{u}
  + \frac{9}{2}(\mathbf{c}_i\cdot\mathbf{u})^2
  - \frac{3}{2}|\mathbf{u}|^2
\Bigr].
\]

The relaxation time \(\tau\) is related to the kinematic viscosity by

\[
\nu = c_s^2\bigl(\tau - \tfrac12\bigr)
= \frac13\bigl(\tau - \tfrac12\bigr)
\]

(\(c_s^2=1/3\) in lattice units).

### Streaming

\[
f_i(\mathbf{x}+\mathbf{c}_i,\,t+1) = f_i^*(\mathbf{x},t).
\]

Solid walls are treated by the simple bounce-back rule  
\(f_{\bar{i}} = f_i\) (where \(\bar{i}\) is the opposite direction).

### Pressure

In the LBM equation of state the static pressure is simply

\[
P = \rho\,c_s^2 = \frac{\rho}{3}.
\]

## Pressure-recovery definitions used in the code

Let  
- \(P_\mathrm{up}\) = average static pressure well upstream,  
- \(P_\mathrm{th}\) = static pressure at the throat (location of maximum speed),  
- \(P_\mathrm{dn}\) = average static pressure well downstream.

Then

\[
\begin{align*}
\text{throat drop} &= P_\mathrm{up}-P_\mathrm{th},\\
\text{recovered} &= P_\mathrm{dn}-P_\mathrm{th},\\
\text{recovery ratio} &= \frac{P_\mathrm{dn}-P_\mathrm{th}}{P_\mathrm{up}-P_\mathrm{th}},\\
\text{net loss} &= P_\mathrm{up}-P_\mathrm{dn}.
\end{align*}
\]

The ideal kinetic-energy recovery (Bernoulli) is

\[
\Delta P_\mathrm{ideal}
= \tfrac12\rho\bigl(u_\mathrm{th}^2 - u_\mathrm{dn}^2\bigr)
\]

and the recovery efficiency is

\[
\eta = \frac{P_\mathrm{dn}-P_\mathrm{th}}{\Delta P_\mathrm{ideal}}.
\]

(Because the exit is open to atmosphere, \(\eta\) is often negative; the useful figure of merit is the net loss.)

## Running the code

```bash
python3 lbm_constricted_flow.py
```

Requires only NumPy and Matplotlib.  
It produces `lbm_constricted_flow.png` and prints the pressure-recovery tables for both geometries.

## Expected behaviour

- **Smooth constriction** – flow stays attached; good static-pressure recovery.  
- **Abrupt constriction** – clear laminar recirculation zones form immediately downstream of the sudden expansion; larger net pressure loss.

## References

- Krüger et al., *The Lattice Boltzmann Method*, Springer, 2017.  
- Succi, *The Lattice Boltzmann Equation*, Oxford, 2001.
