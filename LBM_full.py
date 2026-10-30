#!/usr/bin/env python3
"""
2D Lattice Boltzmann Method (D2Q9) simulation of constricted channel flow.
Real dimensions (L = 10 cm), air properties, Re_L ≈ 338.
Compares smooth (Venturi-like) vs abrupt constriction.
Includes static-pressure recovery analysis.

Author: generated for educational / GitHub use
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Physical parameters (air, real dimensions)
# =============================================================================
L_phys = 0.10          # m  (10 cm)
H_phys = 0.04          # m  (4 cm)
rho_phys = 1.225       # kg/m³
mu_phys = 1.81e-5      # Pa·s
u_phys = 0.05          # m/s
Re_L = rho_phys * u_phys * L_phys / mu_phys
print(f"Physical Re_L = {Re_L:.1f}")

# =============================================================================
# Lattice parameters (chosen for stability and clean separation)
# =============================================================================
nx, ny = 360, 72
u_lb = 0.06
nu_lb = u_lb * (nx - 1) / Re_L
tau = 3.0 * nu_lb + 0.5
if tau < 0.55:
    tau = 0.58
    nu_lb = (tau - 0.5) / 3.0
    Re_L_eff = u_lb * (nx - 1) / nu_lb
    u_phys = Re_L_eff * mu_phys / (rho_phys * L_phys)
    print(f"Stability adjustment: tau={tau:.3f}, Re_L_eff={Re_L_eff:.0f}, u_phys={u_phys:.4f} m/s")
else:
    print(f"tau = {tau:.3f}")

# D2Q9 lattice
c = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
              [1, 1], [-1, 1], [-1, -1], [1, -1]], dtype=np.int32)
w = np.array([4./9, 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36])
opp = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

# =============================================================================
# Geometry
# =============================================================================
def geometry(kind):
    """Return solid mask (True = wall)."""
    solid = np.zeros((nx, ny), dtype=bool)
    solid[:, 0] = solid[:, -1] = True
    cy = ny // 2
    for i in range(nx):
        x = i / (nx - 1.0)
        if kind == "smooth":
            t = np.exp(-((x - 0.5) / 0.16)**2)
            half = int((ny//2 - 4) * (1 - 0.52 * t))
        else:  # abrupt
            if abs(x - 0.5) < 0.07:
                half = 9
            else:
                half = ny//2 - 3
        solid[i, :max(1, cy - half)] = True
        solid[i, min(ny - 1, cy + half):] = True
    return solid

# =============================================================================
# LBM core
# =============================================================================
def init(solid):
    f = np.zeros((9, nx, ny))
    rho = np.ones((nx, ny))
    for k in range(9):
        f[k] = w[k] * rho
    return f

def step(f, solid, tau, u_in):
    # macroscopic moments
    rho = np.maximum(np.sum(f, axis=0), 1e-8)
    ux = np.sum(c[:, 0, None, None] * f, axis=0) / rho
    uy = np.sum(c[:, 1, None, None] * f, axis=0) / rho
    ux[solid] = uy[solid] = 0.0

    # equilibrium
    feq = np.empty_like(f)
    usq = ux**2 + uy**2
    for k in range(9):
        cu = 3.0 * (c[k, 0] * ux + c[k, 1] * uy)
        feq[k] = w[k] * rho * (1.0 + cu + 0.5 * cu**2 - 1.5 * usq)

    # collision (BGK)
    f = f - (f - feq) / tau

    # streaming
    for k in range(9):
        f[k] = np.roll(np.roll(f[k], c[k, 0], axis=0), c[k, 1], axis=1)

    # bounce-back on solids
    for k in range(1, 9):
        f[k, solid] = f[opp[k], solid]

    # Zou-He inlet (left)
    for j in range(1, ny - 1):
        if not solid[0, j]:
            rho_in = (f[0, 0, j] + f[2, 0, j] + f[4, 0, j]
                      + 2.0 * (f[3, 0, j] + f[6, 0, j] + f[7, 0, j])) / (1.0 - u_in)
            f[1, 0, j] = f[3, 0, j] + (2.0 / 3.0) * rho_in * u_in
            f[5, 0, j] = f[7, 0, j] + (1.0 / 6.0) * rho_in * u_in - 0.5 * (f[2, 0, j] - f[4, 0, j])
            f[8, 0, j] = f[6, 0, j] + (1.0 / 6.0) * rho_in * u_in + 0.5 * (f[2, 0, j] - f[4, 0, j])

    # open outlet (zero gradient)
    f[:, -1, :] = f[:, -2, :]
    return f, ux, uy, rho

def simulate(kind, steps=10000):
    print(f"Simulating {kind}...")
    solid = geometry(kind)
    f = init(solid)
    for n in range(steps):
        f, ux, uy, rho = step(f, solid, tau, u_lb)
        if n % 2500 == 0:
            print(f"  step {n}: max|u|_lb = {np.nanmax(np.sqrt(ux**2 + uy**2)):.4f}")
    speed = np.sqrt(ux**2 + uy**2)
    speed[solid] = np.nan
    scale = u_phys / u_lb
    pressure = rho / 3.0          # lattice static pressure
    pressure[solid] = np.nan
    return solid, ux * scale, uy * scale, speed * scale, pressure

# =============================================================================
# Pressure recovery analysis (same definitions used earlier)
# =============================================================================
def pressure_recovery(data, name):
    """Compute upstream / throat / downstream pressures and recovery metrics."""
    solid, u, v, speed, pressure = data
    mid = ny // 2
    # stations
    up_x = 40
    # throat = location of maximum centerline speed
    center_speed = speed[:, mid]
    th_x = int(np.nanargmax(center_speed))
    th_cm = th_x * (L_phys * 100) / (nx - 1)
    dn_x = nx - 40

    def avg_p(x0, width=4):
        vals = []
        for x in range(max(0, x0 - width), min(nx, x0 + width + 1)):
            fluid = ~np.isnan(pressure[x])
            if np.any(fluid):
                vals.append(np.nanmean(pressure[x, fluid]))
        return np.mean(vals) if vals else np.nan

    P_up = avg_p(up_x)
    P_th = avg_p(th_x)
    P_dn = avg_p(dn_x)

    drop = P_up - P_th
    recovered = P_dn - P_th
    ratio = recovered / drop if drop > 1e-8 else np.nan
    net_loss = P_up - P_dn

    # kinetic-energy recovery (centerline)
    u_th = speed[th_x, mid]
    u_dn = speed[dn_x, mid]
    ke_ideal = 0.5 * (u_th**2 - u_dn**2)   # lattice units, rho=1
    actual = P_dn - P_th
    eff = actual / ke_ideal if ke_ideal > 1e-8 else np.nan

    print(f"\n=== {name} pressure recovery ===")
    print(f"  Upstream P   = {P_up:.5f}")
    print(f"  Throat P     = {P_th:.5f}  (x={th_cm:.1f} cm)")
    print(f"  Downstream P = {P_dn:.5f}")
    print(f"  Throat drop  = {drop:.5f}")
    print(f"  Recovered    = {recovered:.5f}")
    print(f"  Recovery ratio = {ratio:.3f}")
    print(f"  Net loss     = {net_loss:.5f}")
    print(f"  Ideal KE recovery = {ke_ideal:.5f}")
    print(f"  Recovery efficiency = {eff:.3f}")
    return dict(P_up=P_up, P_th=P_th, P_dn=P_dn, ratio=ratio, net_loss=net_loss, eff=eff)

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    smooth = simulate("smooth")
    abrupt = simulate("abrupt")

    # pressure analysis
    pressure_recovery(smooth, "Smooth")
    pressure_recovery(abrupt, "Abrupt")

    # plotting
    x_cm = np.linspace(0, L_phys * 100, nx)
    y_cm = np.linspace(0, H_phys * 100, ny)

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True, sharey=True)
    vmax_vel = 0.12
    for col, (data, label) in enumerate([(smooth, "Smooth"), (abrupt, "Abrupt")]):
        solid, u, v, speed, pressure = data
        ax = axes[0, col]
        im_v = ax.imshow(speed.T, origin="lower", cmap="plasma", aspect="auto",
                         vmin=0, vmax=vmax_vel, extent=[0, L_phys*100, 0, H_phys*100])
        X, Y = np.meshgrid(x_cm, y_cm)
        ax.streamplot(X[::4, ::4], Y[::4, ::4],
                      u.T[::4, ::4], v.T[::4, ::4],
                      color="white", density=0.85, linewidth=0.55, arrowsize=0.55)
        ax.grid(True, color="cyan", alpha=0.25, linestyle="--")
        ax.set_title(f"{label}: Velocity [m/s]")

        ax = axes[1, col]
        im_p = ax.imshow(pressure.T, origin="lower", cmap="coolwarm", aspect="auto",
                         extent=[0, L_phys*100, 0, H_phys*100])
        ax.contour(x_cm, y_cm, pressure.T, levels=9, colors="k", alpha=0.25, linewidths=0.4)
        ax.grid(True, color="gray", alpha=0.25, linestyle="--")
        ax.set_title(f"{label}: Lattice pressure (P=ρ/3)")
        ax.set_xlabel("x [cm]")

    axes[0, 0].set_ylabel("y [cm]")
    axes[1, 0].set_ylabel("y [cm]")
    fig.colorbar(im_v, ax=axes[0, :].tolist(), label="|velocity| [m/s]", shrink=0.85)
    fig.colorbar(im_p, ax=axes[1, :].tolist(), label="lattice pressure", shrink=0.85)
    plt.suptitle(f"2D LBM — L={L_phys*100:.0f} cm, air, Re_L≈{Re_L:.0f}\n"
                 f"Inlet {u_phys:.3f} m/s, outlet → atmosphere", fontsize=12)
    plt.tight_layout()
    plt.savefig("lbm_constricted_flow.png", dpi=140, bbox_inches="tight")
    print("\nSaved lbm_constricted_flow.png")
