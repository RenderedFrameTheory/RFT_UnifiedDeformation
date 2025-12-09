import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1. CONSTANTS AND RFT PARAMETERS
# ============================================================

# Gravitational constant in astrophysical units:
# G = 4.3009e-6 kpc (km/s)^2 / Msun
G_KPC = 4.3009e-6

# Unit conversions
MYR_TO_S  = 3.15576e13     # seconds in 1 Myr
KPC_TO_KM = 3.0857e16      # km in 1 kpc

# RFT deformation:
# a0(Gamma_eff) = A0_STAR * Gamma_eff
Gamma_eff = 1.99e-2                 # example effective deformation strength
A0_ref    = 2500.0                  # (km/s)^2 / kpc at reference Gamma
A0_STAR   = A0_ref / 1.99e-2
a0        = A0_STAR * Gamma_eff     # (km/s)^2 / kpc

print(f"Gamma_eff = {Gamma_eff:.3e}")
print(f"a0 (RFT)   = {a0:.2f} (km/s)^2 / kpc")

# Disk + central mass (baryons only)
M_central   = 1.2e11   # Msun, central bulge + inner disk
M_disk      = 6.0e10   # Msun, extended disk
N_particles = 300      # number of disk particles
m_particle  = M_disk / N_particles

# Softening length to avoid divergences
softening_kpc = 0.1


# ============================================================
# 2. GRAVITY LAWS
# ============================================================

def g_newton_central(r_kpc, M=M_central):
    """Newtonian acceleration magnitude from central mass only, in (km/s)^2 / kpc."""
    return G_KPC * M / (r_kpc**2)

def g_rft_central(r_kpc, M=M_central, a0_val=a0):
    """RFT-modified acceleration magnitude from central mass."""
    gN = g_newton_central(r_kpc, M)
    return 0.5 * (gN + np.sqrt(gN**2 + 4.0 * gN * a0_val))


def compute_accel(positions, masses, model="newton"):
    """
    Compute acceleration on each particle (disk + central mass).
    - Disk–disk interactions are always Newtonian (short-range, baryon dominated).
    - Central mass contribution uses Newton OR RFT depending on 'model'.
    Returns accelerations in km/s^2.
    """
    x = positions[:, 0]
    y = positions[:, 1]
    N = positions.shape[0]

    # --- Disk–disk Newtonian gravity (pairwise) ---
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    r2 = dx*dx + dy*dy + softening_kpc**2

    # Remove self-force
    np.fill_diagonal(r2, np.inf)

    inv_r = 1.0 / np.sqrt(r2)
    inv_r3 = inv_r / r2

    ax_dd = -G_KPC * np.sum(masses[None, :] * dx * inv_r3, axis=1)
    ay_dd = -G_KPC * np.sum(masses[None, :] * dy * inv_r3, axis=1)

    # --- Central mass contribution (Newton or RFT) ---
    r2_c = x*x + y*y + softening_kpc**2
    r_c  = np.sqrt(r2_c)

    gN_c = G_KPC * M_central / r2_c   # (km/s)^2 / kpc

    if model == "newton":
        g_c = gN_c
    elif model == "rft":
        g_c = 0.5 * (gN_c + np.sqrt(gN_c**2 + 4.0 * gN_c * a0))
    else:
        raise ValueError("model must be 'newton' or 'rft'")

    ax_c = -g_c * x / r_c
    ay_c = -g_c * y / r_c

    # Total acceleration in (km/s)^2 / kpc
    ax_total = ax_dd + ax_c
    ay_total = ay_dd + ay_c

    # Convert to km/s^2
    ax_km_s2 = ax_total / KPC_TO_KM
    ay_km_s2 = ay_total / KPC_TO_KM

    return np.stack([ax_km_s2, ay_km_s2], axis=1)


# ============================================================
# 3. INITIAL DISK SETUP
# ============================================================

def init_disk(n_particles=N_particles, r_min=2.0, r_max=20.0):
    """
    Build a simple exponential-like disk:
    - radii distributed between r_min and r_max (kpc)
    - angles uniform in [0, 2π)
    - initial velocities set to approximate circular speeds under Newtonian M_enc
    """
    rng = np.random.default_rng(42)

    # Sample r^2 uniformly for more mass at small radii
    u = rng.random(n_particles)
    r = np.sqrt(r_min**2 + (r_max**2 - r_min**2) * u)
    theta = 2.0 * np.pi * rng.random(n_particles)

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    positions = np.stack([x, y], axis=1)

    # Approximate enclosed mass for initial circular speed
    # Inner: central mass dominates, outer: disk contributes gradually
    M_enc = M_central + M_disk * (r / r_max)
    gN = G_KPC * M_enc / r**2
    v_circ = np.sqrt(r * gN)  # km/s

    # Tangential velocities (perpendicular to radius)
    vx = -v_circ * np.sin(theta)
    vy =  v_circ * np.cos(theta)
    velocities = np.stack([vx, vy], axis=1)

    return positions, velocities


# ============================================================
# 4. ORBIT INTEGRATION (LEAPFROG)
# ============================================================

def leapfrog_disk(model="newton", n_steps=3000, dt_myr=0.2):
    """
    Evolve the disk under the chosen gravity model.
    - model: 'newton' or 'rft'
    Returns final (positions, velocities).
    """
    positions, velocities = init_disk()
    masses = np.full(N_particles, m_particle)

    dt_s = dt_myr * MYR_TO_S

    # Kick-drift-kick (leapfrog)
    acc = compute_accel(positions, masses, model=model)
    velocities += 0.5 * acc * dt_s

    for _ in range(n_steps):
        positions += (velocities * dt_s) / KPC_TO_KM
        acc = compute_accel(positions, masses, model=model)
        velocities += acc * dt_s

    return positions, velocities


# ============================================================
# 5. ROTATION CURVE ESTIMATION
# ============================================================

def rotation_curve(positions, velocities, n_bins=15):
    """
    Estimate rotation curve from particle distribution:
    - Compute radius r and tangential speed v_t for each particle.
    - Bin by radius and compute mean v_t per bin.
    """
    x = positions[:, 0]
    y = positions[:, 1]
    vx = velocities[:, 0]
    vy = velocities[:, 1]

    r = np.sqrt(x*x + y*y)

    # Tangential unit vector: (-y/r, x/r)
    # v_t = v · \hat{theta}
    with np.errstate(divide="ignore", invalid="ignore"):
        inv_r = np.where(r > 0, 1.0 / r, 0.0)
        tx = -y * inv_r
        ty =  x * inv_r
        v_t = vx * tx + vy * ty

    # Bin in radius
    r_min = 2.0
    r_max = 40.0
    bins = np.linspace(r_min, r_max, n_bins + 1)
    idx = np.digitize(r, bins)

    r_centers = []
    v_means   = []

    for b in range(1, n_bins + 1):
        mask = idx == b
        if np.any(mask):
            r_centers.append(r[mask].mean())
            v_means.append(np.mean(np.abs(v_t[mask])))

    return np.array(r_centers), np.array(v_means)


# ============================================================
# 6. RUN EXPERIMENT AND PLOT
# ============================================================

def main():
    # Run Newtonian evolution
    print("Integrating disk under Newtonian gravity...")
    posN, velN = leapfrog_disk(model="newton", n_steps=3000, dt_myr=0.2)

    # Run RFT evolution
    print("Integrating disk under RFT gravity...")
    posR, velR = leapfrog_disk(model="rft", n_steps=3000, dt_myr=0.2)

    # Rotation curves
    rN, vN = rotation_curve(posN, velN, n_bins=15)
    rR, vR = rotation_curve(posR, velR, n_bins=15)

    # --- Plot rotation curves ---
    plt.figure(figsize=(7,5))
    plt.plot(rN, vN, "o--", label="Newton (baryons only)")
    plt.plot(rR, vR, "o-",  label="RFT gravity (Gamma_eff)")

    plt.axvspan(10, 40, color="grey", alpha=0.1, label="10–40 kpc")
    plt.axhspan(220, 240, color="orange", alpha=0.1, label="220–240 km/s")

    plt.xlabel("Radius r (kpc)")
    plt.ylabel("v_rot (km/s)")
    plt.title("Toy Galaxy Rotation Curve: Newton vs RFT")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # --- Plot final spatial distribution for both models ---
    fig, axes = plt.subplots(1, 2, figsize=(10,5), sharex=True, sharey=True)

    axes[0].scatter(posN[:,0], posN[:,1], s=4, alpha=0.7)
    axes[0].set_title("Newton disk (final)")
    axes[0].set_xlabel("x (kpc)")
    axes[0].set_ylabel("y (kpc)")
    axes[0].grid(True)
    axes[0].set_aspect("equal", "box")

    axes[1].scatter(posR[:,0], posR[:,1], s=4, alpha=0.7, color="taborange")
    axes[1].set_title("RFT disk (final)")
    axes[1].set_xlabel("x (kpc)")
    axes[1].grid(True)
    axes[1].set_aspect("equal", "box")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
