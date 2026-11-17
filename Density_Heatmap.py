import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
 
rng = np.random.default_rng(0)
 
# ------------------------------------------------------------------ parameters
V0 = -40e3            # cathode potential [V]  (was labelled eV; the field needs volts)
R_prime = 0.02        # cathode radius [m]
R = 0.25              # anode / outer shell radius [m]
 
n_points = 201         # ODD, so the grid contains the y = 0 plane
q_i = 1.602e-19       # [C]
m_i = 3.344e-27       # [kg] deuteron
K_B = 1.381e-23
T = 293.0
P = 0.613             # 4.6 mTorr [Pa]
n0 = P / (K_B * T)    # neutral background density (kept separate from ion density)
 
dt = 4e-9             # s  (ion speed ~2e6 m/s -> ~0.004 m / step, < dx = 0.01 m)
steps = 1400
burn_in = 1000        # steps before density averaging starts
num_m = 100000         # macroparticles
N = 1000              # real ions per macroparticle
 
# ------------------------------------------------------------------ grid & field
axis = np.linspace(-R, R, n_points)
dx = axis[1] - axis[0]                     # linspace spacing is 2R/(n-1), not 2R/n
x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")   # indexing='ij': array[i,j,k] <-> (x_i, y_j, z_k)
 
r = np.sqrt(x**2 + y**2 + z**2)
r_safe = np.where(r == 0, 1.0, r)
 
# phi(r) = V0 (1/r - 1/R) / (1/R' - 1/R)  ->  E_r = -dphi/dr = V0 / (1/R' - 1/R) / r^2
E_r = V0 / (1 / R_prime - 1 / R) / r_safe**2
E_r[r < R_prime] = 0.0                     # field-free inside the cathode
 
# Cartesian components directly (avoids the 0/0 NaNs from cos_phi, sin_phi on the z-axis)
Ex = E_r * x / r_safe
Ey = E_r * y / r_safe
Ez = E_r * z / r_safe
 
 
# ------------------------------------------------------------------ CIC helpers
def cell_weights(px, py, pz):
    """Return flat indices and trilinear weights of the 8 surrounding nodes."""
    fx, fy, fz = (px + R) / dx, (py + R) / dx, (pz + R) / dx
    i0 = np.clip(np.floor(fx).astype(int), 0, n_points - 2)
    j0 = np.clip(np.floor(fy).astype(int), 0, n_points - 2)
    k0 = np.clip(np.floor(fz).astype(int), 0, n_points - 2)
    wx, wy, wz = fx - i0, fy - j0, fz - k0
    idx, w = [], []
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                idx.append(((i0 + di) * n_points + (j0 + dj)) * n_points + (k0 + dk))
                w.append((wx if di else 1 - wx) *
                         (wy if dj else 1 - wy) *
                         (wz if dk else 1 - wz))
    return idx, w
 
 
def gather(field_flat, idx, w):
    return sum(wi * field_flat[ii] for ii, wi in zip(idx, w))
 
 
Ex_f, Ey_f, Ez_f = Ex.ravel(), Ey.ravel(), Ez.ravel()
 
 
# ------------------------------------------------------------------ particles
def spawn(n):
    """Uniform in the shell R' < r < R, thermal velocities."""
    rr = (R_prime**3 + rng.random(n) * (R**3 - R_prime**3)) ** (1 / 3)
    cost = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    sint = np.sqrt(1 - cost**2)
    pos = np.stack([rr * sint * np.cos(phi), rr * sint * np.sin(phi), rr * cost])
    vth = np.sqrt(K_B * T / m_i)
    vel = rng.normal(0, vth, (3, n))
    return pos, vel
 
 
(px, py, pz), (vx, vy, vz) = spawn(num_m)
 
# Leapfrog start: pull velocities back half a step so v lives at t - dt/2
idx, w = cell_weights(px, py, pz)
qm = q_i / m_i
vx -= 0.5 * qm * gather(Ex_f, idx, w) * dt
vy -= 0.5 * qm * gather(Ey_f, idx, w) * dt
vz -= 0.5 * qm * gather(Ez_f, idx, w) * dt
 
density_sum = np.zeros(n_points**3)
n_avg = 0
 
# ------------------------------------------------------------------ main loop
for step in range(steps):          # (np.arange(0, dt*steps) only gave ONE iteration)
    idx, w = cell_weights(px, py, pz)
 
    # Boris pusher with B = 0 reduces to leapfrog: v += (q/m) E dt ; x += v dt
    vx += qm * gather(Ex_f, idx, w) * dt
    vy += qm * gather(Ey_f, idx, w) * dt
    vz += qm * gather(Ez_f, idx, w) * dt
    px += vx * dt
    py += vy * dt
    pz += vz * dt
 
    # Ions that reach the outer shell are lost; re-inject so the population stays constant
    lost = (px**2 + py**2 + pz**2) >= R**2
    if lost.any():
        (npx, npy, npz), (nvx, nvy, nvz) = spawn(lost.sum())
        px[lost], py[lost], pz[lost] = npx, npy, npz
        vx[lost], vy[lost], vz[lost] = nvx, nvy, nvz
 
    # Ion density: cloud-in-cell deposit, time-averaged after burn-in
    if step >= burn_in:
        idx, w = cell_weights(px, py, pz)
        dep = np.zeros(n_points**3)
        for ii, wi in zip(idx, w):
            dep += np.bincount(ii, weights=wi, minlength=n_points**3)
        density_sum += dep
        n_avg += 1
 
n_i = (density_sum / n_avg) * N / dx**3          # [m^-3]
n_i = n_i.reshape(n_points, n_points, n_points)
 
# ------------------------------------------------------------------ x-z cross-section (y = 0)
j0 = n_points // 2

slice_xz = 0.0
for i in range(-2,3):
    slice_xz += n_i[:,j0+i,:]            # axes: (x, z)
 
fig, ax = plt.subplots(figsize=(7, 6))
pos = slice_xz[slice_xz > 0]
mesh = ax.pcolormesh(axis, axis, slice_xz.T,       # transpose: x horizontal, z vertical
                     shading="nearest", cmap="inferno",
                     norm=LogNorm(vmin=pos.min(), vmax=slice_xz.max()))
ax.add_patch(plt.Circle((0, 0), R_prime, fill=False, color="cyan", ls="--", label="cathode"))
ax.add_patch(plt.Circle((0, 0), R, fill=False, color="white", ls=":", label="outer shell"))
ax.set_aspect("equal")
ax.set_xlabel("x [m]")
ax.set_ylabel("z [m]")
ax.set_title("Ion number density, -3dy < y < 3dy cross-section")
ax.legend(loc="upper right")
fig.colorbar(mesh, label=r"$n_i$ [m$^{-3}$]")
fig.tight_layout()
fig.savefig("C:/Users/salaz/GTFE-Fusor/density_xz.png", dpi=150)
print("max density: %.3e m^-3, neutral background: %.3e m^-3" % (n_i.max(), n0))
 

