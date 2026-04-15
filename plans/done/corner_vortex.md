---
jupyter:
  jupytext:
    formats: md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Corner-correcting potential for C-grid interpolation

## Stream function on a C-grid cell

On a NEMO-style C-grid, the velocity inside a cell is linear in one
direction only:

$$u(x, y) = (1 - x) \, U_0 + x \, U_1$$
$$v(x, y) = (1 - y) \, V_0 + y \, V_1$$

with the continuity constraint $U_1 - U_0 + V_1 - V_0 = 0$.

This velocity field is divergence-free by construction. The
corresponding stream function is:

$$\psi(x, y) = U_0 (1-x) y + U_1 x y - V_0 x + C$$

where $V_1$ enters implicitly through continuity. Streamlines are
contours $\psi = \text{const}$.

## Concave corner cell

When a cell has a concave SE corner (land to the east and south),
the land-adjacent edges are masked: $U_1 = 0$, $V_0 = 0$.
Continuity then requires $V_1 = U_0$. The stream function simplifies
to:

$$\psi = U_0 \, (1-x) \, y$$

Streamlines are hyperbolas $(1-x) \, y = K$. All streamlines converge
toward the corner at $(1, 0)$ where $\psi = 0$. The corner is a
stagnation point that traps particles under finite-timestep integration.

A cell can have at most one concave corner, because two concave
corners would require three land neighbors, leaving at most one
non-zero edge — which by continuity forces all edges to zero (no
flow).

```python
import numpy as np
import matplotlib.pyplot as plt
```

```python
x = np.linspace(0, 1, 200)
y = np.linspace(0, 1, 200)
X, Y = np.meshgrid(x, y)

psi_free = (1 - X) * Y

fig, ax = plt.subplots()
ax.contour(X, Y, psi_free, levels=np.linspace(0.005, 0.9, 30), colors="black", linewidths=0.5)
ax.set_aspect("equal")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_title(r"$\psi = (1-x)\,y$ — concave SE corner at (1,0)")
```

## Corner-correcting source potential

We add a correction that pushes particles radially away from the
concave corner. The correction stream function for a source at
$(1, 0)$ is:

$$\psi_{\mathrm{corr}} = \varepsilon \, \mathrm{atan2}(y,\; x - 1)$$

This is the conjugate harmonic of $\ln r$ and gives radial outward
velocities:

$$u_{\mathrm{corr}} = \frac{\partial \psi_{\mathrm{corr}}}{\partial y} = \varepsilon \, \frac{x - 1}{r^2}$$
$$v_{\mathrm{corr}} = -\frac{\partial \psi_{\mathrm{corr}}}{\partial x} = \varepsilon \, \frac{y}{r^2}$$

where $r^2 = (x-1)^2 + y^2$. Inside the cell ($x < 1$, $y > 0$):

- $u_{\mathrm{corr}} < 0$ (westward, away from the east land edge),
- $v_{\mathrm{corr}} > 0$ (northward, away from the south land edge).

The correction scales as $1/r$ while the ambient velocity scales as
$\sim r$, so the correction dominates close to the corner for any
$\varepsilon > 0$.

The total stream function is:

$$\psi_{\mathrm{total}} = U_0 \, (1-x) \, y + \varepsilon \, \mathrm{atan2}(y,\; x - 1)$$

```python
def psi_total(X, Y, epsilon):
    return (1 - X) * Y + epsilon * np.arctan2(Y, X - 1)


fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))

for ax, eps in zip(axes, [0.0, 0.002, 0.01, 0.04]):
    psi = psi_total(X, Y, eps)
    ax.contour(X, Y, psi, levels=30, colors="black", linewidths=0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(rf"$\varepsilon = {eps}$")

fig.suptitle(r"$\psi = (1-x)\,y + \varepsilon\,\mathrm{atan2}(y, x-1)$ — source correction")
fig.tight_layout()
```

## Regularization

The source has a $1/r$ singularity at the corner. To regularize,
replace $r^2$ with $r^2 + \delta^2$ in the velocity (equivalently,
smear the point source over a disk of radius $\delta$):

$$u_{\mathrm{corr}} = \varepsilon \, \frac{x - 1}{r^2 + \delta^2}, \quad v_{\mathrm{corr}} = \varepsilon \, \frac{y}{r^2 + \delta^2}$$

At $r = 0$ the velocity is now finite: $|\mathbf{v}| = \varepsilon /
\delta^2$. For $r \gg \delta$ the behavior is unchanged.

```python
def psi_regularized(X, Y, epsilon, delta):
    r2 = (X - 1)**2 + Y**2
    # Regularized atan2: atan2(y, x-1) ≈ arctan(y / (x-1)) with smoothed r
    # For the stream function plot, use the exact atan2 (regularization only
    # affects velocities, not the stream function topology for visualization)
    return (1 - X) * Y + epsilon * np.arctan2(Y, X - 1)


fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))

eps = 0.01
for ax, delta in zip(axes, [0.0, 0.02, 0.05, 0.1]):
    psi = psi_regularized(X, Y, eps, delta)
    ax.contour(X, Y, psi, levels=30, colors="black", linewidths=0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(rf"$\delta = {delta}$")

fig.suptitle(rf"Source correction ($\varepsilon = {eps}$) — stream function is independent of $\delta$")
fig.tight_layout()
```

## Choosing $\varepsilon$

The correction dominates the ambient flow when $\varepsilon / r$
exceeds the ambient velocity $\sim U_0 \, r$. This gives a crossover
radius:

$$r_{\mathrm{cross}} = \sqrt{\varepsilon / U_0}$$

Inside this radius, the source controls the flow. For
$\varepsilon = 0.01$ and $U_0 = 1$, the crossover is at $r = 0.1$
(10% of the cell). At the opposite corner ($r \approx \sqrt{2}$),
the correction velocity is $\varepsilon / \sqrt{2} \approx 0.007$
— less than 1% of the ambient.

```python
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

epsilons = [0.002, 0.01, 0.04]
for ax, eps in zip(axes, epsilons):
    r_cross = np.sqrt(eps)
    psi = psi_total(X, Y, eps)
    ax.contour(X, Y, psi, levels=30, colors="black", linewidths=0.5)
    theta = np.linspace(0, 2 * np.pi, 100)
    ax.plot(1 + r_cross * np.cos(theta), r_cross * np.sin(theta), color="tab:red", linewidth=0.8)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(rf"$\varepsilon = {eps}$, $r_{{\mathrm{{cross}}}} = {r_cross:.3f}$")

fig.suptitle(r"Crossover radius (red) where source dominates ambient flow ($U_0 = 1$)")
fig.tight_layout()
```

## Particle trajectories: uncorrected vs corrected

Release particles along the $(0,1)$–$(1,0)$ diagonal and integrate
with RK4. Without correction, particles near the corner get trapped.
With the regularized source, they are pushed away from the corner and
exit the cell. The third panel seeds particles directly into the
corner zone.

```python
def velocity(x, y, epsilon=0.0, delta=0.0):
    """C-grid velocity for SE corner cell (U0=1, U1=0, V0=0, V1=1) + source."""
    u = 1.0 - x
    v = y
    if epsilon > 0:
        r2 = (x - 1.0)**2 + y**2 + delta**2
        u += epsilon * (x - 1.0) / r2
        v += epsilon * y / r2
    return u, v


def rk4_step(x, y, dt, **kw):
    u1, v1 = velocity(x, y, **kw)
    u2, v2 = velocity(x + 0.5*dt*u1, y + 0.5*dt*v1, **kw)
    u3, v3 = velocity(x + 0.5*dt*u2, y + 0.5*dt*v2, **kw)
    u4, v4 = velocity(x + dt*u3, y + dt*v3, **kw)
    return x + dt*(u1 + 2*u2 + 2*u3 + u4)/6, y + dt*(v1 + 2*v2 + 2*v3 + v4)/6


def integrate(x0, y0, dt, n_steps, tick_every=None, **kw):
    """Integrate, stop at cell boundary. Return trajectory + tick positions."""
    traj = [(x0, y0)]
    ticks = []
    xp, yp = x0, y0
    for i in range(n_steps):
        xp, yp = rk4_step(xp, yp, dt, **kw)
        if xp < 0 or xp > 1 or yp < 0 or yp > 1:
            break
        traj.append((xp, yp))
        if tick_every and (i + 1) % tick_every == 0:
            ticks.append((xp, yp))
    return np.array(traj), np.array(ticks) if ticks else np.empty((0, 2))


# 30 particles uniform on the (0,1)-(1,0) diagonal
t_starts = np.linspace(0.01, 0.99, 30)
starts = [(t, 1.0 - t) for t in t_starts]

# Particles cross the cell in ~2 time units (speed ~ 0.7 at center).
# dt=0.002, tick every 0.2 time units = 100 steps -> ~10 ticks per transit.
# Many total steps so slow near-corner particles are also resolved.
dt = 0.002
n_steps = 50000
tick_every = 100
eps, delta = 0.01, 0.05

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, (kw, label) in zip(axes, [
    (dict(epsilon=0.0, delta=0.0), "Uncorrected"),
    (dict(epsilon=eps, delta=delta),
     rf"Corrected ($\varepsilon={eps}$, $\delta={delta}$)"),
]):
    psi = psi_total(X, Y, kw.get("epsilon", 0.0))
    ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)
    for x0, y0 in starts:
        traj, ticks = integrate(x0, y0, dt, n_steps, tick_every=tick_every, **kw)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.5)
        ax.plot(x0, y0, ".", color="black", markersize=2)
        if len(ticks) > 0:
            ax.plot(ticks[:, 0], ticks[:, 1], ".", color="black", markersize=2)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.tight_layout()
```

Now take the exit positions on the north edge from the corrected run
and use them as starting positions on the west edge (transpose $x
\leftrightarrow y$). This simulates a particle entering the next cell
at the position where it exited the previous one. The low-$\psi$
particles that were rescued by the correction now start close to the
south edge — exactly where they would get trapped without correction.

```python
def find_exit(x0, y0, dt, n_steps, **kw):
    """Integrate until particle exits [0,1]x[0,1]. Return (exit_edge, position_along_edge)."""
    xp, yp = x0, y0
    for _ in range(n_steps):
        xn, yn = rk4_step(xp, yp, dt, **kw)
        if yn > 1:
            return "north", xp + (xn - xp) * (1 - yp) / (yn - yp)  # interpolate x at y=1
        if xn > 1:
            return "east", yp + (yn - yp) * (1 - xp) / (xn - xp)
        if yn < 0:
            return "south", xp + (xn - xp) * (0 - yp) / (yn - yp)
        if xn < 0:
            return "west", yp + (yn - yp) * (0 - xp) / (xn - xp)
        xp, yp = xn, yn
    return "stuck", None


# Collect exit positions from the UNCORRECTED run — this is what
# actually enters the next cell in the unmodified simulation.
exits = []
for x0, y0 in starts:
    edge, pos = find_exit(x0, y0, dt, n_steps, epsilon=0.0, delta=0.0)
    exits.append((edge, pos))
    
# Map each exit to an entry on the opposite edge of the next cell:
# north exit at x -> west entry at y=x (entering from the left of the cell above)
# east exit at y -> south entry at x=y (entering from below the cell to the right)
starts_next = []
for edge, pos in exits:
    if edge == "north":
        starts_next.append((0.0, 1.0 - pos))  # high x at north -> low y at west (near corner)
    elif edge == "east":
        starts_next.append((1.0 - pos, 0.0))  # high y at east -> low x at south (near corner)
    # stuck particles don't propagate

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, (kw, label) in zip(axes, [
    (dict(epsilon=0.0, delta=0.0), "Uncorrected"),
    (dict(epsilon=eps, delta=delta),
     rf"Corrected ($\varepsilon={eps}$, $\delta={delta}$)"),
]):
    psi = psi_total(X, Y, kw.get("epsilon", 0.0))
    ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)
    for x0, y0 in starts_next:
        traj, ticks = integrate(x0, y0, dt, n_steps, tick_every=tick_every, **kw)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.5)
        ax.plot(x0, y0, ".", color="black", markersize=2)
        if len(ticks) > 0:
            ax.plot(ticks[:, 0], ticks[:, 1], ".", color="black", markersize=2)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.suptitle("Second cell: entering at west edge from previous cell's north-edge exits")
fig.tight_layout()
```

The source correction is **flow-direction-independent**: it always
pushes radially away from the corner. If the ambient flow reverses
($U_0 < 0$, $V_1 < 0$ — entering from north, exiting west), the
corner remains a stagnation point, and the correction still pushes
away from it. The detection condition ($U_1 = 0$ and $V_0 = 0$)
holds regardless of the sign of the non-zero edges.

## Exploring higher-order corrections

The source ($\ln \zeta$) gives $1/r$ velocity but creates asymmetry
at the edges. Higher-order complex potentials concentrate the
correction more tightly at the corner. For $\zeta = (x-1) + iy$:

| Potential $w(\zeta)$ | $\psi = \mathrm{Im}(w)$ | Velocity scaling | Name |
|---|---|---|---|
| $\ln \zeta$ | $\mathrm{atan2}(y, x-1)$ | $1/r$ | Source |
| $-1/\zeta$ | $y / r^2$ | $1/r^2$ | Dipole |
| $-1/(2\zeta^2)$ | $y(x-1)/r^4$ | $1/r^3$ | Quadrupole |

All are harmonic, so all derived velocity fields are divergence-free.
The key question: does the velocity direction work in all parts of
the cell?

```python
def psi_source(X, Y, eps):
    return eps * np.arctan2(Y, X - 1)

def psi_dipole(X, Y, eps):
    r2 = (X - 1)**2 + Y**2
    r2 = np.where(r2 > 1e-20, r2, 1e-20)
    return eps * Y / r2

def psi_quadrupole(X, Y, eps):
    r2 = (X - 1)**2 + Y**2
    r2 = np.where(r2 > 1e-20, r2, 1e-20)
    return eps * Y * (X - 1) / r2**2

def vel_source(x, y, eps, delta=0.0):
    r2 = (x - 1)**2 + y**2 + delta**2
    return eps * (x - 1) / r2, eps * y / r2

def vel_dipole(x, y, eps, delta=0.0):
    dx = x - 1
    r2 = dx**2 + y**2 + delta**2
    u = eps * (r2 - 2*y**2) / r2**2
    v = eps * 2*y*dx / r2**2
    return u, v

def vel_quadrupole(x, y, eps, delta=0.0):
    dx = x - 1
    r2 = dx**2 + y**2 + delta**2
    u = eps * dx * (r2 - 4*y**2) / r2**3
    v = eps * y * (r2 - 4*dx**2) / r2**3
    return u, v


print("Velocity direction at key locations (eps=1, delta=0):")
print(f"  {'Location':<25} {'Source (u,v)':>18} {'Dipole (u,v)':>18} {'Quadrupole (u,v)':>18}")
for label, xp, yp in [
    ("near E (0.99, 0.5)", 0.99, 0.5),
    ("near S (0.5, 0.01)", 0.5, 0.01),
    ("corner (0.95, 0.05)", 0.95, 0.05),
    ("diagonal (0.9, 0.1)", 0.9, 0.1),
    ("far (0.1, 0.9)", 0.1, 0.9),
]:
    us, vs = vel_source(xp, yp, 1.0)
    ud, vd = vel_dipole(xp, yp, 1.0)
    uq, vq = vel_quadrupole(xp, yp, 1.0)
    print(f"  {label:<25} ({us:+.3f},{vs:+.3f})  ({ud:+.3f},{vd:+.3f})  ({uq:+.3f},{vq:+.3f})")
```

```python
# Streamline comparison: ambient + correction for each type
# Scale eps so the visible effect is comparable
fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))

psi_base = (1 - X) * Y
corrections = [
    ("No correction", lambda X, Y, e: np.zeros_like(X), 0),
    ("Source ($\\ln \\zeta$)", psi_source, 0.01),
    ("Dipole ($1/\\zeta$)", psi_dipole, 0.002),
    ("Quadrupole ($1/\\zeta^2$)", psi_quadrupole, 0.0005),
]

for ax, (label, psi_fn, eps) in zip(axes, corrections):
    psi = psi_base + psi_fn(X, Y, eps)
    ax.contour(X, Y, psi, levels=30, colors="black", linewidths=0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(f"{label}\n$\\varepsilon={eps}$")

fig.tight_layout()
```

```python
# Particle trajectories with each correction type
def make_velocity_fn(vel_corr, eps, delta):
    def vfn(x, y):
        u, v = 1.0 - x, y
        if vel_corr is not None:
            uc, vc = vel_corr(x, y, eps, delta)
            u += uc
            v += vc
        return u, v
    return vfn

def integrate_vfn(vfn, x0, y0, dt, n_steps):
    traj = [(x0, y0)]
    xp, yp = x0, y0
    for _ in range(n_steps):
        u1, v1 = vfn(xp, yp)
        u2, v2 = vfn(xp + 0.5*dt*u1, yp + 0.5*dt*v1)
        u3, v3 = vfn(xp + 0.5*dt*u2, yp + 0.5*dt*v2)
        u4, v4 = vfn(xp + dt*u3, yp + dt*v3)
        xn = xp + dt*(u1 + 2*u2 + 2*u3 + u4)/6
        yn = yp + dt*(v1 + 2*v2 + 2*v3 + v4)/6
        if xn < 0 or xn > 1 or yn < 0 or yn > 1:
            break
        traj.append((xn, yn))
        xp, yp = xn, yn
    return np.array(traj)

configs = [
    ("Uncorrected", None, 0, 0),
    ("Source", vel_source, 0.01, 0.05),
    ("Dipole", vel_dipole, 0.002, 0.05),
    ("Quadrupole", vel_quadrupole, 0.0005, 0.05),
]

t_starts = np.linspace(0.01, 0.99, 30)

fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))

for ax, (label, vel_fn, eps_c, delta_c) in zip(axes, configs):
    # Background streamlines
    if vel_fn is not None:
        corr_map = {vel_source: psi_source, vel_dipole: psi_dipole, vel_quadrupole: psi_quadrupole}
        psi = psi_base + corr_map[vel_fn](X, Y, eps_c)
    else:
        psi = psi_base
    ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)

    vfn = make_velocity_fn(vel_fn, eps_c, delta_c)
    for t in t_starts:
        traj = integrate_vfn(vfn, t, 1.0 - t, 0.002, 50000)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.4)
        ax.plot(t, 1.0 - t, ".", color="black", markersize=2)

    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.tight_layout()
```

```python
# Second cell: chain exits from uncorrected first cell into each correction type
# Collect uncorrected exits
vfn_uncorr = make_velocity_fn(None, 0, 0)
exits_uncorr = []
for t in t_starts:
    x0, y0 = t, 1.0 - t
    traj = integrate_vfn(vfn_uncorr, x0, y0, 0.002, 50000)
    xf, yf = traj[-1]
    if yf > 0.99:
        exits_uncorr.append(1.0 - xf)  # flip: high x at north -> low y at west

starts_2nd = [(0.0, yy) for yy in exits_uncorr]

fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))

for ax, (label, vel_fn, eps_c, delta_c) in zip(axes, configs):
    if vel_fn is not None:
        corr_map = {vel_source: psi_source, vel_dipole: psi_dipole, vel_quadrupole: psi_quadrupole}
        psi = psi_base + corr_map[vel_fn](X, Y, eps_c)
    else:
        psi = psi_base
    ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)

    vfn = make_velocity_fn(vel_fn, eps_c, delta_c)
    for x0, y0 in starts_2nd:
        traj = integrate_vfn(vfn, x0, y0, 0.002, 50000)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.4)
        ax.plot(x0, y0, ".", color="black", markersize=2)

    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.suptitle("Second cell: entering from uncorrected first cell's north-edge exits")
fig.tight_layout()
```

## Alternative: coordinate remapping

Instead of adding a correction potential, remap the interpolation
point. A particle at $(x, y)$ evaluates the standard C-grid velocity
at $(x', y')$ where the mapping smoothly pulls the corner inward:

$$x' = x - \delta \, g(x, y), \quad y' = y + \delta \, g(x, y)$$

with $g = \exp\!\left(-(1-x)^2/\sigma^2 - y^2/\sigma^2\right)$
concentrated at the SE corner. The standard interpolation at
$(x', y')$ gives nonzero velocity even at the original corner,
because the remapped point is safely inside the cell.

No correction potential, no divergence concerns from added fields.

```python
def vel_remap(x, y, delta_r=0.1, sigma=0.15):
    """Velocity via coordinate remapping near SE corner."""
    g = np.exp(-((1 - x)**2 + y**2) / sigma**2)
    xp = x - delta_r * g
    yp = y + delta_r * g
    # Clamp to cell
    xp = np.clip(xp, 0, 1)
    yp = np.clip(yp, 0, 1)
    return 1.0 - xp, yp


# Show the remapping: original grid vs remapped grid
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: grid distortion
ax = axes[0]
delta_r, sigma = 0.1, 0.15
xg = np.linspace(0, 1, 21)
yg = np.linspace(0, 1, 21)
Xg, Yg = np.meshgrid(xg, yg)
g = np.exp(-((1 - Xg)**2 + Yg**2) / sigma**2)
Xp = Xg - delta_r * g
Yp = Yg + delta_r * g
for j in range(len(yg)):
    ax.plot(Xg[j, :], Yg[j, :], color="0.8", linewidth=0.5)
    ax.plot(Xp[j, :], Yp[j, :], color="black", linewidth=0.5)
for i in range(len(xg)):
    ax.plot(Xg[:, i], Yg[:, i], color="0.8", linewidth=0.5)
    ax.plot(Xp[:, i], Yp[:, i], color="black", linewidth=0.5)
ax.set_aspect("equal")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_title(f"Grid remapping ($\\delta={delta_r}$, $\\sigma={sigma}$)\ngray=original, black=remapped")

# Right: trajectories
ax = axes[1]
psi = (1 - X) * Y
ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)

def vfn_remap(x, y):
    return vel_remap(x, y, delta_r=delta_r, sigma=sigma)

for t in t_starts:
    traj = integrate_vfn(vfn_remap, t, 1.0 - t, 0.002, 50000)
    ax.plot(traj[:, 0], traj[:, 1], linewidth=0.5)
    ax.plot(t, 1.0 - t, ".", color="black", markersize=2)

ax.set_aspect("equal")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_title("Trajectories with remapping")
fig.tight_layout()
```

```python
# Second cell comparison: uncorrected vs source vs remap
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

remap_configs = [
    ("Uncorrected", make_velocity_fn(None, 0, 0)),
    ("Source", make_velocity_fn(vel_source, 0.01, 0.05)),
    ("Remap", vfn_remap),
]

for ax, (label, vfn) in zip(axes, remap_configs):
    psi = (1 - X) * Y
    ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)
    for x0, y0 in starts_2nd:
        traj = integrate_vfn(vfn, x0, y0, 0.002, 50000)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.5)
        ax.plot(x0, y0, ".", color="black", markersize=2)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.suptitle("Second cell: uncorrected exits → uncorrected / source / remap")
fig.tight_layout()
```

### Edge-preserving remapping

The Gaussian remapping doesn't exactly preserve $\psi$ on the free
edges. A better choice: $g(x, y) = x \, (1-y)$, which vanishes on
the west edge ($x=0$) and north edge ($y=1$) by construction:

$$x' = x - \delta \, x \, (1-y) = x\,(1 - \delta(1-y))$$
$$y' = y + \delta \, x \, (1-y)$$

Note that $g = x(1-y)$ is the "complementary" stream function —
largest at the corner $(1,0)$ where $\psi = (1-x)y$ is smallest.

```python
def remap_edge_preserving(x, y, delta_r=0.1):
    g = x * (1 - y)
    xp = x - delta_r * g
    yp = y + delta_r * g
    return xp, yp


def psi_remapped_ep(X, Y, delta_r=0.1):
    Xp, Yp = remap_edge_preserving(X, Y, delta_r)
    return (1 - Xp) * Yp


def vel_remap_ep(x, y, delta_r=0.1):
    xp, yp = remap_edge_preserving(x, y, delta_r)
    return 1.0 - xp, yp


fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))

for ax, dr in zip(axes, [0.0, 0.05, 0.1, 0.2]):
    psi = psi_remapped_ep(X, Y, dr)
    ax.contourf(X, Y, psi, levels=20, cmap="viridis")
    ax.contour(X, Y, psi, levels=20, colors="black", linewidths=0.3)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(rf"$\delta = {dr}$")

fig.suptitle(r"Edge-preserving remap: $\psi(x', y')$ with $g = x(1-y)$")
fig.tight_layout()
```

```python
# Check edge matching and corner psi
for dr in [0.05, 0.1, 0.2]:
    print(f"delta = {dr}:")
    # West edge
    psi_w = psi_remapped_ep(np.array([0.0]), np.array([0.5]), dr)[0]
    psi_w_orig = (1 - 0.0) * 0.5
    print(f"  west  (0, 0.5): psi={psi_w:.6f}  orig={psi_w_orig:.6f}  match={np.isclose(psi_w, psi_w_orig)}")
    # North edge
    psi_n = psi_remapped_ep(np.array([0.5]), np.array([1.0]), dr)[0]
    psi_n_orig = (1 - 0.5) * 1.0
    print(f"  north (0.5, 1): psi={psi_n:.6f}  orig={psi_n_orig:.6f}  match={np.isclose(psi_n, psi_n_orig)}")
    # Corner
    xc, yc = remap_edge_preserving(1.0, 0.0, dr)
    psi_c = psi_remapped_ep(np.array([1.0]), np.array([0.0]), dr)[0]
    print(f"  corner (1, 0) -> ({xc:.2f}, {yc:.2f}): psi={psi_c:.6f}")
```

The remapped $\psi$ has minimum value $\delta^2$ at the corner.
All streamlines from the free edges have $\psi \geq \delta^2$. The
region where $\psi_{\mathrm{remapped}} < \delta^2$ is unreachable
by advection — a **shadow zone** created by the rounded corner.

Particles starting on the south edge ($y = 0$, $x > 0$) cannot
creep into the shadow zone either: the remap gives $y' = \delta x >
0$, so the velocity there has $v = y' > 0$, lifting the particle off
the edge.

```python
dr = 0.1
psi_orig = (1 - X) * Y
psi_r = psi_remapped_ep(X, Y, dr)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Top row: full cell
for ax, (psi, title) in zip(axes[0], [
    (psi_orig, r"Original $\psi = (1-x)\,y$"),
    (psi_r, rf"Remapped ($\delta = {dr}$)"),
]):
    cs = ax.contourf(X, Y, psi, levels=20, cmap="viridis")
    ax.contour(X, Y, psi, levels=20, colors="black", linewidths=0.3)
    fig.colorbar(cs, ax=ax)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(title)

# Bottom row: zoom on corner
xz = np.linspace(0.6, 1, 300)
yz = np.linspace(0, 0.4, 300)
Xz, Yz = np.meshgrid(xz, yz)
psi_orig_z = (1 - Xz) * Yz
psi_r_z = psi_remapped_ep(Xz, Yz, dr)

lvls = np.linspace(0, 0.05, 25)
for ax, (psi, title) in zip(axes[1], [
    (psi_orig_z, "Original (zoom)"),
    (psi_r_z, "Remapped (zoom)"),
]):
    cs = ax.contourf(Xz, Yz, psi, levels=lvls, cmap="viridis")
    ax.contour(Xz, Yz, psi, levels=lvls, colors="black", linewidths=0.3)
    fig.colorbar(cs, ax=ax)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(title)

# Shadow zone boundary on the zoomed remapped panel
axes[1][1].contour(Xz, Yz, psi_r_z, levels=[dr**2], colors="red", linewidths=1.5)
axes[1][1].set_title(f"Remapped (zoom) — red: shadow zone ($\\psi = {dr**2}$)")

fig.tight_layout()
```

```python
# Velocity on the south edge (y=0): does the remap lift particles?
x_edge = np.linspace(0, 1, 100)
y_edge = np.zeros_like(x_edge)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

for ax, dr in zip(axes, [0.0, 0.1]):
    if dr > 0:
        u, v = vel_remap_ep(x_edge, y_edge, delta_r=dr)
    else:
        u, v = 1.0 - x_edge, y_edge
    ax.plot(x_edge, u, label="$u$")
    ax.plot(x_edge, v, label="$v$")
    ax.set_xlabel("$x$ along south edge ($y=0$)")
    ax.set_ylabel("velocity")
    ax.legend()
    ax.set_title(rf"$\delta = {dr}$")

fig.suptitle("Velocity on south edge: does the remap lift particles off $y=0$?")
fig.tight_layout()
```

```python
# Trajectories: edge-preserving remap
def vfn_remap_ep(x, y):
    return vel_remap_ep(x, y, delta_r=0.1)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

configs_ep = [
    ("Uncorrected", make_velocity_fn(None, 0, 0)),
    ("Source", make_velocity_fn(vel_source, 0.01, 0.05)),
    ("Remap (edge-preserving)", vfn_remap_ep),
]

for ax, (label, vfn) in zip(axes, configs_ep):
    psi = (1 - X) * Y
    ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)
    for t in t_starts:
        traj = integrate_vfn(vfn, t, 1.0 - t, 0.002, 50000)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.5)
        ax.plot(t, 1.0 - t, ".", color="black", markersize=2)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.suptitle("First cell: diagonal starts")
fig.tight_layout()
```

```python
# Second cell with edge-preserving remap
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

for ax, (label, vfn) in zip(axes, configs_ep):
    psi = (1 - X) * Y
    ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)
    for x0, y0 in starts_2nd:
        traj = integrate_vfn(vfn, x0, y0, 0.002, 50000)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.5)
        ax.plot(x0, y0, ".", color="black", markersize=2)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.suptitle("Second cell: uncorrected exits → uncorrected / source / remap")
fig.tight_layout()
```

## The horizon construction

Instead of remapping coordinates, subtract a smooth function that
vanishes on the free edges from $\psi$:

$$\psi_{\mathrm{new}}(x, y) = (1-x)\,y - \varepsilon\,x\,(1-y)$$

The second term $x(1-y)$ is the "complementary" stream function:
maximal at the concave corner $(1,0)$, zero on the free edges.

Properties:
- **West edge** ($x=0$): $\psi_{\mathrm{new}} = y$ — unchanged.
- **North edge** ($y=1$): $\psi_{\mathrm{new}} = 1-x$ — unchanged.
- **Corner** $(1,0)$: $\psi_{\mathrm{new}} = -\varepsilon < 0$.
- **Horizon** ($\psi = 0$): smooth curve from $(0,0)$ to $(1,1)$,
  bulging away from the corner. Beyond it, $\psi < 0$.

The velocity field is:

$$u = \frac{\partial \psi_{\mathrm{new}}}{\partial y} = (1-x) + \varepsilon\,x = 1 - x(1-\varepsilon)$$
$$v = -\frac{\partial \psi_{\mathrm{new}}}{\partial x} = y + \varepsilon(1-y) = \varepsilon + y(1-\varepsilon)$$

This is the standard C-grid interpolation with **the land-edge values
floored at $\varepsilon$ instead of zero**: $U_1 = \varepsilon$,
$V_0 = \varepsilon$. Continuity holds: $\varepsilon - 1 + 1 -
\varepsilon = 0$. No stagnation point. Divergence-free. One
parameter.

```python
def psi_horizon(X, Y, eps):
    return (1 - X) * Y - eps * X * (1 - Y)


def vel_horizon(x, y, eps):
    return 1 - x * (1 - eps), eps + y * (1 - eps)


fig, axes = plt.subplots(2, 4, figsize=(16, 7))

vlim = 0.5
lvls = np.linspace(-vlim, vlim, 30)

for col, eps in enumerate([0.0, 0.02, 0.05, 0.1]):
    # Top: stream function with horizon
    ax = axes[0, col]
    psi = psi_horizon(X, Y, eps)
    cf = ax.contourf(X, Y, psi, levels=lvls, cmap="RdBu_r", extend="both")
    ax.contour(X, Y, psi, levels=lvls, colors="black", linewidths=0.3)
    if eps > 0:
        ax.contour(X, Y, psi, levels=[0], colors="red", linewidths=1.5)
    fig.colorbar(cf, ax=ax)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(rf"$\varepsilon = {eps}$")

    # Bottom: trajectories
    ax = axes[1, col]
    ax.contour(X, Y, psi, levels=20, colors="0.85", linewidths=0.3)
    if eps > 0:
        ax.contour(X, Y, psi, levels=[0], colors="red", linewidths=0.8)

    def vfn(x, y, _eps=eps):
        return vel_horizon(x, y, _eps)

    for t in t_starts:
        traj = integrate_vfn(vfn, t, 1.0 - t, 0.002, 50000)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.4)
        ax.plot(t, 1.0 - t, ".", color="black", markersize=2)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")

axes[0, 0].set_title(r"Original ($\varepsilon = 0$)")
fig.suptitle(r"Horizon construction: $\psi = (1-x)y - \varepsilon\,x(1-y)$" +
             "\nred = horizon ($\\psi = 0$), blue = shadow zone ($\\psi < 0$)")
fig.tight_layout()
```

```python
# Second cell test
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

horizon_configs = [
    ("Uncorrected", lambda x, y: (1.0 - x, y)),
    ("Source", lambda x, y: vel_source_wrap(x, y)),
    (r"Horizon ($\varepsilon=0.05$)", lambda x, y: vel_horizon(x, y, 0.05)),
]

def vel_source_wrap(x, y):
    u, v = 1.0 - x, y
    r2 = (x - 1.0)**2 + y**2 + 0.05**2
    return u + 0.01 * (x - 1.0) / r2, v + 0.01 * y / r2

for ax, (label, vfn) in zip(axes, horizon_configs):
    psi = (1 - X) * Y
    ax.contour(X, Y, psi, levels=30, colors="0.85", linewidths=0.3)
    for x0, y0 in starts_2nd:
        traj = integrate_vfn(vfn, x0, y0, 0.002, 50000)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.5)
        ax.plot(x0, y0, ".", color="black", markersize=2)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.suptitle("Second cell: uncorrected exits → uncorrected / source / horizon")
fig.tight_layout()
```

## Shadow zone repulsion

In the shadow zone ($\psi < 0$), we override the velocity with
$\nabla\psi$, which points perpendicular to streamlines toward the
horizon. This is not divergence-free, but the shadow zone is virtual
land — no physical flow there. It serves as a safety net for
numerical overshoots.

```python
def vel_horizon_with_shadow(x, y, eps, alpha=10.0):
    """Horizon velocity + additive shadow zone repulsion."""
    u = 1 - x * (1 - eps)
    v = eps + y * (1 - eps)
    psi = (1 - x) * y - eps * x * (1 - y)
    if psi < 0:
        # Add grad(psi) scaled by deficit — particle still follows
        # the flow but gets nudged back toward the horizon
        u += -alpha * (-psi) * (eps + y * (1 - eps))
        v += alpha * (-psi) * (1 - x * (1 - eps))
    return u, v


# Seed particles in the shadow zone
eps = 0.05
# Horizon on diagonal: t = 1/(1+sqrt(eps)), so shadow zone is t > that
t_horizon = 1 / (1 + np.sqrt(eps))
shadow_starts = [
    (0.95, 0.03), (0.98, 0.01), (0.99, 0.005),
    (0.90, 0.08), (0.85, 0.05), (0.92, 0.02),
    (1.0 - 0.01, 0.01), (0.97, 0.04),
]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

configs_shadow = [
    ("Horizon only", lambda x, y: vel_horizon(x, y, eps)),
    ("Horizon + shadow repulsion", lambda x, y: vel_horizon_with_shadow(x, y, eps)),
]

for ax, (label, vfn) in zip(axes, configs_shadow):
    psi = psi_horizon(X, Y, eps)
    ax.contourf(X, Y, psi, levels=np.linspace(-0.5, 0.5, 30),
                cmap="RdBu_r", extend="both")
    ax.contour(X, Y, psi, levels=[0], colors="red", linewidths=1)

    for x0, y0 in shadow_starts:
        traj = integrate_vfn(vfn, x0, y0, 0.002, 50000)
        ax.plot(traj[:, 0], traj[:, 1], linewidth=0.8)
        ax.plot(x0, y0, "o", color="black", markersize=3)

    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(label)

fig.suptitle("Particles seeded in shadow zone")
fig.tight_layout()
```

## Transit time correction

The $\varepsilon$-floor speeds up transit slightly. Both ODEs are
linear, so transit times are analytical:

$$T_{\mathrm{orig}} = -\ln y_0, \quad T_{\mathrm{corr}} = \frac{1}{1-\varepsilon} \ln\frac{1}{(1-\varepsilon)\,y_0 + \varepsilon}$$

where $y_0 = \psi$ is the streamline's value at the west edge.

```python
def transit_time_orig(y0):
    return -np.log(y0)

def transit_time_corr(y0, eps):
    a = 1 - eps
    return (1 / a) * np.log(1 / (a * y0 + eps))

y0 = np.linspace(0.001, 0.99, 200)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(y0, transit_time_orig(y0), label="Original")
for eps in [0.02, 0.05, 0.1]:
    ax.plot(y0, transit_time_corr(y0, eps), label=rf"$\varepsilon={eps}$")
ax.set_xlabel("$y_0$ (streamline $\\psi$ value)")
ax.set_ylabel("Transit time")
ax.legend()
ax.set_title("Transit time west → north")

ax = axes[1]
for eps in [0.02, 0.05, 0.1]:
    ratio = transit_time_corr(y0, eps) / transit_time_orig(y0)
    ax.plot(y0, ratio, label=rf"$\varepsilon={eps}$")
ax.axhline(1, color="black", linewidth=0.5)
ax.set_xlabel("$y_0$")
ax.set_ylabel("$T_{\\mathrm{corr}} / T_{\\mathrm{orig}}$")
ax.legend()
ax.set_title("Transit time ratio")
ax.set_ylim(0.8, 1.05)

fig.tight_layout()
```

For most streamlines ($y_0 > 0.1$), the speedup is $O(\varepsilon)$
— a few percent. For near-corner streamlines ($y_0 \to 0$), the
original transit time diverges while the corrected one stays finite:
$T_{\mathrm{corr}} \to \frac{1}{1-\varepsilon}\ln\frac{1}{\varepsilon}$.
This is the desired effect — making trapped particles transit in
finite time.

## Transport budget

The stream function values at the cell corners determine the flux
through each edge. With $U_0 = 1$:

```python
print("Stream function at cell corners and edge fluxes:\n")
print(f"  {'':20s} {'Original':>10s} {'eps=0.02':>10s} {'eps=0.05':>10s} {'eps=0.1':>10s}")
for label, x, y in [("SW (0,0)", 0, 0), ("SE (1,0)", 1, 0), ("NE (1,1)", 1, 1), ("NW (0,1)", 0, 1)]:
    vals = []
    for eps in [0, 0.02, 0.05, 0.1]:
        psi = (1-x)*y - eps*x*(1-y)
        vals.append(psi)
    print(f"  psi{label:15s} {vals[0]:+10.4f} {vals[1]:+10.4f} {vals[2]:+10.4f} {vals[3]:+10.4f}")

print()
print(f"  {'Edge flux':20s} {'Original':>10s} {'eps=0.02':>10s} {'eps=0.05':>10s} {'eps=0.1':>10s}")
for eps in [0, 0.02, 0.05, 0.1]:
    sw = 0
    se = -eps
    ne = 0
    nw = 1
    q_west = nw - sw     # inflow
    q_north = nw - ne    # outflow
    q_east = ne - se      # through land
    q_south = sw - se     # through land
    if eps == 0:
        print(f"  {'West (inflow)':<20s} {q_west:+10.4f}")
        print(f"  {'North (outflow)':<20s} {q_north:+10.4f}")
        print(f"  {'East (land)':<20s} {q_east:+10.4f}")
        print(f"  {'South (land)':<20s} {q_south:+10.4f}")
    else:
        pass

# Print as table
print()
for label, calc in [
    ("West (inflow)", lambda e: 1),
    ("North (outflow)", lambda e: 1),
    ("East → land", lambda e: e),
    ("South ← land", lambda e: e),
]:
    vals = [calc(e) for e in [0, 0.02, 0.05, 0.1]]
    print(f"  {label:<20s} {vals[0]:10.4f} {vals[1]:10.4f} {vals[2]:10.4f} {vals[3]:10.4f}")
```

The west→north transport is exactly $U_0 = 1$ regardless of
$\varepsilon$. The $\varepsilon$-floor creates a small parasitic
circulation: $\varepsilon$ units enter from land at the south and
$\varepsilon$ units leak to land at the east. This is the price of
eliminating the stagnation point. For $\varepsilon = 0.02$, the
parasitic flux is 2% of the main transport.

## Discussion

The horizon construction reduces to a remarkably simple prescription:
**floor the land-edge velocities at $\varepsilon$ instead of zero**.
This is equivalent to subtracting $\varepsilon\,x(1-y)$ from the
stream function, which:

- Preserves $\psi$ exactly on the free edges
- Creates a smooth horizon ($\psi = 0$ curve) from $(0,0)$ to $(1,1)$
- Makes $\psi < 0$ in the shadow zone beyond the horizon
- Eliminates the stagnation point entirely ($u > 0$ and $v > 0$ everywhere)
- Is divergence-free (continuity still holds)
- Has a single free parameter $\varepsilon$

The trade-off: the velocity on the free edges changes slightly.
$V_0$ goes from $0$ to $\varepsilon$ (small northward leak at the
south boundary) and $U_1$ from $0$ to $\varepsilon$ (small eastward
leak at the east boundary). In the C-grid, these edges are shared
with neighboring cells — the correction must be applied inside the
interpolation, not by modifying the stored edge values.: $\varepsilon$ (strength)
and $\delta$ (regularization radius). Key properties:

- **Divergence-free** by construction (radial field from a harmonic
  potential).
- **Correct direction**: pushes away from both land edges
  ($u < 0$ near east edge, $v > 0$ near south edge).
- **Always dominates** near the corner ($1/r$ vs ambient $\sim r$),
  with a smooth cap at $r < \delta$.
- **Small far away**: at the opposite corner, the correction is
  $O(\varepsilon)$ compared to ambient $O(1)$.
- **One correction per cell**: a cell can have at most one concave
  corner.
- **Cell-edge discontinuity**: the correction is applied per-cell, so
  there is a jump at cell boundaries. For small $\varepsilon$, this
  is small.
- **Implementation**: in `spatial_interpolation_UV_c_grid` in
  `parcels.h`, add $\varepsilon \, (x-1) / (r^2 + \delta^2)$ to $U$
  and $\varepsilon \, y / (r^2 + \delta^2)$ to $V$ when the corner
  edge values are zero.
