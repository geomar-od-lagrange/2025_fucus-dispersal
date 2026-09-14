---
jupyter:
  jupytext:
    cell_metadata_filter: tags,-all
    formats: py:percent,md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Build survival-weighted occupancy store

The free-drifting density with beaching **progressively removed**. Standard
occupancy (`024` counts) weights every `(trajectory, obs)` sample equally;
here each sample is weighted by the particle's **surviving (un-beached)
fraction** at that age,

```
S(t) = exp(−A(t)),   A(t) = cumsum(Δt/τ) over in-band steps
```

— the same near-shore beaching rate as `024e_BuildBeaching`
(see [beaching.md](../docs/beaching.md)). `A` only grows inside the
near-shore band, so open-water residence is undiluted; a particle that
lingers in a wave-exposed near-shore band loses weight fast. This is the
deterministic occupancy analogue of 024e's fractional stranding: 024e
records where the weight *leaves* (`beach_hex`), this records where the
still-drifting weight *is* (`target_hex`).

Emits **two** weights per bin so the fold-in is a self-consistent
comparison: `occ` (plain, every sample = 1) and `surv` (Σ S) over the same
window and hexing. `030_SurvivalHeatmaps` maps `surv`, `occ`, and the
surviving fraction `surv/occ`.

**Reads the beaching-forcing sidecar** built by `024d`, never the
trajectory zarrs: the per-(particle, hour) ingredients of the rate
(`w_on`, `dist`, `hex`, `flat`) are cached there, so every rate-model
member is a seconds-per-zarr numpy pass with no raster, no Stokes and no
trajectory I/O. Aggregates with `np.bincount` (occupancy is dense — every
obs, every hex — unlike 024e's sparse in-band deposits). Partitioned per
`(regime, year, month)`; `030` pools. The `024a` key is a hard
prerequisite (its hex-id set defines the dense bincount index).

**Land-seeded particles are excluded.** The sidecar holds water-seeded
("real") drifters only — 024d drops the zero-first-step-displacement
particles, exactly as `024`/`024b` do — so the occupancy sums here run over
real drifters alone, as they did when this notebook read the zarrs itself.

```python
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
```

```python
# Pattern: Fucus_BSH_YYYYMMDDTHHMMSS_{regime}_dt{N}min_seed{S}. The sidecar
# zarrs carry the stem of the trajectory zarr they were built from.
_ZARR_STEM_RE = re.compile(
    r"^Fucus_BSH_(\d{8}T\d{6})_(surface_stokes|surface|bottom)_dt\d+min_seed\d+$"
)


def parse_zarr_stem(path):
    """Parse a trajectory zarr filename into ``(release_time, regime)``."""
    m = _ZARR_STEM_RE.match(Path(path).stem)
    if m is None:
        raise ValueError(
            f"zarr filename does not match expected pattern: {Path(path).name!r}"
        )
    return pd.Timestamp(m.group(1)), m.group(2)
```

# Parameters

```python tags=["parameters"]
# Write root; also holds the BeachingForcing sidecar this reads.
output_root = "../output"

regime = "surface_stokes"
release_year = 2019
# Restrict to releases in this calendar month (1..12); 0 = whole year.
release_month = 8

hex_radius = 6000

# Occupancy horizon (days): how far along each trajectory to accumulate
# residence. Must be ≥ the largest map horizon in 030 and ≤ the sidecar's
# `window_days`. Decoupled from 024e's `max_float_days` — survival occupancy
# runs to the mapping horizon, not the viability cutoff (a Fucus lifetime
# L(t) would enter as an extra factor on S).
occupancy_max_days = 120
# Age-bin granularity (days); matches the counts store.
age_bin_days = 10
# Sidecar cadence (hours).
output_dt_hours = 1

# Rate model: a two-state hazard. Nothing happens below the onshore-Stokes
# threshold `w_c`, and in-band particles strand on a `tau_strong_hours`
# timescale above it, plus an optional forcing-independent background
# `tau_calm_days` — expressing "drifts for a year in ordinary conditions,
# strands within hours in strong waves".

# --- step-hazard parameters ---
# Onshore-Stokes threshold (m/s); THE sweep axis.
w_c = 0.1
# e-folding time (h) while the onshore forcing is above w_c.
tau_strong_hours = 3.0
# Background in-band timescale (d), independent of forcing; 0 turns it off.
tau_calm_days = 0.0
# Linear edge width (m/s) around w_c; 0 gives a hard step.
delta = 0.0

# Near-shore band width (m); must be ≤ the sidecar's `band_max_m`.
band_m = 2000.0

# Shore-type retention weights, DELIBERATELY DEGENERATE (both 1.0) — the trap
# term is wired but currently expresses nothing, so every shore beaches alike
# for a given wave forcing. The `flat`/`wall` classification is carried only
# as the seam for a future real substrate/exposure dataset.
trap_flat = 1.0
trap_wall = 1.0
```

# Derived layout / key + member tag

```python
output_root = Path(output_root)
store_root = output_root / "HexAggregates"
store_root.mkdir(parents=True, exist_ok=True)

key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
if not key_path.exists():
    raise FileNotFoundError(
        f"Key file missing — run 024a_BuildHexKey first.\n  expected: {key_path}"
    )

# The member tag names the rate-model point in the store filename. Consumers
# take it as an opaque string parameter and never parse it.
member = (
    f"step_wc{w_c:g}_ts{tau_strong_hours:g}"
    f"_tc{f'{tau_calm_days:g}' if tau_calm_days > 0 else 'inf'}"
    + (f"_d{delta:g}" if delta > 0 else "")
)
member = member.replace(".", "p")
print(f"member: {member}")

month_suffix = f"_m{release_month:02d}" if release_month else ""
survocc_path = (
    store_root
    / f"HexAgg_survocc_r{hex_radius}m_{regime}_{release_year}{month_suffix}_{member}.parquet"
)
print(f"survocc → {survocc_path}")

# Contiguous hex index for bincount: hex_id → 0..n-1 and back.
hex_ids = pd.read_parquet(key_path, columns=["hex_id"])["hex_id"].astype(int).to_numpy()
hexid_to_idx = pd.Series(np.arange(len(hex_ids)), index=hex_ids)
n_hex = len(hex_ids)
# Ceil, not floor: a partial trailing age bin still holds samples, and a
# floored count would fold them into the next hex's bin 0 in the
# row-major bincount key.
n_agebin = -(-occupancy_max_days // age_bin_days)
print(f"{n_hex:,} hexes, {n_agebin} age bins × {age_bin_days} d")

sidecar_dir = output_root / f"BeachingForcing/{regime}/{release_year}"
```

# Sidecar → survival-weighted occupancy

Per real drifter, weight every in-window occupancy sample by `S = exp(−A)`;
`A` accrues only inside the near-shore band. Aggregate `occ` (=1) and `surv`
(=S) over `(target_hex, age_bin)` with `np.bincount`.

The rate block below is the same ~15 lines as in `024e_BuildBeaching` —
deliberately duplicated rather than shared, as notebooks own their
utilities; the reducers are held equal by reproducing the same store.

```python
def beaching_exponent(w_on, dist, flat):
    """Per-step beaching exponent `a = Δt/τ` from the sidecar ingredients.

    Returns `(in_band, r, a)` — the band gate, the forcing ramp (the share of
    the strong-wave rate that is on), and the per-step exponent, all
    `(trajectory, obs)`.
    """
    in_band = (dist.astype("int16") * 25 < band_m) & (dist != 255)
    w = w_on.astype("float32")
    trap = np.where(flat, trap_flat, trap_wall).astype("float32")
    if delta > 0:
        r = np.clip((w - w_c + delta / 2) / delta, 0.0, 1.0).astype("float32")
    else:
        r = (w >= w_c).astype("float32")
    inv_tau_h = np.float32(
        1.0 / (tau_calm_days * 24.0) if tau_calm_days > 0 else 0.0
    ) + trap * r / np.float32(tau_strong_hours)
    a = np.where(in_band, output_dt_hours * inv_tau_h, np.float32(0.0))
    return in_band, r, a
```

```python
def occupancy_one_sidecar(path, release_doy):
    sc = xr.open_zarr(path)
    assert sc.attrs["band_max_m"] >= band_m, (sc.attrs["band_max_m"], band_m)
    assert sc.attrs["window_days"] >= occupancy_max_days, (
        sc.attrs["window_days"], occupancy_max_days
    )
    assert sc.attrs["hex_radius"] == hex_radius, (sc.attrs["hex_radius"], hex_radius)
    assert sc.sizes["obs"] >= occupancy_max_days * 24, (
        sc.sizes["obs"], occupancy_max_days * 24
    )
    sc = sc.isel(obs=slice(0, occupancy_max_days * 24))
    w_on = sc.w_on.values
    dist = sc.dist.values
    flat = sc.flat.values
    hex_at = sc.hex.values
    ntraj, nobs = w_on.shape

    _, _, a = beaching_exponent(w_on, dist, flat)
    surv = np.exp(-np.cumsum(a, axis=1)).astype("float32")

    hex_idx = hexid_to_idx.reindex(hex_at.ravel()).to_numpy()
    hex_idx = np.where(np.isnan(hex_idx), -1, hex_idx).astype(np.int64).reshape(
        ntraj, nobs
    )
    age_bin = (np.arange(nobs) // (age_bin_days * 24)).astype(np.int64)
    age_bin2d = np.broadcast_to(age_bin[None, :], (ntraj, nobs))
    valid = hex_idx >= 0

    key = hex_idx[valid] * n_agebin + age_bin2d[valid]
    length = n_hex * n_agebin
    assert key.max(initial=-1) < length, (key.max(), length)
    occ = np.bincount(key, minlength=length).astype(np.float64)
    surv_agg = np.bincount(key, weights=surv[valid].astype(np.float64), minlength=length)
    occ = occ.reshape(n_hex, n_agebin)
    surv_agg = surv_agg.reshape(n_hex, n_agebin)

    hi, bi = np.nonzero(occ > 0)
    frame = pd.DataFrame({
        "release_doy": release_doy,
        "age_bin": bi.astype(np.int64),
        "target_hex": hex_ids[hi],
        "occ": occ[hi, bi],
        "surv": surv_agg[hi, bi],
    })
    # Valid-sample count, for the cheap conservation check below.
    return frame, int(valid.sum())
```

```python
sidecars = sorted(sidecar_dir.glob("*.zarr"))
parsed = [(p, *parse_zarr_stem(p)) for p in sidecars]
for p, ts, fn_regime in parsed:
    assert ts.year == release_year, (ts, release_year, p)
    assert fn_regime == regime, (fn_regime, regime, p)
if release_month:
    parsed = [x for x in parsed if x[1].month == release_month]
if not parsed:
    raise FileNotFoundError(
        f"no sidecars at {sidecar_dir}/"
        + (f" for month {release_month}" if release_month else "")
        + " — run 024d_BuildBeachingForcing first"
    )
release_doys = sorted({int(ts.dayofyear) for _, ts, _ in parsed})
print(f"{len(parsed)} sidecars, release_doys "
      f"{release_doys[0]}..{release_doys[-1]} ({len(release_doys)} unique)")
```

```python
t0 = time.time()
frames = []
total_valid = 0
for p, ts, _ in parsed:
    tz = time.time()
    frame, n_valid = occupancy_one_sidecar(p, int(ts.dayofyear))
    frames.append(frame)
    total_valid += n_valid
    o = frame["occ"].sum()
    s = frame["surv"].sum()
    print(f"  {p.name}: occ {o:,.0f}, surv {s:,.0f} "
          f"(drifting {s / max(o, 1):.3f}) [{time.time() - tz:.1f}s]")

survocc = (
    pd.concat(frames, ignore_index=True)
    .groupby(["release_doy", "age_bin", "target_hex"], as_index=False)[["occ", "surv"]].sum()
)
print(f"computed {len(survocc):,} rows in {time.time() - t0:.1f}s")
```

# Validation

Every `target_hex` must be in the key — out-of-key positions carry `hex = -1`
in the sidecar and were dropped above, so a stray id means the store was
built against a different key. Checked before the write, so a mismatched
store never lands on disk.

```python
unseen = set(survocc["target_hex"]) - set(hex_ids.tolist()) - {-1}
if unseen:
    raise ValueError(
        f"{len(unseen)} target_hex not in {key_path.name}: "
        f"{sorted(unseen)[:10]} ..."
    )
print("every target_hex is in the key.")
```

```python
# Survival can never exceed occupancy, and the grouped occ total must match
# the valid-sample count tallied while looping the sidecars.
assert (survocc["surv"] <= survocc["occ"] + 1e-6).all(), "surv exceeds occ somewhere"
occ_total = float(survocc["occ"].sum())
if abs(occ_total - total_valid) > 1e-6 * max(total_valid, 1):
    raise ValueError(
        f"occ total {occ_total:,.0f} != valid-sample count {total_valid:,.0f}"
    )
print(f"surv <= occ everywhere; occ total {occ_total:,.0f} matches "
      f"the valid-sample count.")
```

```python
survocc.to_parquet(survocc_path)
print(f"wrote {survocc_path} ({survocc_path.stat().st_size / 1e6:.2f} MB)")
```

```python
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m, member={member}")
print(f"  params: occupancy_max_days={occupancy_max_days}, band_m={band_m:g}, "
      f"w_c={w_c:g}, tau_strong_hours={tau_strong_hours:g}, "
      f"tau_calm_days={tau_calm_days:g}, delta={delta:g}, "
      f"trap_flat/wall={trap_flat:g}/{trap_wall:g}"
      + (" (degenerate — shore type inert)" if trap_flat == trap_wall else ""))
per_bin = survocc.groupby("age_bin")[["occ", "surv"]].sum()
per_bin["drifting"] = per_bin["surv"] / per_bin["occ"]
print(per_bin.to_string(float_format=lambda v: f"{v:,.3f}"))
```
