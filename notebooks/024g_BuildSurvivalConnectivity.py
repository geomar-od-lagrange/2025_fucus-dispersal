# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: py:percent,md,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Build survival-weighted connectivity store
#
# The subbasin→subbasin residence connectivity of `024c`, with beaching
# **progressively removed**. `024c` re-partitions the 024 counts store and
# weights every `(trajectory, obs)` particle-timestep equally; here each
# timestep is weighted by the particle's surviving (un-beached) fraction at
# that age,
#
# ```
# S(t) = exp(−A(t)),   A(t) = cumsum(Δt/τ) over in-band steps
# ```
#
# — the same two-state near-shore rate as `024e_BuildBeaching` /
# `024f_BuildSurvivalOccupancy` (see [beaching.md](../docs/beaching.md)).
# This is 024f's survival weighting applied to 024c's origin→target axis:
# 024f asks *where the still-drifting weight is*, this asks *which subbasin
# it came from and which it is in*.
#
# Emits **two** weights per bin so the comparison in `028` is self-consistent
# on an identical population: `n_obs` (plain, every sample = 1) and `w_obs`
# (Σ S) over the same window, hexing and subbasin assignment.
#
# **Reads the beaching-forcing sidecar** built by `024d`, never the counts
# store and never the trajectory zarrs — the per-(particle, hour) ingredients
# of the rate (`w_on`, `dist`, `hex`, `flat`) plus `release_hex` are cached
# there, so a rate-model member is a seconds-per-zarr numpy pass. Aggregates
# with `np.bincount` over the dense `(origin, target, age_bin)` index.
# Partitioned per `(regime, year, month, member)`; `028` pools.
#
# **Population caveat.** The sidecar holds water-seeded ("real") drifters
# only — 024d drops the zero-first-step-displacement particles. `024`/`024c`
# instead keep land-seeded particles as `-1` sentinel rows. So a survconn
# partition and the matching 024c partition do **not** describe the same
# particle set; `n_obs` here is the unweighted count over the real drifters
# alone, which is why it is emitted rather than read back from 024c.

# %%
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# %%
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


# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Write root; also holds the BeachingForcing sidecar this reads.
output_root = "../output"

regime = "surface_stokes"
release_year = 2019
# Restrict to releases in this calendar month (1..12); 0 = whole year.
release_month = 8

hex_radius = 6000

# Connectivity horizon (days): how far along each trajectory to accumulate
# residence. Must be ≥ the largest matrix horizon in 028 and ≤ the sidecar's
# `window_days`.
connectivity_max_days = 120
# Age-bin granularity (days); matches the counts / 024c store.
age_bin_days = 10
# Sidecar cadence (hours).
output_dt_hours = 1

# Rate model: a two-state hazard. Nothing happens below the onshore-Stokes
# threshold `w_c`, and in-band particles strand on a `tau_strong_hours`
# timescale above it, plus an optional forcing-independent background
# `tau_calm_days`.

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
# term is wired but currently expresses nothing.
trap_flat = 1.0
trap_wall = 1.0

# %% [markdown]
# # Derived layout / key + member tag

# %%
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
survconn_path = (
    store_root
    / f"HexAgg_survconn_r{hex_radius}m_{regime}_{release_year}{month_suffix}_{member}.parquet"
)
print(f"survconn → {survconn_path}")

# %% [markdown]
# # Hex → subbasin lookup
#
# Same sentinel convention as `024c`: two unnamed states fold to the single
# id `-1` — hexes absent from the key (land-seed / NaN, `hex = -1` in the
# sidecar) and in-key hexes outside every named polygon
# (`helcom_subbasin == -1`, the `_outside` category). The bincount needs a
# contiguous index, so the subbasin ids are mapped to `0..n_sub-1` with the
# `-1` sentinel taking slot 0.

# %%
key = pd.read_parquet(key_path, columns=["hex_id", "helcom_subbasin"])
hex_sub = key["helcom_subbasin"].fillna(-1).astype(int)
sub_ids = np.array(sorted(set(hex_sub.unique().tolist()) | {-1}))
assert sub_ids[0] == -1, sub_ids[:3]
n_sub = len(sub_ids)
subid_to_idx = pd.Series(np.arange(n_sub), index=sub_ids)

# hex_id → contiguous subbasin index; anything not in the key reindexes to
# NaN below and is filled to 0 (the -1 sentinel slot).
hexid_to_subidx = pd.Series(
    subid_to_idx.reindex(hex_sub.to_numpy()).to_numpy(), index=key["hex_id"].astype(int)
)

# The sidecar obs axis is read as hours throughout — the window slice, the
# age-bin divisor and the per-step exponent all convert with a literal 24 or
# with `output_dt_hours` alone — so any other cadence would silently mis-bin
# every age.
assert output_dt_hours == 1, (
    f"output_dt_hours={output_dt_hours}: the obs↔time conversion here is "
    f"hardcoded hourly"
)
assert connectivity_max_days >= age_bin_days, (
    f"connectivity_max_days={connectivity_max_days} d is shorter than one age "
    f"bin ({age_bin_days} d)"
)
n_agebin = -(-connectivity_max_days // age_bin_days)
print(f"{n_sub} subbasin slots (incl. -1 sentinel), {n_agebin} age bins "
      f"× {age_bin_days} d")

sidecar_dir = output_root / f"BeachingForcing/{regime}/{release_year}"


def to_subidx(hex_values):
    """Map `024a` hex ids to the contiguous subbasin index; unknown → 0 (-1)."""
    idx = hexid_to_subidx.reindex(np.asarray(hex_values).ravel()).to_numpy()
    return np.where(np.isnan(idx), 0, idx).astype(np.int64)


# %% [markdown]
# # Sidecar → survival-weighted connectivity
#
# Per real drifter, weight every in-window sample by `S = exp(−A)`; `A`
# accrues only inside the near-shore band. Aggregate `n_obs` (=1) and `w_obs`
# (=S) over `(origin_subbasin, target_subbasin, age_bin)` with `np.bincount`.
#
# The rate block below is the same ~15 lines as in `024e_BuildBeaching` and
# `024f_BuildSurvivalOccupancy` — deliberately duplicated rather than shared,
# as notebooks own their utilities; the reducers are held equal by
# reproducing the same store.

# %%
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


# %%
def connectivity_one_sidecar(path, release_doy):
    sc = xr.open_zarr(path)
    assert sc.attrs["band_max_m"] >= band_m, (sc.attrs["band_max_m"], band_m)
    assert sc.attrs["window_days"] >= connectivity_max_days, (
        sc.attrs["window_days"], connectivity_max_days
    )
    assert sc.attrs["hex_radius"] == hex_radius, (sc.attrs["hex_radius"], hex_radius)
    assert sc.sizes["obs"] >= connectivity_max_days * 24, (
        sc.sizes["obs"], connectivity_max_days * 24
    )
    sc = sc.isel(obs=slice(0, connectivity_max_days * 24))
    w_on = sc.w_on.values
    dist = sc.dist.values
    flat = sc.flat.values
    hex_at = sc.hex.values
    release_hex = sc.release_hex.values
    ntraj, nobs = w_on.shape

    _, _, a = beaching_exponent(w_on, dist, flat)
    surv = np.exp(-np.cumsum(a, axis=1)).astype("float32")

    origin_idx = to_subidx(release_hex)
    target_idx = to_subidx(hex_at).reshape(ntraj, nobs)
    age_bin = (np.arange(nobs) // (age_bin_days * 24)).astype(np.int64)

    # Row-major (origin, target, age_bin) key; origin is per trajectory, so it
    # broadcasts along obs.
    key_arr = (
        (origin_idx[:, None] * n_sub + target_idx) * n_agebin + age_bin[None, :]
    ).ravel()
    length = n_sub * n_sub * n_agebin
    assert key_arr.max(initial=-1) < length, (key_arr.max(), length)
    n = np.bincount(key_arr, minlength=length).astype(np.float64)
    w = np.bincount(
        key_arr, weights=surv.ravel().astype(np.float64), minlength=length
    )
    n = n.reshape(n_sub, n_sub, n_agebin)
    w = w.reshape(n_sub, n_sub, n_agebin)

    oi, ti, bi = np.nonzero(n > 0)
    frame = pd.DataFrame({
        "origin_subbasin": sub_ids[oi],
        "target_subbasin": sub_ids[ti],
        "release_doy": release_doy,
        "age_bin": bi.astype(np.int64),
        "n_obs": n[oi, ti, bi],
        "w_obs": w[oi, ti, bi],
    })
    return frame, ntraj, nobs


# %%
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

# %%
t0 = time.time()
frames = []
total_samples = 0
total_real = 0
for p, ts, _ in parsed:
    tz = time.time()
    frame, ntraj, nobs = connectivity_one_sidecar(p, int(ts.dayofyear))
    frames.append(frame)
    total_samples += ntraj * nobs
    total_real += ntraj
    n = frame["n_obs"].sum()
    w = frame["w_obs"].sum()
    print(f"  {p.name}: n_obs {n:,.0f}, w_obs {w:,.0f} "
          f"(drifting {w / max(n, 1):.3f}) [{time.time() - tz:.1f}s]")

survconn = (
    pd.concat(frames, ignore_index=True)
    .groupby(
        ["origin_subbasin", "target_subbasin", "release_doy", "age_bin"],
        as_index=False,
    )[["n_obs", "w_obs"]]
    .sum()
)
print(f"computed {len(survconn):,} rows in {time.time() - t0:.1f}s")

# %% [markdown]
# # Validation
#
# Checked before the write, so a mismatched store never lands on disk.

# %%
# Survival weight can never exceed the unweighted count.
assert (survconn["w_obs"] <= survconn["n_obs"] + 1e-6).all(), (
    "w_obs exceeds n_obs somewhere"
)

# Every subbasin id must be one the key knows (or the -1 sentinel).
unseen = (
    set(survconn["origin_subbasin"]) | set(survconn["target_subbasin"])
) - set(sub_ids.tolist())
if unseen:
    raise ValueError(f"subbasin ids not in the key: {sorted(unseen)}")

# Conservation: nothing is dropped or duplicated — `-1` targets are kept as
# sentinel rows (as in 024c), so the total is exactly every (trajectory, obs)
# sample of every real drifter in the window.
n_total = float(survconn["n_obs"].sum())
if abs(n_total - total_samples) > 1e-6 * max(total_samples, 1):
    raise ValueError(
        f"n_obs total {n_total:,.0f} != sample count {total_samples:,.0f}"
    )

# The first age bin holds exactly `age_bin_days * 24` samples per released
# real drifter — the per-bin bound, exact because the window is a whole
# number of bins.
bin0 = float(survconn.loc[survconn["age_bin"] == 0, "n_obs"].sum())
expected_bin0 = float(total_real * age_bin_days * 24)
if abs(bin0 - expected_bin0) > 1e-6 * max(expected_bin0, 1):
    raise ValueError(
        f"age_bin 0 n_obs {bin0:,.0f} != released real drifters × obs/bin "
        f"{expected_bin0:,.0f}"
    )
print(
    f"w_obs <= n_obs everywhere; n_obs total {n_total:,.0f} matches the "
    f"sample count; age_bin 0 holds {bin0:,.0f} = {total_real:,} real "
    f"drifters × {age_bin_days * 24} obs."
)

# %%
survconn.to_parquet(survconn_path)
print(f"wrote {survconn_path} ({survconn_path.stat().st_size / 1e6:.2f} MB)")

# %%
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m, member={member}")
print(f"  params: connectivity_max_days={connectivity_max_days}, "
      f"band_m={band_m:g}, w_c={w_c:g}, "
      f"tau_strong_hours={tau_strong_hours:g}, "
      f"tau_calm_days={tau_calm_days:g}, delta={delta:g}, "
      f"trap_flat/wall={trap_flat:g}/{trap_wall:g}"
      + (" (degenerate — shore type inert)" if trap_flat == trap_wall else ""))

named = survconn[
    (survconn["origin_subbasin"] >= 0) & (survconn["target_subbasin"] >= 0)
]
diag = named[named["origin_subbasin"] == named["target_subbasin"]]
print(f"  named-subbasin share of n_obs: "
      f"{named['n_obs'].sum() / n_total:.3f}")
print(f"  within-subbasin diagonal share (named): "
      f"{diag['n_obs'].sum() / max(named['n_obs'].sum(), 1):.3f} unweighted, "
      f"{diag['w_obs'].sum() / max(named['w_obs'].sum(), 1):.3f} weighted")
per_bin = survconn.groupby("age_bin")[["n_obs", "w_obs"]].sum()
per_bin["drifting"] = per_bin["w_obs"] / per_bin["n_obs"]
print(per_bin.to_string(float_format=lambda v: f"{v:,.3f}"))
