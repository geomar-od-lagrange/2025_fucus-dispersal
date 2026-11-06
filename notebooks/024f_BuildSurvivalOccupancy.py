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
# # Build survival-weighted occupancy store
#
# The free-drifting density with beaching **progressively removed**. Standard
# occupancy (`024` counts) weights every `(trajectory, obs)` sample equally;
# here each sample is weighted by the particle's **surviving (un-beached)
# fraction** at that age,
#
# ```
# S(t) = exp(−A(t)),   A(t) = cumsum(Δt/τ) over in-band steps
# ```
#
# — the same near-shore beaching rate as `024e_BuildBeaching`
# (see [beaching.md](../docs/beaching.md)). `A` only grows inside the
# near-shore band, so open-water residence is undiluted; a particle that
# lingers in a wave-exposed near-shore band loses weight fast. This is the
# deterministic occupancy analogue of 024e's fractional stranding: 024e
# records where the weight *leaves* (`beach_hex`), this records where the
# still-drifting weight *is* (`target_hex`).
#
# Emits **two** weights per bin so the fold-in is a self-consistent
# comparison: `occ` (plain, every sample = 1) and `surv` (Σ S) over the same
# window and hexing. `030_SurvivalHeatmaps` maps `surv`, `occ`, and the
# surviving fraction `surv/occ`.
#
# **Origin axis and exact denominator.** Rows are keyed by
# `(release_doy, age_bin, origin_subbasin, target_hex)`: `origin_subbasin` is
# the subbasin of the drifter's release hex (exactly as `024g` assigns it),
# so per-origin survival maps need no second aggregation. Samples with no
# in-key position (`hex = -1` in the sidecar — the drifter left the BSH domain
# or was deleted) are **kept** as `target_hex = -1` sentinel rows, as `024g`
# keeps its `-1` targets: every released drifter then contributes exactly
# `age_bin_days · 24` samples to every full age bin, so `Σocc` per
# `(release_doy, age_bin)` is the released water-seeded drifter count times
# the samples per bin, recoverable from the store alone.
#
# **No lifetime is applied.** `max_age_days` is a reduction window over the
# sidecar's `obs` axis; the store stays resolved per `age_bin`, so a
# biological survival `L(t)` can be multiplied in downstream.
#
# **Reads the beaching-forcing sidecar** built by `024d`, never the
# trajectory zarrs: the per-(particle, hour) ingredients of the rate
# (`w_on`, `dist`, `hex`, `ff`) are cached there, so every rate-model
# member is a seconds-per-zarr numpy pass with no raster, no Stokes and no
# trajectory I/O. Aggregates with `np.bincount` (occupancy is dense — every
# obs, every hex — unlike 024e's sparse in-band deposits). Partitioned per
# `(regime, year, month)`; `030` pools. The `024a` key is a hard
# prerequisite (its hex-id set defines the dense bincount index, its
# `helcom_subbasin` the origin axis).
#
# **Land-seeded particles are excluded.** The sidecar holds water-seeded
# ("real") drifters only — 024d drops the zero-first-step-displacement
# particles, exactly as `024`/`024b` do — so the occupancy sums here run over
# real drifters alone, as they did when this notebook read the zarrs itself.

# %%
import os
import re
import subprocess
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

# Reduction window (days): the half-open `[0, max_age_days)` hours of the
# sidecar's `obs` axis that are accumulated — the same window as 024e/024g.
# Must be ≥ the largest map horizon in 030 and ≤ the sidecar's `window_days`.
max_age_days = 180
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

# Shore-type retention weights: the strong-wave rate is scaled by
# `trap = trap_wall + (trap_flat - trap_wall) * ff/100`, interpolating between
# a fully hard-walled and a fully flat cell along the sidecar's flat fraction
# `ff`. The production point is 1.0/0.0 — retention proportional to the flat
# fraction, nothing retained on pure-wall shore. Equal values make the shore
# type inert.
trap_flat = 1.0
trap_wall = 0.0

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
    # Only a shore-type-sensitive member carries the trap weights, so the
    # 1.0/1.0 production tag stays as it is.
    + (f"_tf{trap_flat:g}_tw{trap_wall:g}"
       if (trap_flat != 1.0 or trap_wall != 1.0) else "")
)
member = member.replace(".", "p")
print(f"member: {member}")

month_suffix = f"_m{release_month:02d}" if release_month else ""
survocc_path = (
    store_root
    / f"HexAgg_survocc_r{hex_radius}m_{regime}_{release_year}{month_suffix}_{member}.parquet"
)
print(f"survocc → {survocc_path}")

sidecar_dir = output_root / f"BeachingForcing/{regime}/{release_year}"

# Provenance for the store attrs: which commit and which allocation wrote it.
# Empty slurm_job_id outside SLURM.
try:
    git_sha = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
except Exception:
    git_sha = "unknown"
slurm_job_id = os.environ.get("SLURM_JOB_ID", "")
print(f"git_sha={git_sha}, slurm_job_id={slurm_job_id or '-'}")

# %% [markdown]
# # Hex and subbasin indices
#
# The bincount needs contiguous indices on both keyed axes, each with a `-1`
# sentinel in slot 0. **Target:** `target_hex` slot 0 is `-1` (no in-key
# position), slots `1..n` the key's hex ids. **Origin:** the subbasin of the
# release hex, with `024g`'s convention — hexes absent from the key and
# in-key hexes outside every named polygon (`helcom_subbasin == -1`) both
# fold to `-1`.

# %%
key = pd.read_parquet(key_path, columns=["hex_id", "helcom_subbasin"])
key_hex_ids = key["hex_id"].astype(int).to_numpy()
assert -1 not in key_hex_ids, "key carries the -1 sentinel id"

# hex_id → 1..n_hex; anything not in the key reindexes to NaN and goes to 0.
hex_ids = np.concatenate([[-1], key_hex_ids])
n_hex = len(hex_ids)
hexid_to_idx = pd.Series(np.arange(1, n_hex), index=key_hex_ids)

hex_sub = key["helcom_subbasin"].fillna(-1).astype(int)
sub_ids = np.array(sorted(set(hex_sub.unique().tolist()) | {-1}))
assert sub_ids[0] == -1, sub_ids[:3]
n_sub = len(sub_ids)
subid_to_idx = pd.Series(np.arange(n_sub), index=sub_ids)
hexid_to_subidx = pd.Series(
    subid_to_idx.reindex(hex_sub.to_numpy()).to_numpy(), index=key_hex_ids
)


def to_hexidx(hex_values):
    """Map `024a` hex ids to the contiguous target index; unknown → 0 (-1)."""
    idx = hexid_to_idx.reindex(np.asarray(hex_values).ravel()).to_numpy()
    return np.where(np.isnan(idx), 0, idx).astype(np.int64)


def to_subidx(hex_values):
    """Map `024a` hex ids to the contiguous subbasin index; unknown → 0 (-1)."""
    idx = hexid_to_subidx.reindex(np.asarray(hex_values).ravel()).to_numpy()
    return np.where(np.isnan(idx), 0, idx).astype(np.int64)


# The sidecar obs axis is read as hours throughout (window slice, age-bin
# divisor, per-step exponent), so any other cadence would mis-bin every age.
assert output_dt_hours == 1, (
    f"output_dt_hours={output_dt_hours}: the obs↔time conversion here is "
    f"hardcoded hourly"
)
# Ceil, not floor: a partial trailing age bin still holds samples, and a
# floored count would fold them into the next hex's bin 0 in the
# row-major bincount key.
n_agebin = -(-max_age_days // age_bin_days)
print(f"{n_hex - 1:,} key hexes (+ -1 sentinel), {n_sub} subbasin slots "
      f"(incl. -1), {n_agebin} age bins × {age_bin_days} d")

# %% [markdown]
# # Sidecar → survival-weighted occupancy
#
# Per real drifter, weight every in-window occupancy sample by `S = exp(−A)`;
# `A` accrues only inside the near-shore band. Aggregate `occ` (=1) and `surv`
# (=S) over `(origin_subbasin, target_hex, age_bin)` with `np.bincount`,
# keeping the `target_hex = -1` samples.
#
# The rate block below is the same ~15 lines as in `024e_BuildBeaching` —
# deliberately duplicated rather than shared, as notebooks own their
# utilities; the reducers are held equal by reproducing the same store.

# %%
def beaching_exponent(w_on, dist, ff):
    """Per-step beaching exponent `a = Δt/τ` from the sidecar ingredients.

    Returns `(in_band, r, a)` — the band gate, the forcing ramp (the share of
    the strong-wave rate that is on), and the per-step exponent, all
    `(trajectory, obs)`.
    """
    in_band = (dist.astype("int16") * 25 < band_m) & (dist != 255)
    w = w_on.astype("float32")
    # Linear in the flat fraction, so a half-flat cell traps halfway between.
    trap = (
        trap_wall + (trap_flat - trap_wall) * ff.astype("float32") / np.float32(100.0)
    ).astype("float32")
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
def occupancy_one_sidecar(path, release_doy):
    sc = xr.open_zarr(path)
    assert sc.attrs["band_max_m"] >= band_m, (sc.attrs["band_max_m"], band_m)
    assert sc.attrs["window_days"] >= max_age_days, (
        sc.attrs["window_days"], max_age_days
    )
    assert sc.attrs["hex_radius"] == hex_radius, (sc.attrs["hex_radius"], hex_radius)
    assert sc.sizes["obs"] >= max_age_days * 24, (
        sc.sizes["obs"], max_age_days * 24
    )
    sc = sc.isel(obs=slice(0, max_age_days * 24))
    w_on = sc.w_on.values
    dist = sc.dist.values
    ff = sc.ff.values
    assert ff.min() >= 0 and ff.max() <= 100, (ff.min(), ff.max())
    hex_at = sc.hex.values
    release_hex = sc.release_hex.values
    ntraj, nobs = w_on.shape

    _, _, a = beaching_exponent(w_on, dist, ff)
    surv = np.exp(-np.cumsum(a, axis=1)).astype("float32")

    origin_idx = to_subidx(release_hex)
    target_idx = to_hexidx(hex_at).reshape(ntraj, nobs)
    age_bin = (np.arange(nobs) // (age_bin_days * 24)).astype(np.int64)

    # Row-major (origin, target, age_bin) key; origin is per trajectory, so it
    # broadcasts along obs. Every sample is kept — no validity mask.
    key_arr = (
        (origin_idx[:, None] * n_hex + target_idx) * n_agebin + age_bin[None, :]
    ).ravel()
    length = n_sub * n_hex * n_agebin
    assert key_arr.max(initial=-1) < length, (key_arr.max(), length)
    occ = np.bincount(key_arr, minlength=length).astype(np.float64)
    surv_agg = np.bincount(
        key_arr, weights=surv.ravel().astype(np.float64), minlength=length
    )
    occ = occ.reshape(n_sub, n_hex, n_agebin)
    surv_agg = surv_agg.reshape(n_sub, n_hex, n_agebin)

    oi, hi, bi = np.nonzero(occ > 0)
    frame = pd.DataFrame({
        "release_doy": release_doy,
        "age_bin": bi.astype(np.int64),
        "origin_subbasin": sub_ids[oi].astype(np.int64),
        "target_hex": hex_ids[hi].astype(np.int64),
        "occ": occ[oi, hi, bi],
        "surv": surv_agg[oi, hi, bi],
    })
    # Released drifters and samples, for the exact-denominator check below.
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
# (release_doy, released water-seeded drifters, obs per drifter) per sidecar.
released = []
for p, ts, _ in parsed:
    tz = time.time()
    frame, ntraj, nobs = occupancy_one_sidecar(p, int(ts.dayofyear))
    frames.append(frame)
    released.append((int(ts.dayofyear), ntraj, nobs))
    o = frame["occ"].sum()
    s = frame["surv"].sum()
    print(f"  {p.name}: occ {o:,.0f}, surv {s:,.0f} "
          f"(drifting {s / max(o, 1):.3f}) [{time.time() - tz:.1f}s]")

survocc = (
    pd.concat(frames, ignore_index=True)
    .groupby(["release_doy", "age_bin", "origin_subbasin", "target_hex"],
             as_index=False)[["occ", "surv"]].sum()
)
print(f"computed {len(survocc):,} rows in {time.time() - t0:.1f}s")

# %% [markdown]
# # Validation
#
# Every `target_hex` must be in the key or the `-1` sentinel (out-of-key
# positions carry `hex = -1` in the sidecar), and every `origin_subbasin` a
# subbasin the key knows — a stray id means the store was built against a
# different key. Checked before the write, so a mismatched store never lands
# on disk.

# %%
unseen = set(survocc["target_hex"]) - set(hex_ids.tolist())
if unseen:
    raise ValueError(
        f"{len(unseen)} target_hex not in {key_path.name}: "
        f"{sorted(unseen)[:10]} ..."
    )
unseen_sub = set(survocc["origin_subbasin"]) - set(sub_ids.tolist())
if unseen_sub:
    raise ValueError(f"origin_subbasin ids not in the key: {sorted(unseen_sub)}")
print("every target_hex is in the key (or -1), every origin_subbasin too.")

# %%
# Survival can never exceed occupancy.
assert (survocc["surv"] <= survocc["occ"] + 1e-6).all(), "surv exceeds occ somewhere"

# Exact denominator: with the -1 sentinel rows kept, Σocc per
# (release_doy, age_bin) is the released water-seeded drifters of that doy
# times the samples the bin holds — age_bin_days * 24 for a full bin, the
# remainder for a partial trailing one.
bin_hours = age_bin_days * 24
expected = (
    pd.DataFrame(
        [
            (doy, b, ntraj * min(bin_hours, nobs - b * bin_hours))
            for doy, ntraj, nobs in released
            for b in range(-(-nobs // bin_hours))
        ],
        columns=["release_doy", "age_bin", "expected"],
    )
    .groupby(["release_doy", "age_bin"])["expected"].sum()
)
got = survocc.groupby(["release_doy", "age_bin"])["occ"].sum()
assert got.index.equals(expected.index), (got.index, expected.index)
assert np.array_equal(got.to_numpy(), expected.to_numpy().astype(np.float64)), (
    (got - expected).abs().max()
)
n_released = sum(ntraj for _, ntraj, _ in released)
print(f"surv <= occ everywhere; Σocc per (release_doy, age_bin) = released "
      f"drifters × samples per bin exactly ({n_released:,} drifters, "
      f"{len(got)} (doy, bin) cells).")

# %% [markdown]
# # Write
#
# The window, the binning and the rate-model point travel with the partition
# as `df.attrs` (pandas stores them in the parquet metadata), so a consumer
# pooling partitions can assert they agree instead of trusting filenames.

# %%
survocc.attrs = {
    "builder": "024f_BuildSurvivalOccupancy",
    "regime": regime,
    "release_year": int(release_year),
    "release_month": int(release_month),
    "hex_radius": int(hex_radius),
    "member": member,
    "max_age_days": int(max_age_days),
    "age_bin_days": int(age_bin_days),
    "output_dt_hours": int(output_dt_hours),
    "w_c": float(w_c),
    "tau_strong_hours": float(tau_strong_hours),
    "tau_calm_days": float(tau_calm_days),
    "delta": float(delta),
    "trap_flat": float(trap_flat),
    "trap_wall": float(trap_wall),
    "band_m": float(band_m),
    "git_sha": git_sha,
    "slurm_job_id": slurm_job_id,
}
survocc.to_parquet(survocc_path)
print(f"wrote {survocc_path} ({survocc_path.stat().st_size / 1e6:.2f} MB)")

# %%
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m, member={member}")
print(f"  params: max_age_days={max_age_days}, band_m={band_m:g}, "
      f"w_c={w_c:g}, tau_strong_hours={tau_strong_hours:g}, "
      f"tau_calm_days={tau_calm_days:g}, delta={delta:g}, "
      f"trap_flat/wall={trap_flat:g}/{trap_wall:g}")
# Drifting fraction per age bin over all released drifters (sentinel rows
# included, so surv/occ is the mean survival of the whole release), and the
# share of samples with no in-key position.
per_bin = survocc.groupby("age_bin")[["occ", "surv"]].sum()
per_bin["drifting"] = per_bin["surv"] / per_bin["occ"]
per_bin["out_of_key"] = (
    survocc[survocc["target_hex"] == -1].groupby("age_bin")["occ"].sum()
    .reindex(per_bin.index, fill_value=0.0) / per_bin["occ"]
)
print(per_bin.to_string(float_format=lambda v: f"{v:,.3f}"))
print(f"  origin subbasins: {survocc['origin_subbasin'].nunique()} "
      f"({sorted(survocc['origin_subbasin'].unique().tolist())})")
