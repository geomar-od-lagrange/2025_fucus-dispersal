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
# # Build beaching store partition (weighted)
#
# Reduce the **beaching-forcing sidecar** built by `024d` into the stranding
# store for one rate-model member (see
# [beaching_sidecar.md](../plans/done/beaching_sidecar.md)). Each real drifter
# carries unit surviving (free-drifting) weight; at every near-shore step a
# fraction of that weight *strands* at the current coastal hex and leaves the
# drifting pool. This is the **fractional / weighted** scheme — the
# deterministic expectation of a per-step first-stranding process, with no
# Monte-Carlo noise. Writes a flat
# `(release_hex, release_doy, beach_hex, beach_age_bin, shore_type, disp_bin)
# → weight` table — additive across `release_doy`/year like the other `024x`
# stores.
#
# **Why weighted, not stochastic.** Releases carry only ~100 particles per
# seeding cell, so a per-particle random-strand rule gives noisy coverage at
# the high-age tail exactly where the free-drifting/beached split matters.
# The weighted deposit `dep[t,h] = S(t,h)·(1−e^{−Δt/τ})` equals the
# probability that a first-stranding process strands at step `h`
# (`S` = survival), so summed over the ensemble it reproduces the stranding
# field as a smooth expectation. It also composes multiplicatively with a
# Fucus **lifetime** `L(t)` — survival is `exp(−∫dt/τ)·L(t)`, today's
# `max_float_days` being a step-function `L`.
#
# **This notebook does no I/O beyond the sidecar and the key.** The
# per-(particle, hour) ingredients of the rate — onshore Stokes `w_on`,
# distance to the BSH H0 coast `dist`, the `024a` hex id, the tidal-flat
# shore flag, and the crow-flies displacement from release — are cached by
# `024d`, so a rate-model member costs seconds per zarr instead of the ~90 s
# the raster + WAM sampling used to cost.
#
# **Rate model: a two-state hazard.**
#
# ```
# 1/τ = 1/τ_calm + trap · r(w_on) / τ_strong     in band, else 0
# r(w) = clip((w − w_c + δ/2) / δ, 0, 1)         step at w_c, optional
#                                                linear edge of width δ
# ```
#
# Nothing in the data constrains the shape of `r` between calm and strong
# wave, so a shape parameter is pure nuisance; the step is the limit of any
# such form and its one knob, `w_c`, means what it says. Once
# `τ_strong ≪` strong-wave-episode duration the model degenerates gracefully
# into an **exposure model**: the beached fraction is the share of particles
# ever in band during a strong-wave hour inside the viability window.
#
# **Land-seeded particles** never enter: the sidecar holds water-seeded
# ("real") drifters only.

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

# One (regime, release_year) per run. surface_stokes is the baseline — the
# beaching driver is the wave field those runs actually felt; surface/bottom
# are sensitivity variants (Stokes was not in their drift).
regime = "surface_stokes"
release_year = 2019
# Restrict to releases in this calendar month (1..12); 0 = whole year.
release_month = 8

# Hex radius (must match an existing key file built by 024a).
hex_radius = 6000

# --- step-hazard parameters ---
# Onshore-Stokes threshold (m/s); THE sweep axis. Measured in-band forcing:
# p75 0.057, p90 0.088, p99 0.15 m/s.
w_c = 0.1
# e-folding time (h) while the onshore forcing is above w_c. Below ~6 h
# everything in band during a Baltic strong-wave episode strands and the
# value stops mattering.
tau_strong_hours = 3.0
# Background in-band timescale (d), independent of forcing; 0 turns it off.
tau_calm_days = 0.0
# Linear edge width (m/s) around w_c; 0 gives a hard step.
delta = 0.0

# Near-shore band width (m); must be ≤ the sidecar's `band_max_m`.
band_m = 2000.0

# Shore-type retention weights, DELIBERATELY DEGENERATE (both 1.0) — the trap
# term is wired but currently expresses nothing, so every shore beaches alike
# for a given wave forcing. The only shore typing available is the BSH
# `H0 <= 0` tidal-flat flag, which is not a retentiveness proxy for *Baltic*
# shores: the basin is effectively tide-free, so the flag fires only in the
# German Bight. The plumbing stays so a real substrate/exposure
# classification can drive it later — set these apart to enable it.
trap_flat = 1.0
trap_wall = 1.0

# Viability / float window: cap each trajectory's contributing age (days)
# before scoring beaching (Rothäusler et al. 2019: weeks to a few months).
# A step-function lifetime; a smooth L(t) would multiply the survival. Must
# be ≤ the sidecar's `window_days`. The scored window is the half-open
# `[0, max_float_days)` in hours — the same convention as 024f's
# `occupancy_max_days`, so the last age bin is never a one-hour sliver.
max_float_days = 60
# Age-bin granularity for the deposition age (days); matches the counts store.
age_bin_days = 10
# Crow-flies travel-distance bin width (km) for the `disp_bin` axis.
disp_bin_km = 10
# Sidecar cadence (hours).
output_dt_hours = 1

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
beaching_path = (
    store_root
    / f"HexAgg_beaching_r{hex_radius}m_{regime}_{release_year}{month_suffix}_{member}.parquet"
)
print(f"beaching → {beaching_path}")

sidecar_dir = output_root / f"BeachingForcing/{regime}/{release_year}"

# %% [markdown]
# # Sidecar → beaching deposition
#
# Per real drifter, deposit the fractional weight that strands at each
# near-shore step (`dep = S_before − S_after`, a telescoping survival
# difference), binned by the coastal hex, elapsed-age bin, shore type and
# travel-distance bin of that step. The surviving weight remaining at the
# window's end is the never-beached residual, recorded once per source hex
# (`beach_hex = -1`, `beach_age_bin = -1`, `disp_bin = -1`,
# `shore_type = "none"`), so deposits + residual per `release_hex` sum to the
# released real drifters and downstream can form the beached fraction against
# the full release pool.
#
# `disp_bin = floor(disp_km / disp_bin_km)` is a **diagnostic axis, never a
# mask**: strong-wave stranding of freshly released material at home is a
# real outcome. The sidecar saturates `disp` at 255 km, so the top bin
# (`255 // disp_bin_km`) pools everything from `disp_bin_km · (255 //
# disp_bin_km)` km outwards.
#
# The rate block below is the same ~15 lines as in
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
# Exposure statistics, one entry per sidecar; the deposited weight per hour is
# accumulated across sidecars so the stranding-age median is hour-resolved.
EXPOSURE = []
DEP_BY_HOUR = np.zeros(max_float_days * 24, dtype=np.float64)
QS = [10, 25, 50, 75, 90, 95, 99, 99.9, 100]
TOP_DISP_BIN = 255 // disp_bin_km


def deposit_one_sidecar(path, release_doy):
    sc = xr.open_zarr(path)
    assert sc.attrs["band_max_m"] >= band_m, (sc.attrs["band_max_m"], band_m)
    assert sc.attrs["window_days"] >= max_float_days, (
        sc.attrs["window_days"], max_float_days
    )
    assert sc.attrs["hex_radius"] == hex_radius, (sc.attrs["hex_radius"], hex_radius)
    assert sc.sizes["obs"] >= max_float_days * 24, (
        sc.sizes["obs"], max_float_days * 24
    )
    # Half-open window [0, max_float_days) in hours, as in 024f: an inclusive
    # end hour would open a one-hour age bin `max_float_days // age_bin_days`.
    sc = sc.isel(obs=slice(0, max_float_days * 24))
    w_on = sc.w_on.values
    dist = sc.dist.values
    flat = sc.flat.values
    hex_at = sc.hex.values
    disp = sc.disp.values
    release_hex = sc.release_hex.values.astype(np.int64)
    ntraj, nobs = w_on.shape

    in_band, r, a = beaching_exponent(w_on, dist, flat)
    A = np.cumsum(a, axis=1)
    dep = np.exp(-(A - a)) - np.exp(-A)
    residual = np.exp(-A[:, -1])

    # Exposure statistics: whether a member is physically sensible is read off
    # the realised forcing and the stranding-age distribution, not off the
    # beached total.
    active = in_band & (r > 0)
    row = {
        "n_traj": ntraj,
        "in_band_steps": int(in_band.sum()),
        "active_steps": int(active.sum()),
        "n_ever_active": int(active.any(axis=1).sum()),
        "beached": float(dep.sum()),
        "total": float(dep.sum() + residual.sum()),
    }
    if in_band.any():
        for q, wq in zip(QS, np.percentile(w_on[in_band].astype("float32"), QS)):
            row[f"w_on_p{q:g}"] = float(wq)
    EXPOSURE.append(row)
    DEP_BY_HOUR[:nobs] += dep.sum(axis=0)

    # Deposition rows, aggregated in age-bin chunks (bounded memory).
    bin_hours = age_bin_days * 24
    frames = []
    for b in range((nobs + bin_hours - 1) // bin_hours):
        sl = slice(b * bin_hours, (b + 1) * bin_hours)
        d = dep[:, sl]
        m = d > 1e-12
        if not m.any():
            continue
        rel = np.broadcast_to(release_hex[:, None], d.shape)[m]
        g = (
            pd.DataFrame({
                "release_hex": rel,
                "beach_hex": hex_at[:, sl][m],
                "shore_type": np.where(flat[:, sl][m], "flat", "wall"),
                "disp_bin": (disp[:, sl][m] // disp_bin_km).astype(np.int64),
                "weight": d[m],
            })
            .groupby(["release_hex", "beach_hex", "shore_type", "disp_bin"],
                     as_index=False)["weight"].sum()
        )
        g["beach_age_bin"] = b
        frames.append(g)
    deposits = (
        pd.concat(frames, ignore_index=True) if frames
        else pd.DataFrame(columns=["release_hex", "beach_hex", "shore_type",
                                   "disp_bin", "weight", "beach_age_bin"])
    )

    residual_rows = (
        pd.DataFrame({"release_hex": release_hex, "weight": residual})
        .loc[lambda df: df["weight"] > 0]
        .groupby("release_hex", as_index=False)["weight"].sum()
        .assign(beach_hex=-1, beach_age_bin=-1, disp_bin=-1, shore_type="none")
    )
    cols = ["release_hex", "release_doy", "beach_hex", "beach_age_bin",
            "shore_type", "disp_bin", "weight"]
    out = pd.concat([deposits, residual_rows], ignore_index=True)
    out["release_doy"] = release_doy
    out = out[cols].astype(
        {"release_hex": "int64", "beach_hex": "int64", "beach_age_bin": "int64",
         "disp_bin": "int64"}
    )
    # Real-drifter count per release_hex, for the conservation check below —
    # collected while the sidecar is open rather than re-derived from `out`,
    # so it counts drifters (weight 1 each), not accumulated fractional weight.
    counts = pd.Series(release_hex).value_counts()
    return out, counts


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
release_counts = pd.Series(dtype=np.int64)
for p, ts, _ in parsed:
    tz = time.time()
    frame, counts = deposit_one_sidecar(p, int(ts.dayofyear))
    frames.append(frame)
    release_counts = release_counts.add(counts, fill_value=0)
    beached = float(frame.loc[frame["beach_hex"] >= 0, "weight"].sum())
    total = float(frame["weight"].sum())
    print(f"  {p.name}: {total:,.0f} drifters, "
          f"{beached:,.0f} beached ({100 * beached / max(total, 1):.1f}%) "
          f"[{time.time() - tz:.1f}s]")

beaching = (
    pd.concat(frames, ignore_index=True)
    .groupby(["release_hex", "release_doy", "beach_hex", "beach_age_bin",
              "shore_type", "disp_bin"], as_index=False)["weight"].sum()
)
print(f"computed {len(beaching):,} rows in {time.time() - t0:.1f}s")

# %% [markdown]
# # Exposure statistics
#
# The plan's decision criteria for a member: how much of the near-shore
# residence is "strong-wave" residence, how many particles ever see a
# strong-wave hour in band, the forcing distribution the threshold is
# cutting, and where the stranding-age mass sits.

# %%
ex = pd.DataFrame(EXPOSURE)
strong_share = ex["active_steps"].sum() / max(ex["in_band_steps"].sum(), 1)
ever_share = ex["n_ever_active"].sum() / max(ex["n_traj"].sum(), 1)
beached_frac = ex["beached"].sum() / max(ex["total"].sum(), 1)
cum = np.cumsum(DEP_BY_HOUR)
median_age_h = (
    int(np.searchsorted(cum, 0.5 * cum[-1])) if cum[-1] > 0 else -1
)
print(f"exposure statistics [member={member}, band_m={band_m:g}, "
      f"max_float_days={max_float_days}]")
print(f"  {'in-band strong-wave hours (r > 0)':42s} {100 * strong_share:8.2f} %")
print(f"  {'particles ever in band, strong wave':42s} {100 * ever_share:8.2f} %")
print(f"  {'beached fraction':42s} {100 * beached_frac:8.2f} %")
age_txt = f"{median_age_h / 24:8.2f} d" if median_age_h >= 0 else f"{'n/a':>8s}"
print(f"  {'median stranding age':42s} {age_txt}")
print(f"  w_onshore over in-band hours (mean of per-sidecar quantiles, m/s):")
for q in QS:
    col = f"w_on_p{q:g}"
    if col in ex:
        print(f"    {'p' + f'{q:g}':>8s} {ex[col].mean():10.4f}")

# %% [markdown]
# # Validation
#
# Every hex id must be in the key — land-seeded / out-of-key positions carry
# `hex = -1` in the sidecar (and the residual rows use `-1` by construction),
# so a stray id means the store was built against a different key. Checked
# before the write, so a mismatched store never lands on disk.

# %%
key_ids = set(pd.read_parquet(key_path, columns=["hex_id"])["hex_id"].astype(int))
seen = set(beaching["release_hex"]) | set(beaching["beach_hex"])
unseen = seen - key_ids - {-1}
if unseen:
    raise ValueError(
        f"{len(unseen)} hex_ids not in {key_path.name}: "
        f"{sorted(unseen)[:10]} ..."
    )
print(f"every release_hex/beach_hex is in {key_path.name} (or -1).")

# %%
# A deposit at an unlabelled hex would masquerade as residual (beach_hex=-1
# is reserved for the never-beached remainder, recorded once per release_hex
# with beach_age_bin=-1 by construction).
bad_residual = beaching[
    (beaching["beach_hex"] == -1) & (beaching["beach_age_bin"] >= 0)
]
if len(bad_residual):
    raise ValueError(
        f"{len(bad_residual)} rows have beach_hex=-1 with beach_age_bin>=0 "
        "— a real deposit masquerading as residual"
    )
print("no beach_hex=-1 row carries a real beach_age_bin.")

# %%
# Conservation: Σweight per release_hex (deposits + residual) must equal the
# real-drifter count released from that hex, across the whole partition.
release_totals = beaching.groupby("release_hex")["weight"].sum()
idx = release_totals.index.union(release_counts.index)
release_totals = release_totals.reindex(idx, fill_value=0.0)
counts = release_counts.reindex(idx, fill_value=0.0)
rel_err = ((release_totals - counts).abs() / counts.replace(0, np.nan)).fillna(0.0)
bad = rel_err[rel_err > 1e-6]
if len(bad):
    raise ValueError(
        f"conservation violated for {len(bad)} release_hex "
        "(Σweight vs. real-drifter count differs by > 1e-6 relative): "
        f"{bad.head(10).to_dict()}"
    )
print(f"Σweight per release_hex conserves the released drifter count "
      f"for all {len(idx)} release_hex (max rel. err "
      f"{rel_err.max() if len(rel_err) else 0:.2e}).")

# %%
beaching.to_parquet(beaching_path)
print(f"wrote {beaching_path} ({beaching_path.stat().st_size / 1e6:.2f} MB)")

# %%
total = float(beaching["weight"].sum())
beached = float(beaching.loc[beaching["beach_hex"] >= 0, "weight"].sum())
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m, member={member}")
print(f"  params: max_float_days={max_float_days}, band_m={band_m:g}, "
      f"w_c={w_c:g}, tau_strong_hours={tau_strong_hours:g}, "
      f"tau_calm_days={tau_calm_days:g}, delta={delta:g}, "
      f"trap_flat/wall={trap_flat:g}/{trap_wall:g}"
      + (" (degenerate — shore type inert)" if trap_flat == trap_wall else ""))
print(f"  drifters (Σweight): {total:,.0f}")
print(f"  beached:           {beached:,.0f} ({100 * beached / max(total, 1):.1f}%)")
print(f"  release_doys:      {beaching['release_doy'].nunique()} "
      f"({beaching['release_doy'].min()}..{beaching['release_doy'].max()})")
print(f"  beach hexes:       {beaching.loc[beaching['beach_hex'] >= 0, 'beach_hex'].nunique():,}")
beach_bins = beaching.loc[beaching["beach_age_bin"] >= 0]
if len(beach_bins):
    print(f"  beach age bins:    {beach_bins['beach_age_bin'].min()}.."
          f"{beach_bins['beach_age_bin'].max()} (× {age_bin_days} d)")
    print(f"  disp bins:         {beach_bins['disp_bin'].min()}.."
          f"{beach_bins['disp_bin'].max()} (× {disp_bin_km} km, "
          f"top bin {TOP_DISP_BIN} saturated)")
