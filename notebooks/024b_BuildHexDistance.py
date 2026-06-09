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
# # Build hex distance histogram partition
#
# Aggregate one `(regime, release_year)` worth of trajectory zarrs into a
# flat `(release_hex, release_doy, distance_bin) → n_traj` table: the
# distribution of per-trajectory crow-flies **final displacement** binned
# per source hex.
#
# Structural sibling of `024_BuildHexAggregates.md` — 024 aggregates
# per-obs occupancy; this aggregates a per-trajectory scalar. The key file
# from `024a_BuildHexKey.md` is a hard prerequisite (its sidecar carries
# the `HexProj` parameters so release-hex labels match the key one-for-one).
#
# Distance bins are additive across years/doys, so downstream consumers
# (027) pool Aug/Sep across years by summing partitions, then derive
# quantiles from the pooled histogram.

# %%
import json
import os
import re
import time
from pathlib import Path

import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr
from dask.distributed import Client

from hextraj import HexProj

# %%
# Pattern: Fucus_BSH_YYYYMMDDTHHMMSS_{regime}_dt{N}min_seed{S}.
# `surface_stokes` must precede `surface` so the alternation matches the
# longer form first.
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


# %%
def release_hex_id(lon0, lat0, hp):
    # Label obs=0 release positions onto the hex lattice, lazy per
    # (trajectory,) chunk. NaN positions (land-seeded) map to the -1
    # sentinel; only valid points are passed to hp.label. Same NaN-safe
    # idiom as 023's release_subbasin.
    def _label(lon, lat):
        out = np.full(lon.shape, -1, dtype=np.int64)
        valid = ~(np.isnan(lon) | np.isnan(lat))
        if valid.any():
            out[valid] = hp.label(lon[valid], lat[valid])
        return out

    return xr.apply_ufunc(
        _label, lon0, lat0,
        dask="parallelized", output_dtypes=[np.int64],
    )

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of trajectory zarrs and write root for the distance partition.
output_root = "../output"

# One (regime, release_year) per run.
regime = "surface"
release_year = 2019

# Hex radius (must match an existing key file built by 024a).
hex_radius = 6000

# Distance histogram bin width (km). The set of bins emerges from the data
# — no upper cap; downstream derives quantiles from the histogram.
distance_bin_km = 1.0

# %% [markdown]
# # Derived layout / key + projection

# %%
output_root = Path(output_root)
store_root = output_root / "HexAggregates"
store_root.mkdir(parents=True, exist_ok=True)

key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
meta_path = key_path.with_suffix(".json")
distance_path = (
    store_root / f"HexAgg_distance_r{hex_radius}m_{regime}_{release_year}.parquet"
)

if not key_path.exists() or not meta_path.exists():
    raise FileNotFoundError(
        f"Key file or sidecar missing — run 024a_BuildHexKey.md first.\n"
        f"  expected: {key_path}\n  expected: {meta_path}"
    )

meta = json.loads(meta_path.read_text())
hp = HexProj(**meta["hex_proj"])
print(f"HexProj: {meta['hex_proj']}")
print(f"distance → {distance_path}")

# %% [markdown]
# # Dask cluster

# %%
scheduler_file = os.environ.get("SCHEDULER_FILE")
if scheduler_file:
    for _ in range(60):
        if os.path.exists(scheduler_file):
            break
        time.sleep(1)
    client = Client(scheduler_file=scheduler_file)
else:
    client = Client(ip="0.0.0.0")
client

# %% [markdown]
# # Trajectory zarrs → distance histogram

# %%
zarrs = sorted(
    (output_root / f"Trajectories/{regime}/{release_year}").glob("*.zarr")
)
parsed = [(p, *parse_zarr_stem(p)) for p in zarrs]
for p, ts, fn_regime in parsed:
    assert ts.year == release_year, (ts, release_year, p)
    assert fn_regime == regime, (fn_regime, regime, p)

if not parsed:
    raise FileNotFoundError(
        f"no zarrs at {output_root}/Trajectories/{regime}/{release_year}/"
    )
release_doys = sorted({int(ts.dayofyear) for _, ts, _ in parsed})
print(f"{len(parsed)} zarrs, release_doys "
      f"{release_doys[0]}..{release_doys[-1]} ({len(release_doys)} unique)")


# %%
def zarr_to_distance_frame(path, release_doy):
    """Lazy (release_hex, release_doy, distance_bin) frame for one zarr.
    One row per trajectory: its source hex and the bin of its crow-flies
    final displacement (equirectangular 111 km/deg metric, see
    docs/distance_calculation.md).

    Two NaN paths, both handled: land-seeded particles (zero first-step
    displacement) are masked to NaN throughout, so their release_lon → -1
    hex and NaN distance, and they are dropped downstream; domain-exited
    particles (NaN tail mid-run) carry their last in-domain distance
    forward via ffill, which is the intended last-valid-position value."""
    ds = xr.open_zarr(path)
    on_land = (
        (ds.lon.diff("obs").isel(obs=0, drop=True) == 0)
        & (ds.lat.diff("obs").isel(obs=0, drop=True) == 0)
    )
    ds = ds[["lon", "lat"]].where(~on_land)

    # NaN for land-seeded trajectories (masked above); release_hex_id maps
    # NaN → -1, so those rows drop out at the release_hex >= 0 filter.
    release_lon = ds.lon.isel(obs=0, drop=True)
    release_lat = ds.lat.isel(obs=0, drop=True)
    dlat = ds.lat - release_lat
    dlon = (ds.lon - release_lon) * np.cos(np.deg2rad(release_lat))
    dist = 111.0 * np.sqrt(dlat ** 2 + dlon ** 2)
    # Final displacement = crow-flies distance at the last valid obs.
    final_dist = dist.ffill("obs").isel(obs=-1, drop=True)

    frame = xr.Dataset(
        dict(
            release_hex=release_hex_id(release_lon, release_lat, hp),
            distance_bin=(final_dist // distance_bin_km),
        )
    )
    ddf = frame.to_dask_dataframe()
    ddf = ddf[ddf["release_hex"] >= 0].dropna(subset=["distance_bin"])
    ddf = ddf.astype({"release_hex": "int64", "distance_bin": "int64"})
    ddf["release_doy"] = release_doy
    return ddf[["release_hex", "release_doy", "distance_bin"]]


t0 = time.time()
counts = (
    dd.concat([zarr_to_distance_frame(p, ts.dayofyear) for p, ts, _ in parsed])
    .groupby(["release_hex", "release_doy", "distance_bin"])
    .size().rename("n_traj").reset_index()
    .compute()
    .reset_index(drop=True)  # collapse the per-partition RangeIndex stack
)
print(f"computed {len(counts):,} rows in {time.time() - t0:.1f}s")

# %%
counts.to_parquet(distance_path)
print(f"wrote {distance_path} ({distance_path.stat().st_size / 1e6:.2f} MB)")

# %% [markdown]
# # Validation

# %%
key_ids = set(pd.read_parquet(key_path, columns=["hex_id"])["hex_id"].astype(int))
unseen = set(counts["release_hex"].astype(int)) - key_ids
if unseen:
    print(f"WARNING: {len(unseen)} release_hex not in {key_path.name}: "
          f"{sorted(unseen)[:10]} ...")
else:
    print(f"every release_hex is in {key_path.name}.")

# %%
print(f"regime={regime}, release_year={release_year}, hex_radius={hex_radius} m")
print(f"  distance_bin_km: {distance_bin_km}")
print(f"  rows:            {len(counts):,}")
print(f"  trajectories:    {int(counts['n_traj'].sum()):,}")
print(f"  release_doys:    {counts['release_doy'].nunique()} "
      f"({counts['release_doy'].min()}..{counts['release_doy'].max()})")
print(f"  source hexes:    {counts['release_hex'].nunique():,}")
print(f"  distance_bins:   {counts['distance_bin'].min()}.."
      f"{counts['distance_bin'].max()} "
      f"({counts['distance_bin'].max() * distance_bin_km:g} km max)")
