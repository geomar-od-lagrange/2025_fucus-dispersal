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
# — the same near-shore beaching rate `τ = τ0/(trap·g(w_onshore))` as
# `024d_BuildBeaching` (see [beaching.md](../docs/beaching.md)). `A` only
# grows inside the near-shore band, so open-water residence is undiluted; a
# particle that lingers near a retentive shore loses weight fast. This is the
# deterministic occupancy analogue of 024d's fractional stranding: 024d
# records where the weight *leaves* (`beach_hex`), this records where the
# still-drifting weight *is* (`target_hex`).
#
# Emits **two** weights per bin so the fold-in is a self-consistent
# comparison: `occ` (plain, every sample = 1) and `surv` (Σ S) over the same
# window and hexing. `030_SurvivalHeatmaps` maps `surv`, `occ`, and the
# surviving fraction `surv/occ`.
#
# Reuses 024d's raster + Stokes machinery; aggregates with `np.bincount`
# (occupancy is dense — every obs, every hex — unlike 024d's sparse in-band
# deposits). Partitioned per `(regime, year, month)`; `030` pools. The
# `024a` key is a hard prerequisite (its `HexProj` + hex-id set).

# %%
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Transformer
from scipy import ndimage

from hextraj import HexProj

# %%
# Pattern: Fucus_BSH_YYYYMMDDTHHMMSS_{regime}_dt{N}min_seed{S}.
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
data_root = "../data"
output_root = "../output"

regime = "surface_stokes"
release_year = 2019
# Restrict to releases in this calendar month (1..12); 0 = whole year.
release_month = 8

hex_radius = 6000

# Occupancy horizon (days): how far along each trajectory to accumulate
# residence. Must be ≥ the largest map horizon in 030. Decoupled from 024d's
# `max_float_days` — survival occupancy runs to the mapping horizon, not the
# viability cutoff (a Fucus lifetime L(t) would enter as an extra factor on S).
occupancy_max_days = 120
# Age-bin granularity (days); matches the counts store.
age_bin_days = 10
# Zarr output cadence (hours).
output_dt_hours = 1

# Rate-model parameters — must match 024d for a consistent survival field.
band_m = 2000.0
tau0_hours = 24.0
trap_flat = 2.0
trap_wall = 1.0
w_half = 0.05

# Distance-to-coast raster resolution (m, EPSG:3035).
raster_dx_m = 500.0

# %% [markdown]
# # Derived layout / key + projection

# %%
data_root = Path(data_root)
output_root = Path(output_root)
store_root = output_root / "HexAggregates"
store_root.mkdir(parents=True, exist_ok=True)

key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
meta_path = key_path.with_suffix(".json")
if not key_path.exists() or not meta_path.exists():
    raise FileNotFoundError(
        f"Key file or sidecar missing — run 024a_BuildHexKey.md first.\n"
        f"  expected: {key_path}\n  expected: {meta_path}"
    )

month_suffix = f"_m{release_month:02d}" if release_month else ""
survocc_path = (
    store_root
    / f"HexAgg_survocc_r{hex_radius}m_{regime}_{release_year}{month_suffix}.parquet"
)

meta = json.loads(meta_path.read_text())
hp = HexProj(**meta["hex_proj"])
print(f"HexProj: {meta['hex_proj']}")
print(f"survocc → {survocc_path}")

# Contiguous hex index for bincount: hex_id → 0..n-1 and back.
hex_ids = pd.read_parquet(key_path, columns=["hex_id"])["hex_id"].astype(int).to_numpy()
hexid_to_idx = pd.Series(np.arange(len(hex_ids)), index=hex_ids)
n_hex = len(hex_ids)
n_agebin = occupancy_max_days // age_bin_days
print(f"{n_hex:,} hexes, {n_agebin} age bins × {age_bin_days} d")

stokes_dir = output_root / "stokes" / "baltic_highres" / str(release_year)

# %% [markdown]
# # Beaching geometry raster (+ dense hex index per cell)
#
# Same field as 024d — distance to the BSH H0 coast, seaward normal rotated
# to geographic east/north, wall/flat shore type — plus, here, the dense hex
# **index** (into the key's hex-id order) of each raster cell, so per-position
# hexing is a lookup ready for `np.bincount`.

# %%
def _h0_nearest_sampler(h0):
    lon = h0.lon.values
    lat = h0.lat.values
    values = h0.H0.values
    lon0, dlon = lon[0], lon[1] - lon[0]
    lat0, dlat = lat[0], lat[1] - lat[0]

    def sample(lon_q, lat_q):
        i = np.round((lon_q - lon0) / dlon).astype(np.int64)
        j = np.round((lat_q - lat0) / dlat).astype(np.int64)
        ok = (i >= 0) & (i < len(lon)) & (j >= 0) & (j < len(lat))
        out = np.full(lon_q.shape, np.nan)
        out[ok] = values[j[ok], i[ok]]
        return out

    return sample


def build_survival_raster(data_root, dx_m, hp, hexid_to_idx):
    h0_fine = xr.open_dataset(
        data_root / "bsh_hbmnoku_static/static_file_fine/H0_file_fine.nc"
    )
    h0_coarse = xr.open_dataset(
        data_root / "bsh_hbmnoku_static/static_file_coarse/H0_file_coarse.nc"
    )
    sample_fine = _h0_nearest_sampler(h0_fine)
    sample_coarse = _h0_nearest_sampler(h0_coarse)
    fine_bbox = (
        float(h0_fine.lon.min()), float(h0_fine.lon.max()),
        float(h0_fine.lat.min()), float(h0_fine.lat.max()),
    )

    to_3035 = Transformer.from_crs(4326, 3035, always_xy=True)
    from_3035 = Transformer.from_crs(3035, 4326, always_xy=True)

    lon2d, lat2d = np.meshgrid(h0_coarse.lon.values, h0_coarse.lat.values)
    x, y = to_3035.transform(lon2d.ravel(), lat2d.ravel())
    pad = 5000.0
    xmin, xmax = x.min() - pad, x.max() + pad
    ymin, ymax = y.min() - pad, y.max() + pad
    xs = np.arange(xmin, xmax + dx_m, dx_m)
    ys = np.arange(ymax, ymin - dx_m, -dx_m)
    ncol, nrow = len(xs), len(ys)

    gx, gy = np.meshgrid(xs, ys)
    lon_g, lat_g = from_3035.transform(gx.ravel(), gy.ravel())
    in_fine = (
        (lon_g >= fine_bbox[0]) & (lon_g <= fine_bbox[1])
        & (lat_g >= fine_bbox[2]) & (lat_g <= fine_bbox[3])
    )
    h0r = sample_coarse(lon_g, lat_g)
    h0r[in_fine] = sample_fine(lon_g[in_fine], lat_g[in_fine])
    h0r = h0r.reshape(nrow, ncol)

    water = np.isfinite(h0r)
    flat = water & (h0r <= 0)

    dist_cells, (jy, jx) = ndimage.distance_transform_edt(water, return_indices=True)
    dist_m = (dist_cells * dx_m).astype("float32")
    flat_fronted_land = (~water) & ndimage.binary_dilation(flat)
    nearest_flat = flat_fronted_land[jy, jx]

    grad_row, grad_col = np.gradient(dist_m.astype("float64"), dx_m)
    proj_e = grad_col
    proj_n = -grad_row
    d = 1e-3
    xe1, ye1 = to_3035.transform(lon_g + d, lat_g)
    xe0, ye0 = to_3035.transform(lon_g - d, lat_g)
    e_x, e_y = xe1 - xe0, ye1 - ye0
    e_norm = np.hypot(e_x, e_y)
    xn1, yn1 = to_3035.transform(lon_g, lat_g + d)
    xn0, yn0 = to_3035.transform(lon_g, lat_g - d)
    n_x, n_y = xn1 - xn0, yn1 - yn0
    n_norm = np.hypot(n_x, n_y)
    with np.errstate(invalid="ignore", divide="ignore"):
        e_x, e_y = (e_x / e_norm).reshape(nrow, ncol), (e_y / e_norm).reshape(nrow, ncol)
        n_x, n_y = (n_x / n_norm).reshape(nrow, ncol), (n_y / n_norm).reshape(nrow, ncol)
    n_east = proj_e * e_x + proj_n * e_y
    n_north = proj_e * n_x + proj_n * n_y
    mag = np.hypot(n_east, n_north)
    with np.errstate(invalid="ignore", divide="ignore"):
        n_out_x = np.where(mag > 0, n_east / mag, 0.0).astype("float32")
        n_out_y = np.where(mag > 0, n_north / mag, 0.0).astype("float32")

    # Dense hex index of every cell centre (vectorised map; -1 outside key).
    lon_good = np.isfinite(lon_g) & np.isfinite(lat_g)
    hex_id_flat = np.full(lon_g.shape, -1, dtype=np.int64)
    hex_id_flat[lon_good] = hp.label(lon_g[lon_good], lat_g[lon_good])
    hex_idx = (
        hexid_to_idx.reindex(hex_id_flat).to_numpy()
    )
    hex_idx = np.where(np.isnan(hex_idx), -1, hex_idx).astype(np.int64).reshape(nrow, ncol)

    def to_rowcol(lon_q, lat_q):
        x, y = to_3035.transform(lon_q, lat_q)
        col = np.round((x - xmin) / dx_m)
        row = np.round((ymax - y) / dx_m)
        ok = np.isfinite(x) & (col >= 0) & (col < ncol) & (row >= 0) & (row < nrow)
        row = np.clip(np.nan_to_num(row), 0, nrow - 1).astype(np.int64)
        col = np.clip(np.nan_to_num(col), 0, ncol - 1).astype(np.int64)
        return row, col, ok

    return dict(
        dist_m=dist_m, n_out_x=n_out_x, n_out_y=n_out_y, nearest_flat=nearest_flat,
        hex_idx=hex_idx, to_rowcol=to_rowcol, water_cells=int(water.sum()),
        shape=(nrow, ncol),
    )


t0 = time.time()
rast = build_survival_raster(data_root, raster_dx_m, hp, hexid_to_idx)
print(f"raster {rast['shape']} water={rast['water_cells']:,} in {time.time() - t0:.1f}s")

# %% [markdown]
# # Onshore Stokes sampler (as 024d)

# %%
class OnshoreStokes:
    """Nearest-neighbour raw Stokes sampler with a one-day file cache."""

    def __init__(self, stokes_dir):
        self.stokes_dir = Path(stokes_dir)
        self._day = None
        self._cube = None
        first = sorted(self.stokes_dir.glob("stokes_*.nc"))
        if not first:
            raise FileNotFoundError(f"no raw Stokes under {self.stokes_dir}")
        g = xr.open_dataset(first[0])
        self.lon = g.longitude.values
        self.lat = g.latitude.values
        self.missing_days = 0

    def _load(self, day):
        if day == self._day:
            return self._cube
        path = self.stokes_dir / f"stokes_{day.strftime('%Y%m%d')}.nc"
        if not path.exists():
            self.missing_days += 1
            self._day, self._cube = day, None
            return None
        g = xr.open_dataset(path)
        self._cube = (g.VSDX.values, g.VSDY.values, g.time.values)
        self._day = day
        return self._cube

    def onshore(self, lon_q, lat_q, when, n_out_x, n_out_y):
        cube = self._load(pd.Timestamp(when).normalize())
        if cube is None:
            return np.zeros(lon_q.shape, dtype="float32")
        vsdx, vsdy, times = cube
        t = np.argmin(np.abs(times - np.datetime64(when)))
        i = np.round((lon_q - self.lon[0]) / (self.lon[1] - self.lon[0]))
        j = np.round((lat_q - self.lat[0]) / (self.lat[1] - self.lat[0]))
        ok = (
            np.isfinite(lon_q)
            & (i >= 0) & (i < len(self.lon)) & (j >= 0) & (j < len(self.lat))
        )
        i = np.clip(np.nan_to_num(i), 0, len(self.lon) - 1).astype(np.int64)
        j = np.clip(np.nan_to_num(j), 0, len(self.lat) - 1).astype(np.int64)
        sx = vsdx[t][j, i]
        sy = vsdy[t][j, i]
        onsh = -(sx * n_out_x + sy * n_out_y)
        onsh = np.where(ok & np.isfinite(onsh), onsh, 0.0)
        return np.maximum(0.0, onsh).astype("float32")


# %% [markdown]
# # Trajectory zarrs → survival-weighted occupancy
#
# Per real drifter, weight every in-window occupancy sample by `S = exp(−A)`;
# `A` accrues only inside the near-shore band. Aggregate `occ` (=1) and `surv`
# (=S) over `(target_hex, age_bin)` with `np.bincount`.

# %%
def occupancy_one_zarr(path, release_doy, rast, stokes):
    ds = xr.open_zarr(path).isel(obs=slice(0, occupancy_max_days * 24))
    lon = ds.lon.values.astype("float32")
    lat = ds.lat.values.astype("float32")
    ntraj, nobs = lon.shape

    real = (
        np.isfinite(lon[:, 0])
        & ~((lon[:, 1] - lon[:, 0] == 0) & (lat[:, 1] - lat[:, 0] == 0))
    )

    row, col, ok = rast["to_rowcol"](lon.ravel(), lat.ravel())
    dist_at = rast["dist_m"][row, col].reshape(ntraj, nobs)
    nx_at = rast["n_out_x"][row, col].reshape(ntraj, nobs)
    ny_at = rast["n_out_y"][row, col].reshape(ntraj, nobs)
    flat_at = rast["nearest_flat"][row, col].reshape(ntraj, nobs)
    hex_idx = rast["hex_idx"][row, col].reshape(ntraj, nobs)
    ok = ok.reshape(ntraj, nobs)
    in_band = (dist_at < band_m) & ok

    release_time = np.datetime64(ds.time.isel(trajectory=0, obs=0).values)
    abs_time = release_time + np.arange(nobs).astype("timedelta64[h]")
    w_on = np.zeros((ntraj, nobs), dtype="float32")
    for h in range(nobs):
        m = in_band[:, h]
        if not m.any():
            continue
        w_on[m, h] = stokes.onshore(
            lon[m, h], lat[m, h], abs_time[h], nx_at[m, h], ny_at[m, h]
        )

    trap = np.where(flat_at, trap_flat, trap_wall).astype("float32")
    g_w = w_on / (w_on + w_half)
    a = np.where(in_band, output_dt_hours / (tau0_hours / (trap * np.maximum(g_w, 1e-6))), 0.0)
    surv = np.exp(-np.cumsum(a, axis=1)).astype("float32")

    age_bin = (np.arange(nobs) // (age_bin_days * 24)).astype(np.int64)
    age_bin2d = np.broadcast_to(age_bin[None, :], (ntraj, nobs))
    valid = ok & (hex_idx >= 0) & real[:, None]

    key = hex_idx[valid] * n_agebin + age_bin2d[valid]
    length = n_hex * n_agebin
    occ = np.bincount(key, minlength=length).astype(np.float64)
    surv_agg = np.bincount(key, weights=surv[valid].astype(np.float64), minlength=length)
    occ = occ.reshape(n_hex, n_agebin)
    surv_agg = surv_agg.reshape(n_hex, n_agebin)

    hi, bi = np.nonzero(occ > 0)
    return pd.DataFrame({
        "release_doy": release_doy,
        "age_bin": bi.astype(np.int64),
        "target_hex": hex_ids[hi],
        "occ": occ[hi, bi],
        "surv": surv_agg[hi, bi],
    })


# %%
zarrs = sorted(
    (output_root / f"Trajectories/{regime}/{release_year}").glob("*.zarr")
)
parsed = [(p, *parse_zarr_stem(p)) for p in zarrs]
for p, ts, fn_regime in parsed:
    assert ts.year == release_year, (ts, release_year, p)
    assert fn_regime == regime, (fn_regime, regime, p)
if release_month:
    parsed = [x for x in parsed if x[1].month == release_month]
if not parsed:
    raise FileNotFoundError(
        f"no zarrs at {output_root}/Trajectories/{regime}/{release_year}/"
        + (f" for month {release_month}" if release_month else "")
    )
release_doys = sorted({int(ts.dayofyear) for _, ts, _ in parsed})
print(f"{len(parsed)} zarrs, release_doys "
      f"{release_doys[0]}..{release_doys[-1]} ({len(release_doys)} unique)")

# %%
stokes = OnshoreStokes(stokes_dir)

t0 = time.time()
frames = []
for p, ts, _ in parsed:
    tz = time.time()
    frame = occupancy_one_zarr(p, int(ts.dayofyear), rast, stokes)
    frames.append(frame)
    o = frame["occ"].sum()
    s = frame["surv"].sum()
    print(f"  {p.name}: occ {o:,.0f}, surv {s:,.0f} "
          f"(drifting {s / max(o, 1):.3f}) [{time.time() - tz:.1f}s]")

survocc = (
    pd.concat(frames, ignore_index=True)
    .groupby(["release_doy", "age_bin", "target_hex"], as_index=False)[["occ", "surv"]].sum()
)
print(f"computed {len(survocc):,} rows in {time.time() - t0:.1f}s; "
      f"missing Stokes days: {stokes.missing_days}")

# %%
survocc.to_parquet(survocc_path)
print(f"wrote {survocc_path} ({survocc_path.stat().st_size / 1e6:.2f} MB)")

# %% [markdown]
# # Validation

# %%
key_ids = set(hex_ids.tolist())
unseen = set(survocc["target_hex"]) - key_ids
if unseen:
    print(f"WARNING: {len(unseen)} target_hex not in key: {sorted(unseen)[:10]} ...")
else:
    print("every target_hex is in the key.")

# %%
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m")
print(f"  params: occupancy_max_days={occupancy_max_days}, band_m={band_m:g}, "
      f"tau0_hours={tau0_hours:g}, trap_flat/wall={trap_flat:g}/{trap_wall:g}")
per_bin = survocc.groupby("age_bin")[["occ", "surv"]].sum()
per_bin["drifting"] = per_bin["surv"] / per_bin["occ"]
print(per_bin.to_string(float_format=lambda v: f"{v:,.3f}"))
