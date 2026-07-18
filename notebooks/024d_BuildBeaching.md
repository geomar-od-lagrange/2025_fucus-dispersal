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

# Build beaching store partition (weighted)

Post-simulation beaching pass over one `(regime, release_year)` worth of
trajectory zarrs (see [beaching.md](../plans/done/beaching.md)). Each real
drifter carries unit surviving (free-drifting) weight; at every near-shore
step a fraction of that weight *strands* at the current coastal hex and
leaves the drifting pool, with the remainder continuing to drift. This is
the **fractional / weighted** scheme — the deterministic expectation of a
per-step first-stranding process, with no Monte-Carlo noise. Writes a flat
`(release_hex, release_doy, beach_hex, beach_age_bin, shore_type) → weight`
table — additive across `release_doy`/year like the other `024x` stores.

**Why weighted, not stochastic.** Releases carry only ~100 particles per
seeding cell, so a per-particle random-strand rule gives noisy coverage at
the high-age tail exactly where the free-drifting/beached split matters.
The weighted deposit `dep[t,h] = S(t,h)·(1−e^{−Δt/τ})` equals the
probability that a first-stranding process strands at step `h`
(`S` = survival), so summed over the ensemble it reproduces the stranding
field as a smooth expectation. It also composes multiplicatively with a
Fucus **lifetime** `L(t)` — survival is `exp(−∫dt/τ)·L(t)`, today's
`max_float_days` being a step-function `L`.

**Three-ingredient rate model.** `τ = τ0 / (trap(shore_type)·g(w_onshore))`:

- **distance to shore** — a rasterised distance-to-coast field built from
  the **BSH H0 land-sea mask** (finite `H0` = water, NaN = land, `H0 ≤ 0`
  = tidal flat), fine-over-coarse, in EPSG:3035. This is the mask the
  particles were advected on; the coastline geojson polygons miss a large
  fraction of genuine water positions and are not used here.
- **shore type** — the nearest land cell fronted by a tidal-flat
  (`H0 ≤ 0`) cell reads as `flat` (dissipative, retentive), else `wall`.
- **onshore wave forcing** — the onshore component of the raw
  `baltic_highres` Stokes drift (`VSDX/VSDY`), i.e. the cross-shore
  transport the `surface_stokes` runs masked at blocked faces, sampled
  here *before* that mask.

**Land-seeded particles are dropped** (zero first-step displacement),
exactly as `024`/`024b` do. The key file from `024a_BuildHexKey.md` is a
hard prerequisite — its sidecar carries the `HexProj` used to label hexes.

```python
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
```

```python
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
```

# Parameters

```python tags=["parameters"]
# Read root of the data twin (BSH static H0) and of the trajectory zarrs +
# raw Stokes; write root for the beaching partition.
data_root = "../data"
output_root = "../output"

# One (regime, release_year) per run. surface_stokes is the baseline — the
# beaching driver is the wave field those runs actually felt; surface/bottom
# are sensitivity variants (Stokes was not in their drift).
regime = "surface_stokes"
release_year = 2019

# Restrict to releases in this calendar month (1..12); 0 = whole year. When
# set, the output filename gets a `_mMM` suffix (a partial-year store).
release_month = 8

# Hex radius (must match an existing key file built by 024a).
hex_radius = 6000

# Viability / float window: cap each trajectory's contributing age (days)
# before scoring beaching (Rothäusler et al. 2019: weeks to a few months).
# A step-function lifetime; a smooth L(t) would multiply the survival.
max_float_days = 60
# Age-bin granularity for the deposition age (days); matches the counts store.
age_bin_days = 10
# Zarr output cadence (hours).
output_dt_hours = 1

# Rate-model parameters (see beaching.md "Open questions"). Sweep + report.
band_m = 2000.0        # near-shore band width (m)
tau0_hours = 24.0      # base e-folding beaching timescale (h)
trap_flat = 2.0        # retention weight, dissipative flat shore
trap_wall = 1.0        # retention weight, reflective wall shore
w_half = 0.05          # onshore-Stokes half-saturation (m/s) in the g ramp

# Distance-to-coast raster resolution (m, EPSG:3035).
raster_dx_m = 500.0
```

# Derived layout / key + projection

```python
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
beaching_path = (
    store_root
    / f"HexAgg_beaching_r{hex_radius}m_{regime}_{release_year}{month_suffix}.parquet"
)

meta = json.loads(meta_path.read_text())
hp = HexProj(**meta["hex_proj"])
print(f"HexProj: {meta['hex_proj']}")
print(f"beaching → {beaching_path}")

stokes_dir = output_root / "stokes" / "baltic_highres" / str(release_year)
```

# Beaching geometry raster

Build the near-shore geometry field once: a regular EPSG:3035 raster
carrying, per cell, the distance to the nearest BSH land cell, the seaward
unit normal `n_out` (∇distance, rotated into geographic east/north), the
`shore_type`, and the `024a` hex id of the cell centre (so per-position hex
labels are a lookup, not a per-point projection). Particle positions are
sampled against it by nearest cell.

```python
def _h0_nearest_sampler(h0):
    """Nearest-cell H0 lookup for a regular BSH grid; NaN outside the grid."""
    lon = h0.lon.values
    lat = h0.lat.values
    values = h0.H0.values
    lon0, dlon = lon[0], lon[1] - lon[0]
    lat0, dlat = lat[0], lat[1] - lat[0]  # dlat < 0 (descending)

    def sample(lon_q, lat_q):
        i = np.round((lon_q - lon0) / dlon).astype(np.int64)
        j = np.round((lat_q - lat0) / dlat).astype(np.int64)
        ok = (i >= 0) & (i < len(lon)) & (j >= 0) & (j < len(lat))
        out = np.full(lon_q.shape, np.nan)
        out[ok] = values[j[ok], i[ok]]
        return out

    return sample


def build_beaching_raster(data_root, dx_m, hp):
    """Distance-to-coast field + seaward normal + shore type + hex id on a
    3035 raster.

    Water = finite H0 (fine grid over coarse), land = NaN, tidal flat =
    finite H0 ≤ 0. Returns a dict with the raster arrays and the affine
    parameters + a lon/lat→(row, col) mapper for sampling.
    """
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

    # Raster extent = coarse-grid footprint in 3035, padded one coarse cell.
    lon2d, lat2d = np.meshgrid(h0_coarse.lon.values, h0_coarse.lat.values)
    x, y = to_3035.transform(lon2d.ravel(), lat2d.ravel())
    pad = 5000.0
    xmin, xmax = x.min() - pad, x.max() + pad
    ymin, ymax = y.min() - pad, y.max() + pad
    xs = np.arange(xmin, xmax + dx_m, dx_m)
    ys = np.arange(ymax, ymin - dx_m, -dx_m)  # descending → row 0 is north
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
    # shore_type: nearest land cell adjacent to any tidal-flat cell → flat.
    flat_fronted_land = (~water) & ndimage.binary_dilation(flat)
    nearest_flat = flat_fronted_land[jy, jx]

    # Seaward normal = ∇distance (distance grows into open water), in the
    # projected (3035) plane: +col = +easting, row increases southward.
    grad_row, grad_col = np.gradient(dist_m.astype("float64"), dx_m)
    proj_e = grad_col               # component along projected easting
    proj_n = -grad_row              # component along projected northing

    # Rotate the normal from the 3035 grid frame into geographic east/north
    # so it dots consistently with the geographic Stokes components
    # (VSDX = east, VSDY = north). LAEA meridian convergence rotates grid
    # north away from true north by up to ~17° in the eastern Baltic.
    d = 1e-3  # degrees (~100 m); step for the local basis
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

    # Hex id of every cell centre (labelled once; per-position hex is then a
    # lookup, not a per-point projection).
    hex_id = np.full(lon_g.shape, -1, dtype=np.int64)
    good = np.isfinite(lon_g) & np.isfinite(lat_g)
    hex_id[good] = hp.label(lon_g[good], lat_g[good])
    hex_id = hex_id.reshape(nrow, ncol)

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
        hex_id=hex_id, to_rowcol=to_rowcol, water_cells=int(water.sum()),
        flat_cells=int(flat.sum()), shape=(nrow, ncol),
    )


t0 = time.time()
rast = build_beaching_raster(data_root, raster_dx_m, hp)
print(f"raster {rast['shape']} water={rast['water_cells']:,} "
      f"flat={rast['flat_cells']:,} in {time.time() - t0:.1f}s")
```

# Onshore Stokes sampler

Sample the raw `baltic_highres` `VSDX/VSDY` by nearest hour and cell,
caching one day-file at a time. Positions outside the Baltic wave grid
(e.g. the German Bight strip < 9 °E) return zero onshore transport.

```python
class OnshoreStokes:
    """Nearest-neighbour raw Stokes sampler with a one-day file cache."""

    def __init__(self, stokes_dir):
        self.stokes_dir = Path(stokes_dir)
        self._day = None
        self._cube = None  # (VSDX, VSDY, times)
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
        """Onshore Stokes magnitude max(0, -(VSDX,VSDY)·n_out) at each point."""
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
```

# Trajectory zarrs → beaching deposition

Per real drifter, deposit the fractional weight that strands at each
near-shore step (`dep = S_before − S_after`, a telescoping survival
difference), binned by the coastal hex, elapsed-age bin, and shore type of
that step. The surviving weight remaining at the window's end is the
never-beached residual, recorded once per source hex
(`beach_hex = -1`, `beach_age_bin = -1`, `shore_type = "none"`), so
downstream can form the beached fraction against the full release pool.

```python
def deposit_one_zarr(path, release_doy, rast, stokes):
    ds = xr.open_zarr(path).isel(obs=slice(0, max_float_days * 24 + 1))
    lon = ds.lon.values.astype("float32")
    lat = ds.lat.values.astype("float32")
    ntraj, nobs = lon.shape

    # Land-seeded = zero first-step displacement (as 024/024b).
    real = (
        np.isfinite(lon[:, 0])
        & ~((lon[:, 1] - lon[:, 0] == 0) & (lat[:, 1] - lat[:, 0] == 0))
    )

    row, col, ok = rast["to_rowcol"](lon.ravel(), lat.ravel())
    dist_at = rast["dist_m"][row, col].reshape(ntraj, nobs)
    nx_at = rast["n_out_x"][row, col].reshape(ntraj, nobs)
    ny_at = rast["n_out_y"][row, col].reshape(ntraj, nobs)
    flat_at = rast["nearest_flat"][row, col].reshape(ntraj, nobs)
    hex_at = rast["hex_id"][row, col].reshape(ntraj, nobs)
    in_band = (dist_at < band_m) & ok.reshape(ntraj, nobs)

    # Onshore Stokes at in-band positions, hour by hour (day-file cached).
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

    # Per-step beaching exponent a = Δt/τ (0 outside the band), then the
    # telescoping survival deposit dep[t,h] = e^{-A_before} - e^{-A_after}.
    trap = np.where(flat_at, trap_flat, trap_wall).astype("float32")
    g_w = w_on / (w_on + w_half)
    a = np.where(in_band, output_dt_hours / (tau0_hours / (trap * np.maximum(g_w, 1e-6))), 0.0)
    A = np.cumsum(a, axis=1)
    dep = np.exp(-(A - a)) - np.exp(-A)
    dep[~real] = 0.0
    residual = np.exp(-A[:, -1]) * real

    release_hex = np.full(ntraj, -1, dtype=np.int64)
    good0 = real & np.isfinite(lon[:, 0]) & np.isfinite(lat[:, 0])
    release_hex[good0] = hp.label(lon[good0, 0], lat[good0, 0])

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
                "shore_type": np.where(flat_at[:, sl][m], "flat", "wall"),
                "weight": d[m],
            })
            .groupby(["release_hex", "beach_hex", "shore_type"], as_index=False)["weight"].sum()
        )
        g["beach_age_bin"] = b
        frames.append(g)
    deposits = (
        pd.concat(frames, ignore_index=True) if frames
        else pd.DataFrame(columns=["release_hex", "beach_hex", "shore_type", "weight", "beach_age_bin"])
    )

    residual_rows = (
        pd.DataFrame({"release_hex": release_hex, "weight": residual})
        .loc[lambda df: df["weight"] > 0]
        .groupby("release_hex", as_index=False)["weight"].sum()
        .assign(beach_hex=-1, beach_age_bin=-1, shore_type="none")
    )
    cols = ["release_hex", "release_doy", "beach_hex", "beach_age_bin", "shore_type", "weight"]
    out = pd.concat([deposits, residual_rows], ignore_index=True)
    out["release_doy"] = release_doy
    return out[cols].astype(
        {"release_hex": "int64", "beach_hex": "int64", "beach_age_bin": "int64"}
    )
```

```python
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
```

```python
stokes = OnshoreStokes(stokes_dir)

t0 = time.time()
frames = []
for p, ts, _ in parsed:
    tz = time.time()
    frame = deposit_one_zarr(p, int(ts.dayofyear), rast, stokes)
    frames.append(frame)
    beached = float(frame.loc[frame["beach_hex"] >= 0, "weight"].sum())
    total = float(frame["weight"].sum())
    print(f"  {p.name}: {total:,.0f} drifters, "
          f"{beached:,.0f} beached ({100 * beached / max(total, 1):.1f}%) "
          f"[{time.time() - tz:.1f}s]")

beaching = (
    pd.concat(frames, ignore_index=True)
    .groupby(["release_hex", "release_doy", "beach_hex", "beach_age_bin", "shore_type"],
             as_index=False)["weight"].sum()
)
print(f"computed {len(beaching):,} rows in {time.time() - t0:.1f}s; "
      f"missing Stokes days: {stokes.missing_days}")
```

```python
beaching.to_parquet(beaching_path)
print(f"wrote {beaching_path} ({beaching_path.stat().st_size / 1e6:.2f} MB)")
```

# Validation

```python
key_ids = set(pd.read_parquet(key_path, columns=["hex_id"])["hex_id"].astype(int))
seen = set(beaching["release_hex"]) | set(beaching["beach_hex"])
unseen = seen - key_ids - {-1}
if unseen:
    print(f"WARNING: {len(unseen)} hex_ids not in {key_path.name}: "
          f"{sorted(unseen)[:10]} ...")
else:
    print(f"every release_hex/beach_hex is in {key_path.name} (or -1).")
```

```python
total = float(beaching["weight"].sum())
beached = float(beaching.loc[beaching["beach_hex"] >= 0, "weight"].sum())
by_type = beaching.groupby("shore_type")["weight"].sum()
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m")
print(f"  params: max_float_days={max_float_days}, band_m={band_m:g}, "
      f"tau0_hours={tau0_hours:g}, trap_flat/wall={trap_flat:g}/{trap_wall:g}, "
      f"w_half={w_half:g}")
print(f"  drifters (Σweight): {total:,.0f}")
print(f"  beached:           {beached:,.0f} ({100 * beached / max(total, 1):.1f}%)")
print(f"  release_doys:      {beaching['release_doy'].nunique()} "
      f"({beaching['release_doy'].min()}..{beaching['release_doy'].max()})")
print(f"  beach hexes:       {beaching.loc[beaching['beach_hex'] >= 0, 'beach_hex'].nunique():,}")
print(f"  shore_type split:  "
      + ", ".join(f"{k}={float(v):,.0f}" for k, v in by_type.items()))
beach_bins = beaching.loc[beaching["beach_age_bin"] >= 0]
if len(beach_bins):
    print(f"  beach age bins:    {beach_bins['beach_age_bin'].min()}.."
          f"{beach_bins['beach_age_bin'].max()} (× {age_bin_days} d)")
```
