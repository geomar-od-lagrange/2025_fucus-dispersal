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

# Build the beaching forcing sidecar

Cache, per real drifter and per hour, the five quantities any beaching or
survival rate model needs — and nothing else (see
[beaching_sidecar.md](../plans/done/beaching_sidecar.md)). One sidecar zarr per
trajectory zarr **at one hex radius** — the `024a` hex ids are baked in, so
a different `hex_radius` needs the sidecar rebuilt — at
`output_root/BeachingForcing/<regime>/<year>/<trajectory zarr stem>.zarr`:

| variable | dtype | meaning |
|---|---|---|
| `w_on` | float16 | onshore Stokes (m/s), 0 outside the sampled band |
| `dist` | uint8 | distance to BSH land in 25 m steps, 255 = beyond `band_max_m` / no position |
| `hex` | int32 | `024a` hex id of the position, everywhere; -1 = no position |
| `flat` | bool | nearest land is fronted by a tidal flat |
| `disp` | uint8 | crow-flies displacement from release (km), 255 = saturated |

Every one of these is **parameter-free with respect to the rate model**:
`τ`, the functional form of `s(w)`, the band width, the viability window,
the trap weights and the age binning are all reductions over these arrays,
which is what makes a rate-model sweep a seconds-per-zarr numpy pass
instead of a re-read of the trajectories and the raw wave field.

**What is sampled.**

- **distance to shore** — a rasterised distance-to-coast field built from
  the **BSH H0 land-sea mask** (finite `H0` = water, NaN = land, `H0 ≤ 0`
  = tidal flat), fine-over-coarse, in EPSG:3035. This is the mask the
  particles were advected on; the coastline geojson polygons miss a large
  fraction of genuine water positions and are not used here.
- **shore type** — the nearest land cell fronted by a tidal-flat
  (`H0 ≤ 0`) cell reads as `flat`, else `wall`. Carried as the seam for a
  real substrate classification; the BSH tidal-flat flag is not itself a
  retentiveness proxy for Baltic shores (the basin is tide-free and the
  flag fires only in the German Bight).
- **onshore wave forcing** — the onshore component of the raw
  `baltic_highres` Stokes drift (`VSDX/VSDY`), i.e. the cross-shore
  transport the `surface_stokes` runs masked at blocked faces, sampled
  here *before* that mask, and only where the position is within
  `band_max_m` of land (elsewhere it is exactly 0 and no rate model can
  use it).

**Land-seeded particles are dropped** (zero first-step displacement),
exactly as `024`/`024b` do, and the sidecar's `trajectory` coordinate
records which source trajectories survived. The key file from
`024a_BuildHexKey.md` is a hard prerequisite — its sidecar carries the
`HexProj` used to label hexes.

```python
import json
import re
import subprocess
import time
from pathlib import Path

import numcodecs
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
    r"^Fucus_BSH_(\d{8}T\d{6})_(surface_stokes|surface|bottom)_dt(\d+)min_seed\d+$"
)


def parse_zarr_stem(path):
    """Parse a trajectory zarr filename into ``(release_time, regime, dt_min)``."""
    m = _ZARR_STEM_RE.match(Path(path).stem)
    if m is None:
        raise ValueError(
            f"zarr filename does not match expected pattern: {Path(path).name!r}"
        )
    return pd.Timestamp(m.group(1)), m.group(2), int(m.group(3))
```

# Parameters

```python tags=["parameters"]
# Read root of the data twin (BSH static H0) and of the trajectory zarrs +
# raw Stokes; write root for the sidecar store.
data_root = "../data"
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

# Sampling band: onshore Stokes is sampled, and `dist` resolved, only within
# this distance of BSH land (m). 5 km is the widest band the coarse grid can
# resolve and covers 81% of particle-hours; beyond it w_on is exactly 0 and
# dist reads 255.
band_max_m = 5000.0

# Hours cached per trajectory, as 0..window_days*24 inclusive. 120 d covers
# the longest downstream horizon; shorter viability windows are a slice of
# the `obs` axis in the reducers.
window_days = 120

# Zarr output cadence (hours). The trajectory zarrs are dt60min.
output_dt_hours = 1

# Distance-to-coast raster resolution (m, EPSG:3035).
raster_dx_m = 500.0

# Max rounds of geodesic (through-water) propagation when extrapolating the
# WAM field onto BSH water, in WAM cells (~1.6 km each). Caps how far a
# sheltered cell may import wave conditions from; beyond it, w_onshore = 0.
stokes_fill_max_cells = 32

# Rebuild sidecars that already exist with matching sampling parameters.
overwrite = False
# Stop after this many zarrs of the partition (0 = all); a test knob.
max_zarrs = 0
```

# Derived layout / key + projection

```python
data_root = Path(data_root)
output_root = Path(output_root)

key_path = output_root / "HexAggregates" / f"HexAgg_key_r{hex_radius}m.parquet"
meta_path = key_path.with_suffix(".json")
if not key_path.exists() or not meta_path.exists():
    raise FileNotFoundError(
        f"Key file or sidecar missing — run 024a_BuildHexKey.md first.\n"
        f"  expected: {key_path}\n  expected: {meta_path}"
    )

meta = json.loads(meta_path.read_text())
hp = HexProj(**meta["hex_proj"])
print(f"HexProj: {meta['hex_proj']}")

forcing_root = output_root / "BeachingForcing" / regime / str(release_year)
forcing_root.mkdir(parents=True, exist_ok=True)
print(f"sidecars → {forcing_root}")

stokes_dir = output_root / "stokes" / "baltic_highres" / str(release_year)

nobs = window_days * 24 + 1
try:
    git_sha = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
except Exception:
    git_sha = "unknown"
print(f"nobs={nobs} (hourly), git_sha={git_sha}")
```

# Beaching geometry raster

Build the near-shore geometry field once: a regular EPSG:3035 raster
carrying, per cell, the distance to the nearest BSH land cell, the seaward
unit normal `n_out` (∇distance, rotated into geographic east/north), the
shore type, and the `024a` hex id of the cell centre (so per-position hex
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
    # shore type: nearest land cell adjacent to any tidal-flat cell → flat.
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
    hex_id = hex_id.reshape(nrow, ncol).astype(np.int32)

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

# Onshore Stokes sampler (land-extrapolated)

Sample the raw `baltic_highres` `VSDX/VSDY` by nearest hour and cell,
caching one day-file at a time.

**The WAM grid does not cover the BSH water mask**, and a bare nearest-cell
lookup returns NaN→0 there — which for beaching means *rate zero*, i.e. a
coastline that cannot strand under any rate model. That is a structural
bias, not a parameter choice: ~20 % (coarse) / ~34 % (fine) of near-shore
BSH water cells inside the WAM bbox sit on WAM **static land**, and WAM's
bbox (lon ≥ 9.01°E) excludes the German Bight strip entirely. So the field
is **extrapolated to full BSH coverage** by nearest-wet lookup.

The subtlety is that WAM NaN is *not* a land mask: it is land **or ice**
(the wet-cell count varies hour to hour, ~4.4 % of the grid is seasonally
ice-blanked in the Bothnian Bay / Gulf of Finland). Extrapolating across ice
would be wrong — ice genuinely suppresses waves, so zero forcing there is
the physical answer. The two are separated by a **static water mask**
(`ever_wet`: finite in *any* hour of a seasonal sample), and only static
land is filled. Transient ice keeps `w_onshore = 0`.

Fill distance is recorded per sample so the extrapolation is auditable:
most of it is trivial (median ~1 cell — the two coastlines simply disagree
by a cell), but a thin tail reaches tens of km into lagoons and fjords that
WAM does not represent at all.

```python
def bsh_water_on_grid(data_root, lon, lat):
    """BSH water mask (fine H0 over coarse) on a lon/lat grid, by FOOTPRINT.

    A cell counts as water if *any* of a 3x3 subsample spanning its footprint
    is BSH water. Sampling the centre only is wrong here: WAM cells are ~1.6 km
    and coastal ones are part water, part land, so a centre that happens to
    land on BSH land excludes a cell that particles legitimately occupy —
    which then never receives a donor and silently forces w_onshore = 0. That
    misclassified 27.9 % of particle positions (56 % of in-band samples).
    """
    h0_fine = xr.open_dataset(
        data_root / "bsh_hbmnoku_static/static_file_fine/H0_file_fine.nc"
    )
    h0_coarse = xr.open_dataset(
        data_root / "bsh_hbmnoku_static/static_file_coarse/H0_file_coarse.nc"
    )
    sample_fine = _h0_nearest_sampler(h0_fine)
    sample_coarse = _h0_nearest_sampler(h0_coarse)
    lon2, lat2 = np.meshgrid(lon, lat)
    dlon, dlat = lon[1] - lon[0], lat[1] - lat[0]
    lon_lo, lon_hi = float(h0_fine.lon.min()), float(h0_fine.lon.max())
    lat_lo, lat_hi = float(h0_fine.lat.min()), float(h0_fine.lat.max())

    water = np.zeros(lon2.shape, dtype=bool)
    for ox in (-0.5, 0.0, 0.5):
        for oy in (-0.5, 0.0, 0.5):
            lo = (lon2 + ox * dlon).ravel()
            la = (lat2 + oy * dlat).ravel()
            h0 = sample_coarse(lo, la)
            in_fine = (lo >= lon_lo) & (lo <= lon_hi) & (la >= lat_lo) & (la <= lat_hi)
            h0[in_fine] = sample_fine(lo[in_fine], la[in_fine])
            water |= np.isfinite(h0).reshape(lon2.shape)
    return water


class OnshoreStokes:
    """Raw Stokes sampler, nearest-hour/cell, extrapolated over WAM static
    land so it covers the whole BSH water mask. One-day file cache."""

    def __init__(self, stokes_dir, bsh_water_fn, fill_max_cells=32,
                 mask_sample_day=15):
        # Static mask + grid come from the release year's directory; day
        # files are resolved by the sample's own calendar year, so a window
        # running past New Year keeps its forcing (layout: <root>/<YYYY>/).
        self.stokes_dir = Path(stokes_dir)
        self.stokes_root = self.stokes_dir.parent
        self._day = None
        self._cube = None  # (VSDX, VSDY, times)
        files = sorted(self.stokes_dir.glob("stokes_*.nc"))
        if not files:
            raise FileNotFoundError(f"no raw Stokes under {self.stokes_dir}")
        g = xr.open_dataset(files[0])
        self.lon = g.longitude.values
        self.lat = g.latitude.values
        self.missing_days = 0

        # Static water mask: wet in ANY hour of a monthly sample spanning the
        # seasonal cycle, so seasonal ice does not read as land.
        by_month = {}
        for f in files:
            stem = f.stem.split("_")[-1]
            if int(stem[6:8]) == mask_sample_day:
                by_month.setdefault(stem[4:6], f)
        sample = sorted(by_month.values()) or files[:1]
        ever_wet = None
        for f in sample:
            wet = np.isfinite(xr.open_dataset(f).VSDX.values).any(axis=0)
            ever_wet = wet if ever_wet is None else (ever_wet | wet)
        self.ever_wet = ever_wet
        self.mask_files = len(sample)

        # Geodesic donor map: which WAM-wet cell each BSH-water cell reads.
        #
        # Built by breadth-first propagation of the *flat source index* — one
        # masked 4-neighbour dilation per round, so the front advances one
        # cell (~1.6 km) at a time and BSH land blocks it. Propagating indices
        # rather than values means this runs once on the static masks; per-hour
        # sampling stays a single gather.
        #
        # 4-neighbour, not 8: a 3x3 dilation squeezes between diagonally
        # touching land cells, which is exactly the thin-barrier bridging the
        # surface_stokes N=5 Stokes spread is faulted for.
        nrow, ncol = ever_wet.shape
        bsh_water_on_wam = bsh_water_fn(self.lon, self.lat)
        donor = np.where(ever_wet, np.arange(ever_wet.size).reshape(nrow, ncol), -1)
        rounds = np.full((nrow, ncol), -1, dtype="int16")
        rounds[ever_wet] = 0
        # Propagate only through BSH water; land (and non-BSH cells) block.
        allowed = bsh_water_on_wam & ~ever_wet
        for r in range(1, fill_max_cells + 1):
            todo = (donor < 0) & allowed
            if not todo.any():
                break
            src = donor
            # A plain np.roll wraps the vacated edge in from the far side of
            # the grid, which is a land-blind jump across the domain —
            # exactly what propagating through water exists to prevent. Fill
            # the vacated row/column with -1 (absent) instead, so an edge
            # cell never adopts a donor from the opposite edge.
            down = np.full_like(src, -1)
            down[1:, :] = src[:-1, :]
            up = np.full_like(src, -1)
            up[:-1, :] = src[1:, :]
            right = np.full_like(src, -1)
            right[:, 1:] = src[:, :-1]
            left = np.full_like(src, -1)
            left[:, :-1] = src[:, 1:]
            for shifted in (down, up, right, left):
                take = todo & (donor < 0) & (shifted >= 0)
                donor = np.where(take, shifted, donor)
                rounds = np.where(take, r, rounds)
        self.donor = donor
        self.fill_max_cells = fill_max_cells
        self.rounds = rounds
        self.n_unreachable = int(((donor < 0) & bsh_water_on_wam).sum())
        self.fill_rounds_used = int(rounds.max())

        dy_km = abs(self.lat[1] - self.lat[0]) * 111.32
        dx_km = abs(self.lon[1] - self.lon[0]) * 111.32 * np.cos(
            np.radians(float(np.mean(self.lat)))
        )
        # Path length along the propagation, not crow-flies distance.
        self.fill_km = (
            np.maximum(rounds, 0) * 0.5 * (dx_km + dy_km)
        ).astype("float32")
        # Diagnostics accumulated over all sampled positions.
        self.n_sampled = 0
        self.n_filled = 0
        self.n_unreachable_samples = 0
        self.n_outside_bbox = 0
        self.fill_km_sum = 0.0
        self.fill_km_max = 0.0

    def _load(self, day):
        if day == self._day:
            return self._cube
        path = self.stokes_root / str(day.year) / f"stokes_{day.strftime('%Y%m%d')}.nc"
        if not path.exists():
            self.missing_days += 1
            self._day, self._cube = day, None
            return None
        g = xr.open_dataset(path)
        self._cube = (g.VSDX.values, g.VSDY.values, g.time.values)
        self._day = day
        return self._cube

    def onshore(self, lon_q, lat_q, when, n_out_x, n_out_y):
        """Onshore Stokes magnitude max(0, -(VSDX,VSDY)·n_out) at each point.

        Positions on WAM static land (or outside the WAM bbox, which clips to
        the edge) read the nearest static-water cell. Positions on water that
        is NaN *this hour* — ice — stay at zero.
        """
        cube = self._load(pd.Timestamp(when).normalize())
        if cube is None:
            return np.zeros(lon_q.shape, dtype="float32")
        vsdx, vsdy, times = cube
        t = np.argmin(np.abs(times - np.datetime64(when)))
        i = np.round((lon_q - self.lon[0]) / (self.lon[1] - self.lon[0]))
        j = np.round((lat_q - self.lat[0]) / (self.lat[1] - self.lat[0]))
        finite = np.isfinite(lon_q) & np.isfinite(lat_q)
        inside = (
            finite & (i >= 0) & (i < len(self.lon)) & (j >= 0) & (j < len(self.lat))
        )
        i = np.clip(np.nan_to_num(i), 0, len(self.lon) - 1).astype(np.int64)
        j = np.clip(np.nan_to_num(j), 0, len(self.lat) - 1).astype(np.int64)

        # Redirect static-land (and clipped out-of-bbox) samples along the
        # geodesic donor map; ice cells are left where they are (their value
        # is NaN this hour and falls through to zero below).
        need = ~self.ever_wet[j, i]
        km = np.where(need, self.fill_km[j, i], 0.0)
        donor = self.donor[j, i]
        unreachable = need & (donor < 0)
        flat = np.where(need & ~unreachable, donor, j * len(self.lon) + i)
        sx = vsdx[t].ravel()[flat]
        sy = vsdy[t].ravel()[flat]
        onsh = -(sx * n_out_x + sy * n_out_y)
        onsh = np.where(finite & ~unreachable & np.isfinite(onsh), onsh, 0.0)

        self.n_sampled += int(finite.sum())
        self.n_filled += int((need & ~unreachable & finite).sum())
        self.n_unreachable_samples += int((unreachable & finite).sum())
        self.n_outside_bbox += int((finite & ~inside).sum())
        self.fill_km_sum += float(km[finite].sum())
        self.fill_km_max = max(self.fill_km_max, float(km[finite].max(initial=0.0)))
        return np.maximum(0.0, onsh).astype("float32")
```

# Trajectory zarr → sidecar zarr

Raster lookups run in trajectory blocks so the transient `(n, nobs)` index
arrays stay bounded; the Stokes loop then walks the hours, sampling only
the in-band positions of that hour. `dist` is quantised to 25 m steps
(`≪` the 500 m raster) with 255 reserved for "beyond the band or no
position", and `w_on` is stored as float16 — exact for the zeros that
dominate it, and well inside the error the 1.6 km WAM grid already carries.

```python
_DIST_STEP_M = 25.0
_TRAJ_BLOCK = 10000


def build_one_zarr(path, release_time, rast, stokes):
    """Write the forcing sidecar for one trajectory zarr; return its stats."""
    ds = xr.open_zarr(path)
    n_traj_source = ds.sizes["trajectory"]
    n_obs_source = min(ds.sizes["obs"], nobs)
    lon = ds.lon.isel(obs=slice(0, n_obs_source)).values.astype("float32")
    lat = ds.lat.isel(obs=slice(0, n_obs_source)).values.astype("float32")

    # Land-seeded = zero first-step displacement (as 024/024b). They sit
    # motionless on BSH land at distance 0, read as in-band every hour, and
    # dominate both the Stokes cost and the extrapolation diagnostics.
    real = (
        np.isfinite(lon[:, 0])
        & ~((lon[:, 1] - lon[:, 0] == 0) & (lat[:, 1] - lat[:, 0] == 0))
    )
    traj_index = np.flatnonzero(real).astype("int32")
    n_real = traj_index.size
    lon = lon[real]
    lat = lat[real]

    dist = np.full((n_real, nobs), 255, dtype="uint8")
    hex_at = np.full((n_real, nobs), -1, dtype="int32")
    flat_at = np.zeros((n_real, nobs), dtype=bool)
    disp = np.full((n_real, nobs), 255, dtype="uint8")
    w_on = np.zeros((n_real, nobs), dtype="float16")
    in_band = np.zeros((n_real, nobs), dtype=bool)
    n_out_x = np.zeros((n_real, n_obs_source), dtype="float32")
    n_out_y = np.zeros((n_real, n_obs_source), dtype="float32")

    cos_lat0 = np.cos(np.radians(lat[:, 0].astype("float64")))
    for b0 in range(0, n_real, _TRAJ_BLOCK):
        sl = slice(b0, min(b0 + _TRAJ_BLOCK, n_real))
        lo, la = lon[sl], lat[sl]
        shape = lo.shape
        row, col, ok = rast["to_rowcol"](lo.ravel(), la.ravel())
        ok = ok.reshape(shape)
        d_m = rast["dist_m"][row, col].reshape(shape)
        band = (d_m < band_max_m) & ok
        in_band[sl, :n_obs_source] = band
        d_idx = np.minimum(d_m / _DIST_STEP_M, 254.0).astype("uint8")
        dist[sl, :n_obs_source] = np.where(band, d_idx, 255)
        hex_at[sl, :n_obs_source] = np.where(
            ok, rast["hex_id"][row, col].reshape(shape), -1
        )
        flat_at[sl, :n_obs_source] = rast["nearest_flat"][row, col].reshape(shape) & ok
        n_out_x[sl] = np.where(band, rast["n_out_x"][row, col].reshape(shape), 0.0)
        n_out_y[sl] = np.where(band, rast["n_out_y"][row, col].reshape(shape), 0.0)
        # Equirectangular crow-flies displacement from the release position.
        dlat = la - la[:, :1]
        dlon = (lo - lo[:, :1]) * cos_lat0[sl, None]
        km = 111.0 * np.hypot(dlat, dlon)
        disp[sl, :n_obs_source] = np.where(
            np.isfinite(km), np.minimum(km, 255.0), 255.0
        ).astype("uint8")

    # Onshore Stokes at in-band positions, hour by hour (day-file cached).
    release_hours = np.datetime64(release_time) + np.arange(
        n_obs_source
    ) * np.timedelta64(output_dt_hours, "h")
    for h in range(n_obs_source):
        m = in_band[:, h]
        if not m.any():
            continue
        w_on[m, h] = stokes.onshore(
            lon[m, h], lat[m, h], release_hours[h], n_out_x[m, h], n_out_y[m, h]
        )

    release_hex = np.full(n_real, -1, dtype="int32")
    good0 = np.isfinite(lon[:, 0]) & np.isfinite(lat[:, 0])
    release_hex[good0] = hp.label(lon[good0, 0], lat[good0, 0]).astype("int32")

    out = xr.Dataset(
        {
            "w_on": (("trajectory", "obs"), w_on,
                     {"long_name": "onshore Stokes drift", "units": "m s-1"}),
            "dist": (("trajectory", "obs"), dist,
                     {"long_name": "distance to BSH land",
                      "units": f"{_DIST_STEP_M:g} m",
                      "comment": "255 = beyond band_max_m or no position"}),
            "hex": (("trajectory", "obs"), hex_at,
                    {"long_name": "024a hex id of the position",
                     "comment": "-1 = no position or outside the raster"}),
            "flat": (("trajectory", "obs"), flat_at,
                     {"long_name": "nearest land is tidal-flat-fronted"}),
            "disp": (("trajectory", "obs"), disp,
                     {"long_name": "crow-flies displacement from release",
                      "units": "km", "comment": "255 = saturated or no position"}),
        },
        coords={
            "trajectory": ("trajectory", traj_index,
                           {"long_name": "index into the source zarr trajectory dim"}),
            "release_hex": ("trajectory", release_hex,
                            {"long_name": "024a hex id of the release position"}),
        },
        attrs={
            "release_time": pd.Timestamp(release_time).isoformat(),
            "release_doy": int(pd.Timestamp(release_time).dayofyear),
            "regime": regime,
            "release_year": release_year,
            "band_max_m": float(band_max_m),
            "window_days": int(window_days),
            "raster_dx_m": float(raster_dx_m),
            "stokes_fill_max_cells": int(stokes_fill_max_cells),
            "hex_radius": int(hex_radius),
            "n_traj_source": int(n_traj_source),
            "n_obs_source": int(n_obs_source),
            "builder": "024d_BuildBeachingForcing",
            "git_sha": git_sha,
        },
    )
    # Every array gets Zstd explicitly — xarray would otherwise fall back to
    # zarr's default Blosc, which measured 2-3x larger on all of these.
    zstd = numcodecs.Zstd(level=5)
    encoding = {
        v: {"compressor": zstd, "dtype": out[v].dtype,
            "chunks": (_TRAJ_BLOCK, nobs)}
        for v in out.data_vars
    }
    encoding.update({
        c: {"compressor": zstd, "dtype": out[c].dtype, "chunks": (_TRAJ_BLOCK,)}
        for c in out.coords
    })
    target = forcing_root / f"{Path(path).stem}.zarr"
    out.to_zarr(target, mode="w", encoding=encoding)

    ncell = n_real * nobs
    return dict(
        target=target,
        n_real=n_real,
        n_traj_source=n_traj_source,
        frac_2km=float((dist < 2000.0 / _DIST_STEP_M).sum()) / ncell,
        frac_band=float((dist < 255).sum()) / ncell,
        frac_forced=float((w_on > 0).sum()) / max(float((dist < 255).sum()), 1.0),
        mb=sum(f.stat().st_size for f in target.rglob("*") if f.is_file()) / 1e6,
    )


def sidecar_is_current(target):
    """True if `target` exists and was built at today's sampling parameters."""
    if not target.exists():
        return False
    attrs = xr.open_zarr(target).attrs
    return all(
        attrs.get(k) == v
        for k, v in (
            ("band_max_m", float(band_max_m)),
            ("window_days", int(window_days)),
            ("hex_radius", int(hex_radius)),
            ("raster_dx_m", float(raster_dx_m)),
            ("stokes_fill_max_cells", int(stokes_fill_max_cells)),
        )
    )
```

```python
zarrs = sorted(
    (output_root / f"Trajectories/{regime}/{release_year}").glob("*.zarr")
)
parsed = [(p, *parse_zarr_stem(p)) for p in zarrs]
for p, ts, fn_regime, dt_min in parsed:
    assert ts.year == release_year, (ts, release_year, p)
    assert fn_regime == regime, (fn_regime, regime, p)
    assert dt_min == output_dt_hours * 60, (dt_min, output_dt_hours, p)
if release_month:
    parsed = [x for x in parsed if x[1].month == release_month]
if not parsed:
    raise FileNotFoundError(
        f"no zarrs at {output_root}/Trajectories/{regime}/{release_year}/"
        + (f" for month {release_month}" if release_month else "")
    )
if max_zarrs:
    parsed = parsed[:max_zarrs]
release_doys = sorted({int(ts.dayofyear) for _, ts, _, _ in parsed})
print(f"{len(parsed)} zarrs, release_doys "
      f"{release_doys[0]}..{release_doys[-1]} ({len(release_doys)} unique)")
```

```python
stokes = OnshoreStokes(
    stokes_dir,
    lambda lon, lat: bsh_water_on_grid(data_root, lon, lat),
    fill_max_cells=stokes_fill_max_cells,
)
print(f"WAM donor map: {stokes.fill_rounds_used} of {stokes.fill_max_cells} "
      f"propagation rounds used; {stokes.n_unreachable:,} BSH-water cells "
      f"unreachable through water (they keep w_onshore = 0)")
```

```python
t0 = time.time()
stats = []
n_skipped = 0
for p, ts, _, _ in parsed:
    if not overwrite and sidecar_is_current(forcing_root / f"{p.stem}.zarr"):
        n_skipped += 1
        print(f"  {p.stem}: sidecar up to date, skipped")
        continue
    tz = time.time()
    s = build_one_zarr(p, ts, rast, stokes)
    stats.append(s)
    print(f"  {p.stem}: {s['n_real']:,} real of {s['n_traj_source']:,}, "
          f"in-band frac {s['frac_2km']:.2f} (2 km) {s['frac_band']:.2f} (5 km), "
          f"w_on>0 frac {s['frac_forced']:.2f}, "
          f"{s['mb']:.1f} MB, [{time.time() - tz:.1f}s]")
```

```python
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m")
print(f"  band_max_m={band_max_m:g}, window_days={window_days} (nobs={nobs}), "
      f"raster_dx_m={raster_dx_m:g}, stokes_fill_max_cells={stokes_fill_max_cells}")
print(f"  {len(stats)} sidecars written, {n_skipped} up to date, "
      f"in {time.time() - t0:.1f}s; missing Stokes days: {stokes.missing_days}")
# A missing day file would silently read as zero forcing for every particle
# in band that day, so treat it as a data-coverage error, not a warning.
assert stokes.missing_days == 0, f"{stokes.missing_days} Stokes day files missing"
if stats:
    print(f"  drifters:     {sum(s['n_real'] for s in stats):,} real of "
          f"{sum(s['n_traj_source'] for s in stats):,}")
    print(f"  on disk:      {sum(s['mb'] for s in stats):,.1f} MB "
          f"({sum(s['mb'] for s in stats) / len(stats):.1f} MB per zarr)")
    print(f"  in-band frac: {np.mean([s['frac_2km'] for s in stats]):.3f} (2 km), "
          f"{np.mean([s['frac_band'] for s in stats]):.3f} (5 km)")
if stokes.n_sampled:
    print(f"WAM extrapolation ({stokes.mask_files} files in the static-water mask): "
          f"{stokes.n_filled:,} / {stokes.n_sampled:,} in-band samples filled "
          f"({100 * stokes.n_filled / stokes.n_sampled:.1f}%), "
          f"{stokes.n_unreachable_samples:,} unreachable through water, "
          f"{stokes.n_outside_bbox:,} outside the WAM bbox "
          f"({100 * stokes.n_outside_bbox / stokes.n_sampled:.1f}%); "
          f"mean fill {stokes.fill_km_sum / stokes.n_sampled:.2f} km, "
          f"max {stokes.fill_km_max:.1f} km")
```
