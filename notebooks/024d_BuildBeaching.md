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

**Rate model.** `τ = τ0 / (trap(shore_type)·s(w_onshore))`:

- **distance to shore** — a rasterised distance-to-coast field built from
  the **BSH H0 land-sea mask** (finite `H0` = water, NaN = land, `H0 ≤ 0`
  = tidal flat), fine-over-coarse, in EPSG:3035. This is the mask the
  particles were advected on; the coastline geojson polygons miss a large
  fraction of genuine water positions and are not used here.
- **shore type** — the nearest land cell fronted by a tidal-flat
  (`H0 ≤ 0`) cell reads as `flat` (dissipative, retentive), else `wall`.
  **Currently degenerate**: `trap_flat == trap_wall == 1.0`, so this term
  contributes nothing and the rate is uniform along the coast for a given
  wave forcing. The classification is still computed and carried into the
  store as `shore_type` — the seam for a real substrate classification, not
  an active model term. It is deliberately **not reported** by this notebook
  or by any consumer while `trap` is degenerate: a label that expresses
  nothing about the model invites over-reading. See the parameters cell for
  why the H0 flag is not a usable Baltic retentiveness proxy.
- **onshore wave forcing** — the onshore component of the raw
  `baltic_highres` Stokes drift (`VSDX/VSDY`), i.e. the cross-shore
  transport the `surface_stokes` runs masked at blocked faces, sampled
  here *before* that mask. With `trap` degenerate this is the *only* term
  that modulates the rate.

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
w_half = 0.05          # onshore-Stokes half-saturation (m/s) in the s ramp

# Shore-type retention weights, DELIBERATELY DEGENERATE (both 1.0) — the trap
# term is wired but currently expresses nothing, so every shore beaches alike
# for a given wave forcing. The only shore typing available here is the BSH
# `H0 <= 0` tidal-flat flag, which is not a retentiveness proxy for *Baltic*
# shores: the basin is effectively tide-free, so the flag fires only in the
# German Bight (outside the wave grid, hence zero forcing anyway). A two-class
# split would read as resolved coastal morphology while resolving none. The
# plumbing stays so a real substrate/exposure classification can drive it
# later without rebuilding the per-step lookup — set these apart to enable it.
trap_flat = 1.0
trap_wall = 1.0

# Max rounds of geodesic (through-water) propagation when extrapolating the
# WAM field onto BSH water, in WAM cells (~1.6 km each). Caps how far a
# sheltered cell may import wave conditions from; beyond it, w_onshore = 0.
stokes_fill_max_cells = 32

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
# w_half is the swept rate parameter, so it is part of the partition identity
# — sweep members must not collide. 0.05 -> "wh0p05".
wh_suffix = f"_wh{w_half:g}".replace(".", "p")
beaching_path = (
    store_root
    / f"HexAgg_beaching_r{hex_radius}m_{regime}_{release_year}{month_suffix}{wh_suffix}.parquet"
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

# Onshore Stokes sampler (land-extrapolated)

Sample the raw `baltic_highres` `VSDX/VSDY` by nearest hour and cell,
caching one day-file at a time.

**The WAM grid does not cover the BSH water mask**, and a bare nearest-cell
lookup returns NaN→0 there — which for beaching means *rate zero*, i.e. a
coastline that cannot strand at any `τ0`/`w_half`. That is a structural
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
        self.stokes_dir = Path(stokes_dir)
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
            for shifted in (
                np.roll(src, 1, 0), np.roll(src, -1, 0),
                np.roll(src, 1, 1), np.roll(src, -1, 1),
            ):
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

# Trajectory zarrs → beaching deposition

Per real drifter, deposit the fractional weight that strands at each
near-shore step (`dep = S_before − S_after`, a telescoping survival
difference), binned by the coastal hex, elapsed-age bin, and shore type of
that step. The surviving weight remaining at the window's end is the
never-beached residual, recorded once per source hex
(`beach_hex = -1`, `beach_age_bin = -1`, `shore_type = "none"`), so
downstream can form the beached fraction against the full release pool.

```python
# Realized rate diagnostics, one entry per zarr (see the s(w) block below).
RATE_STATS = []


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
    # Land-seeded particles are excluded here, not just at deposition: they sit
    # motionless on BSH land at distance 0, so they read as in-band every hour,
    # dominate the extrapolation diagnostics, and cost a third of the Stokes
    # loop — while contributing nothing (dep[~real] = 0 below).
    in_band = (dist_at < band_m) & ok.reshape(ntraj, nobs) & real[:, None]

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
    # s(w) = w/(w + w_half), the saturating wave-forcing factor. Named `s`
    # (saturation), not `g`: `g` collides with gravitational acceleration and
    # understates that this is the model's one term with no precedent in the
    # cited beaching literature (see docs/beaching.md).
    s_w = w_on / (w_on + w_half)
    a = np.where(in_band, output_dt_hours / (tau0_hours / (trap * np.maximum(s_w, 1e-6))), 0.0)
    # Realized rate diagnostics over in-band steps: the sweep reports totals,
    # but whether a member is physically sensible is read off the timescale
    # and forcing distributions, so record them rather than inferring.
    _ib = in_band & (w_on > 0)
    if _ib.any():
        _tau = tau0_hours / np.maximum(s_w[_ib], 1e-6)
        _QS = [10, 25, 50, 75, 90, 95, 99, 99.9, 100]
        _wq = np.percentile(w_on[_ib], _QS)
        _row = {'in_band_steps': int(in_band.sum()), 'forced_steps': int(_ib.sum())}
        for _q, _w in zip(_QS, _wq):
            _row[f'w_on_p{_q:g}'] = _w
            # tau AT this w quantile, not the quantile of tau: tau is monotonically
            # decreasing in w, so independent quantiles would pair opposite tails.
            _row[f'tau_h_p{_q:g}'] = tau0_hours / max(_w / (_w + w_half), 1e-6)
        RATE_STATS.append(_row)
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
stokes = OnshoreStokes(
    stokes_dir,
    lambda lon, lat: bsh_water_on_grid(data_root, lon, lat),
    fill_max_cells=stokes_fill_max_cells,
)
print(f"WAM donor map: {stokes.fill_rounds_used} of {stokes.fill_max_cells} "
      f"propagation rounds used; {stokes.n_unreachable:,} BSH-water cells "
      f"unreachable through water (they keep w_onshore = 0)")

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
if RATE_STATS:
    _rs = pd.DataFrame(RATE_STATS)
    _forced = _rs["forced_steps"].sum() / max(_rs["in_band_steps"].sum(), 1)
    print(f"realized rate over in-band steps ({100 * _forced:.1f}% with w_onshore > 0)"
          f"   [tau0={tau0_hours:g} h, w_half={w_half:g} m/s]")
    print(f"  {'quantile':>10s} {'w_onshore m/s':>15s} {'tau (h)':>12s} {'tau (d)':>10s}")
    for _c in [c for c in _rs.columns if c.startswith("w_on_p")]:
        _q = _c[len("w_on_p"):]
        _w = _rs[_c].mean(); _t = _rs[f"tau_h_p{_q}"].mean()
        print(f"  {'p' + _q:>10s} {_w:15.4f} {_t:12,.0f} {_t / 24:10,.1f}")
print(f"WAM extrapolation ({stokes.mask_files} files in the static-water mask): "
      f"{stokes.n_filled:,} / {stokes.n_sampled:,} in-band samples filled "
      f"({100 * stokes.n_filled / max(stokes.n_sampled, 1):.1f}%), "
      f"{stokes.n_unreachable_samples:,} unreachable through water, "
      f"{stokes.n_outside_bbox:,} outside the WAM bbox "
      f"({100 * stokes.n_outside_bbox / max(stokes.n_sampled, 1):.1f}%); "
      f"mean fill {stokes.fill_km_sum / max(stokes.n_sampled, 1):.2f} km, "
      f"max {stokes.fill_km_max:.1f} km")
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
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m")
print(f"  params: max_float_days={max_float_days}, band_m={band_m:g}, "
      f"tau0_hours={tau0_hours:g}, trap_flat/wall={trap_flat:g}/{trap_wall:g}"
      + (" (degenerate — shore type inert)" if trap_flat == trap_wall else "")
      + f", w_half={w_half:g}")
print(f"  drifters (Σweight): {total:,.0f}")
print(f"  beached:           {beached:,.0f} ({100 * beached / max(total, 1):.1f}%)")
print(f"  release_doys:      {beaching['release_doy'].nunique()} "
      f"({beaching['release_doy'].min()}..{beaching['release_doy'].max()})")
print(f"  beach hexes:       {beaching.loc[beaching['beach_hex'] >= 0, 'beach_hex'].nunique():,}")
beach_bins = beaching.loc[beaching["beach_age_bin"] >= 0]
if len(beach_bins):
    print(f"  beach age bins:    {beach_bins['beach_age_bin'].min()}.."
          f"{beach_bins['beach_age_bin'].max()} (× {age_bin_days} d)")
```
