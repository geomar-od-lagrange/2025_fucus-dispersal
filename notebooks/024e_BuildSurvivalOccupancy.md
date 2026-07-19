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

— the same near-shore beaching rate `τ = τ0/(trap·g(w_onshore))` as
`024d_BuildBeaching` (see [beaching.md](../docs/beaching.md)), with `trap`
currently degenerate (`trap_flat == trap_wall == 1.0`), so wave forcing is
the only term that modulates it. `A` only grows inside the near-shore band,
so open-water residence is undiluted; a particle that lingers in a
wave-exposed near-shore band loses weight fast. This is the
deterministic occupancy analogue of 024d's fractional stranding: 024d
records where the weight *leaves* (`beach_hex`), this records where the
still-drifting weight *is* (`target_hex`).

Emits **two** weights per bin so the fold-in is a self-consistent
comparison: `occ` (plain, every sample = 1) and `surv` (Σ S) over the same
window and hexing. `030_SurvivalHeatmaps` maps `surv`, `occ`, and the
surviving fraction `surv/occ`.

Reuses 024d's raster + Stokes machinery; aggregates with `np.bincount`
(occupancy is dense — every obs, every hex — unlike 024d's sparse in-band
deposits). Partitioned per `(regime, year, month)`; `030` pools. The
`024a` key is a hard prerequisite (its `HexProj` + hex-id set).

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
w_half = 0.05

# Shore-type retention weights, DELIBERATELY DEGENERATE (both 1.0) — see the
# 024d parameters cell. The trap term is wired but expresses nothing, so the
# beaching rate is uniform along the coast for a given wave forcing; the
# `flat`/`wall` classification is computed and carried only as the seam for a
# future real substrate/exposure dataset.
trap_flat = 1.0
trap_wall = 1.0

# Max rounds of geodesic (through-water) propagation when extrapolating the
# WAM field onto BSH water, in WAM cells (~1.6 km each).
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
# w_half is a rate parameter, so it is part of the partition identity — runs
# at different w_half must not collide. 0.05 -> "wh0p05".
wh_suffix = f"_wh{w_half:g}".replace(".", "p")
survocc_path = (
    store_root
    / f"HexAgg_survocc_r{hex_radius}m_{regime}_{release_year}{month_suffix}{wh_suffix}.parquet"
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
```

# Beaching geometry raster (+ dense hex index per cell)

Same field as 024d — distance to the BSH H0 coast, seaward normal rotated
to geographic east/north, wall/flat shore type — plus, here, the dense hex
**index** (into the key's hex-id order) of each raster cell, so per-position
hexing is a lookup ready for `np.bincount`.

```python
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

# Trajectory zarrs → survival-weighted occupancy

Per real drifter, weight every in-window occupancy sample by `S = exp(−A)`;
`A` accrues only inside the near-shore band. Aggregate `occ` (=1) and `surv`
(=S) over `(target_hex, age_bin)` with `np.bincount`.

```python
# Realized rate diagnostics, one entry per zarr (see the s(w) block below).
RATE_STATS = []


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
    # Land-seeded particles are excluded here, not just at aggregation: they
    # sit motionless on BSH land at distance 0, so they read as in-band every
    # hour, dominate the extrapolation diagnostics, and cost a third of the
    # Stokes loop — while being dropped from the occupancy sums anyway.
    in_band = (dist_at < band_m) & ok & real[:, None]

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
        _wq = np.percentile(w_on[_ib], [10, 50, 90])
        _tq = np.percentile(_tau, [10, 50, 90])
        RATE_STATS.append({
            'in_band_steps': int(in_band.sum()),
            'forced_steps': int(_ib.sum()),
            'w_on_p10': _wq[0], 'w_on_p50': _wq[1], 'w_on_p90': _wq[2],
            'tau_h_p10': _tq[0], 'tau_h_p50': _tq[1], 'tau_h_p90': _tq[2],
        })
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
if RATE_STATS:
    _rs = pd.DataFrame(RATE_STATS)
    _forced = _rs["forced_steps"].sum() / max(_rs["in_band_steps"].sum(), 1)
    print(f"realized rate over in-band steps ({100 * _forced:.1f}% with w_onshore > 0):")
    print(f"  w_onshore m/s  p10/p50/p90: {_rs['w_on_p10'].mean():.4f} / "
          f"{_rs['w_on_p50'].mean():.4f} / {_rs['w_on_p90'].mean():.4f}")
    print(f"  tau hours      p10/p50/p90: {_rs['tau_h_p10'].mean():,.0f} / "
          f"{_rs['tau_h_p50'].mean():,.0f} / {_rs['tau_h_p90'].mean():,.0f}"
          f"   (tau0={tau0_hours:g} h, w_half={w_half:g} m/s)")
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
survocc.to_parquet(survocc_path)
print(f"wrote {survocc_path} ({survocc_path.stat().st_size / 1e6:.2f} MB)")
```

# Validation

```python
key_ids = set(hex_ids.tolist())
unseen = set(survocc["target_hex"]) - key_ids
if unseen:
    print(f"WARNING: {len(unseen)} target_hex not in key: {sorted(unseen)[:10]} ...")
else:
    print("every target_hex is in the key.")
```

```python
print(f"regime={regime}, release_year={release_year}"
      + (f", month={release_month}" if release_month else "")
      + f", hex_radius={hex_radius} m")
print(f"  params: occupancy_max_days={occupancy_max_days}, band_m={band_m:g}, "
      f"tau0_hours={tau0_hours:g}, trap_flat/wall={trap_flat:g}/{trap_wall:g}"
      + (" (degenerate — shore type inert)" if trap_flat == trap_wall else ""))
per_bin = survocc.groupby("age_bin")[["occ", "surv"]].sum()
per_bin["drifting"] = per_bin["surv"] / per_bin["occ"]
print(per_bin.to_string(float_format=lambda v: f"{v:,.3f}"))
```
