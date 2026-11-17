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

# Wall-velocity and land fields

The time-constant fields `010` needs for its wall-velocity kernel (see
[wall_velocity.md](../docs/wall_velocity.md)), one NetCDF per nest for one
regime at `output_root/wall_fields/walls_<nest>_<regime>.nc`:

| variable | dims | meaning |
|---|---|---|
| `uvel`, `vvel` | (lat, lon) | unit wall velocity: ±1 on signed wall faces, pointing away from land; 0 elsewhere |
| `land` | (lat, lon) | 1 on non-wet T-cells, 0 on wet ones |

Same lon/lat values and order as the regime's `c_file`s, so Parcels builds
identical grids and picks the same nest for `UV`, `UVwall` and `Land`.
`010` multiplies `uvel`/`vvel` by ε as the fields' scaling factor.

**Wall face** = a C-grid face whose velocity is 0 or NaN (Parcels fills
NaN with 0) at **every** 15-min step of every `c_file` of the regime that
`010` loads (fine ∩ coarse file stems). A face that carries flow at any
step — e.g. a tidal flat that floods now and then — is not a wall.

**Wet T-cell** = at least one non-wall face (faces beyond the array count
as walls); `land` is its complement.

**Sign.** A wall face between a wet and a non-wet cell points into the wet
cell. A wall face with wet cells on both sides (a thin wall: sills between
shallow cells that flood from opposite sides) points into the deeper cell
by BSH H0; on an H0 tie into the cell with more open faces; on a second
tie it stays 0. Wall faces with no wet side stay 0.

**Index convention** on lat-ascending arrays (T-index `(j, i)` = Parcels
cell `(yi, xi) + 1`, `gridindexingtype="nemo"`): `U[j, i]` separates
`T(j, i)` (west) from `T(j, i+1)` (east), `V[j, i]` separates `T(j, i)`
(south) from `T(j+1, i)` (north). Positive U = east, positive V = north.

```python
import os
import time
from pathlib import Path

import dask.bag as db
import numpy as np
import xarray as xr
from dask.distributed import Client
```

# Parameters

```python tags=["parameters"]
# Data twin (BSH static H0) and the output root holding 2d_fields/; the
# wall fields are written to output_root/wall_fields/.
data_root = "../data"
output_root = "../output"

# 2D-field set to build the walls for ("surface_stokes" or "bottom"). 010 reads the walls of its `regime`, so `bottom_ub49` uses
# the bottom walls.
regime = "surface_stokes"

# c_files per Dask partition of the zero-face reduction.
files_per_partition = 20
```

# Setup

```python
data_root = Path(data_root)
output_root = Path(output_root)
path_2d_fields = output_root / "2d_fields"
path_wall_fields = output_root / "wall_fields"
path_wall_fields.mkdir(parents=True, exist_ok=True)
nests = ["fine", "coarse"]
```

# c_files

The file selection `010` makes: every `c_file` of the regime, restricted
to the stems both nests have.

```python
files = {n: sorted(path_2d_fields.glob(f"c_file_{n}_*_{regime}.nc")) for n in nests}
common = set.intersection(*({f.name.split("_")[3] for f in files[n]} for n in nests))
files = {n: [f for f in files[n] if f.name.split("_")[3] in common] for n in nests}
for n in nests:
    print(f"{n}: {len(files[n])} files, {files[n][0].name} .. {files[n][-1].name}")
```

# Dask cluster

```python
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
```

# Zero faces over all steps

Per file: faces 0 or NaN at all 24 steps; across files: logical AND, as a
Dask bag fold (per-partition reduce, then a tree combine).

```python
def zero_faces(path):
    """(U, V) faces that are 0 or NaN at every step of one c_file, on-disk order."""
    with xr.open_dataset(path) as ds:
        uv = ds[["uvel", "vvel"]]
        z = ((uv == 0) | uv.isnull()).all("time")
        return z.uvel.values, z.vvel.values


def both_and(a, b):
    return a[0] & b[0], a[1] & b[1]


zero = {}
for n in nests:
    t0 = time.time()
    zero[n] = (
        db.from_sequence(files[n], partition_size=files_per_partition)
        .map(zero_faces)
        .fold(both_and)
        .compute()
    )
    print(f"{n}: zero faces from {len(files[n])} files in {time.time() - t0:.0f} s")
```

# H0 at the T-cells

The H0 file carries one row and column more than the `c_file`s; T-cell
`(j, i)` is the H0 point nearest to `(lon[i] − Δlon/2, lat[j] − Δlat/2)`.
Used only to orient thin walls.

```python
def h0_tcells(nest, lon, lat):
    """H0 at T(j, i) on the lat-ascending (lat, lon) grid of the c_files."""
    h0 = xr.open_dataset(
        data_root / "bsh_hbmnoku_static" / f"static_file_{nest}" / f"H0_file_{nest}.nc"
    ).H0
    dlon, dlat = float(np.diff(lon).mean()), float(np.diff(lat).mean())
    hs = h0.sel(
        lon=xr.DataArray(lon - dlon / 2, dims="lon"),
        lat=xr.DataArray(lat - dlat / 2, dims="lat"),
        method="nearest",
    )
    assert np.abs(hs.lon.values - (lon - dlon / 2)).max() < abs(dlon) / 4
    assert np.abs(hs.lat.values - (lat - dlat / 2)).max() < abs(dlat) / 4
    return hs.values
```

# Signed walls

Per T-cell the four face flags (faces beyond the array are walls); T-row
and T-column 0 are not Parcels cells and never wet.

```python
def signed_walls(zu, zv, h0):
    """Unit wall velocity (uw, vw), wet mask and counts from lat-ascending zero-face masks."""
    zu_w = np.pad(zu, ((0, 0), (1, 0)), constant_values=True)  # zu_w[j, i] = zu[j, i-1]
    zv_s = np.pad(zv, ((1, 0), (0, 0)), constant_values=True)  # zv_s[j, i] = zv[j-1, i]
    w_west, w_east = zu_w[:, :-1], zu_w[:, 1:]
    w_south, w_north = zv_s[:-1, :], zv_s[1:, :]
    wet = ~(w_west & w_east & w_south & w_north)
    wet[0, :] = False
    wet[:, 0] = False
    n_open = 4 - (w_west.astype(int) + w_east + w_south + w_north)
    wet_east = np.pad(wet[:, 1:], ((0, 0), (0, 1)), constant_values=False)  # T(j, i+1)
    wet_north = np.pad(wet[1:, :], ((0, 1), (0, 0)), constant_values=False)  # T(j+1, i)

    uw = np.zeros(zu.shape)
    vw = np.zeros(zv.shape)
    uw[zu & wet_east & ~wet] = +1.0
    uw[zu & wet & ~wet_east] = -1.0
    vw[zv & wet_north & ~wet] = +1.0
    vw[zv & wet & ~wet_north] = -1.0

    # Thin walls: side 1 = T(j, i), side 2 = the east / north neighbour.
    counts = dict(thin_h0=0, thin_open=0, thin_tie=0)
    for arr, thin, shift in [
        (uw, zu & wet & wet_east, (0, 1)),
        (vw, zv & wet & wet_north, (1, 0)),
    ]:
        for j, i in zip(*np.nonzero(thin)):
            j2, i2 = j + shift[0], i + shift[1]
            if h0[j2, i2] != h0[j, i]:
                arr[j, i] = 1.0 if h0[j2, i2] > h0[j, i] else -1.0
                counts["thin_h0"] += 1
            elif n_open[j2, i2] != n_open[j, i]:
                arr[j, i] = 1.0 if n_open[j2, i2] > n_open[j, i] else -1.0
                counts["thin_open"] += 1
            else:
                counts["thin_tie"] += 1

    counts.update(
        wall_faces=int(zu.sum() + zv.sum()),
        signed=int((uw != 0).sum() + (vw != 0).sum()),
        no_wet_side=int((zu & ~wet & ~wet_east).sum() + (zv & ~wet & ~wet_north).sum()),
        wet_cells=int(wet.sum()),
    )
    assert (
        counts["signed"] + counts["thin_tie"] + counts["no_wet_side"] == counts["wall_faces"]
    )
    return uw, vw, wet, counts
```

# Write

```python
summary = {}
for n in nests:
    with xr.open_dataset(files[n][0]) as ref:
        coords = {"lat": ref.lat.values, "lon": ref.lon.values}
    zu, zv = zero[n]
    z = xr.Dataset({"zu": (("lat", "lon"), zu), "zv": (("lat", "lon"), zv)}, coords=coords)
    z = z.sortby("lat")
    uw, vw, wet, counts = signed_walls(
        z.zu.values, z.zv.values, h0_tcells(n, z.lon.values, z.lat.values)
    )
    out = xr.Dataset(
        {
            "uvel": (("lat", "lon"), uw, {"long_name": "unit wall velocity, eastward",
                                         "units": "m s-1"}),
            "vvel": (("lat", "lon"), vw, {"long_name": "unit wall velocity, northward",
                                         "units": "m s-1"}),
            "land": (("lat", "lon"), (~wet).astype("float32"),
                     {"long_name": "land T-cell (no face carries flow at any step)"}),
        },
        coords={"lat": z.lat.values, "lon": z.lon.values},
    ).reindex(lat=coords["lat"])  # back to the c_files' on-disk lat order
    out.attrs.update(
        regime=regime,
        source=f"{len(files[n])} c_files, {files[n][0].name} .. {files[n][-1].name}",
        description="wall faces = 0 or NaN at every step; +-1 into the wet cell; "
        "thin walls into the deeper cell (H0), tie -> more open faces",
    )
    path = path_wall_fields / f"walls_{n}_{regime}.nc"
    out.to_netcdf(path)
    summary[n] = counts
    print(f"{n:6s} {regime}: {counts} -> {path}")
```

# Summary

```python
for n, c in summary.items():
    print(
        f"{n:6s} {regime}: {c['wall_faces']:,} wall faces, {c['signed']:,} signed "
        f"(thin walls by H0 {c['thin_h0']}, by open faces {c['thin_open']}, "
        f"left 0 {c['thin_tie']}), {c['no_wet_side']:,} with no wet side, "
        f"{c['wet_cells']:,} wet T-cells"
    )
```
