---
jupyter:
  jupytext:
    formats: md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: min_data (pixi)
    language: python
    name: min_data
---

# Import

```python
import parcels
from parcels import FieldSet, JITParticle, ParticleSet, NestedField
from zarr.storage import MemoryStore

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from pathlib import Path
import geopandas as gpd

import warnings

warnings.filterwarnings("ignore", category=xr.SerializationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
```

# Parameters

```python tags=["parameters"]
# RNG seed for reproducible release-point sampling within each cell.
RNG_seed = 123

# Particle release time (ISO format). Particles enter the simulation at
# this instant; the fieldset must contain this timestamp.
start_time = "2019-01-01"

# Simulation end time (ISO format). pset.execute terminates here. The
# fieldset must contain timestamps up to end_time unless
# allow_time_extrapolation is True.
end_time = "2019-01-11"

# Release regime / kernel forcing variant: "surface", "bottom", or
# "surface_stokes". Selects which 2D-field set under
# output_root/2d_fields/ the fieldset reads from.
regime = "surface"

# Internal RK4 timestep in minutes; particle position is integrated at
# this resolution.
calc_dt_mins = 5

# Trajectory write cadence in minutes; each particle writes one
# (lon, lat) every output_dt_mins of model time.
output_dt_mins = 60

# Particles released per Fucus shapefile cell. Total particle count is
# particles_per_cell * (number of release cells).
particles_per_cell = 100

# Allow Parcels to extrapolate beyond the fieldset's time bounds. False
# (production default) raises if a particle reaches a time outside the
# loaded fields. Set True only for short verification runs against the
# one-day demo subset under data/bsh_minimal/.
allow_time_extrapolation = False

# Roots: data twin checkout (inputs) and heavy-output store (outputs).
data_root = "../data"
output_root = "../output"

# Output zarr chunking along (trajectory, obs).
chunk_traj = 10000
chunk_obs = 1000
```

# RNG

```python
rng = np.random.default_rng(RNG_seed)
print(f"RNG seed: {RNG_seed}")
```

# Setup

```python
data_root = Path(data_root)
output_root = Path(output_root)

path_2d_fields = output_root / "2d_fields"
path_trajectories = output_root / "Trajectories"

start_time = datetime.fromisoformat(start_time)
end_time = datetime.fromisoformat(end_time)
release_date_str = start_time.strftime("%Y%m%d")

calc_dt = timedelta(minutes=calc_dt_mins)
output_dt = timedelta(minutes=output_dt_mins)

print(f"start: {start_time} | end: {end_time}")
```

# Load 2D velocity fields

```python
current_files_fine = sorted(path_2d_fields.glob(f"c_file_fine_*_{regime}.nc"))
current_files_coarse = sorted(path_2d_fields.glob(f"c_file_coarse_*_{regime}.nc"))

print(
    f"Fine files: {len(current_files_fine)} | "
    f"Coarse files: {len(current_files_coarse)}"
)
```

```python
# Restrict to the time intersection of fine and coarse file stems.
fine_stems = {f.name.split("_")[3] for f in current_files_fine}
coarse_stems = {f.name.split("_")[3] for f in current_files_coarse}
common = fine_stems & coarse_stems

if common != fine_stems or common != coarse_stems:
    n_dropped = len(fine_stems) + len(coarse_stems) - 2 * len(common)
    warnings.warn(
        f"Fine/coarse timesteps not aligned: dropping {n_dropped} files "
        f"(fine: {len(fine_stems)} -> {len(common)}, "
        f"coarse: {len(coarse_stems)} -> {len(common)})"
    )

current_files_fine = sorted(
    f for f in current_files_fine if f.name.split("_")[3] in common
)
current_files_coarse = sorted(
    f for f in current_files_coarse if f.name.split("_")[3] in common
)
print(f"Timesteps: {len(common)} (fine and coarse aligned)")
```

```python
def get_timestamp_from_file(fname):
    """Derive timestamps from BSH c_file filename."""
    YYYYMMDDHH = fname.name.split("_")[3]
    y = YYYYMMDDHH[:4]
    m = YYYYMMDDHH[4:6]
    d = YYYYMMDDHH[6:8]
    h = YYYYMMDDHH[-2:]
    t0 = np.datetime64(f"{y}-{m}-{d}T{h}:15:00")
    return list(t0 + np.arange(4 * 6) * np.timedelta64(15 * 60, "s"))


timestamps_fine = [get_timestamp_from_file(f) for f in current_files_fine]
timestamps_coarse = [get_timestamp_from_file(f) for f in current_files_coarse]
```

# Fieldset

```python
def make_fieldset(data_files, timestamps, allow_time_extrapolation):
    return FieldSet.from_netcdf(
        timestamps=timestamps,
        filenames={
            "U": data_files,
            "V": data_files,
        },
        variables={
            "U": "uvel",
            "V": "vvel",
        },
        dimensions={
            "U": {"lon": "lon", "lat": "lat", "time": "time"},
            "V": {"lon": "lon", "lat": "lat", "time": "time"},
        },
        interp_method={
            "U": "cgrid_velocity",
            "V": "cgrid_velocity",
        },
        allow_time_extrapolation=allow_time_extrapolation,
        gridindexingtype="nemo",
    )


current_fieldset_fine = make_fieldset(
    data_files=current_files_fine,
    timestamps=timestamps_fine,
    allow_time_extrapolation=allow_time_extrapolation,
)
current_fieldset_coarse = make_fieldset(
    data_files=current_files_coarse,
    timestamps=timestamps_coarse,
    allow_time_extrapolation=allow_time_extrapolation,
)
```

```python
# FieldSet detects the NestedField case and auto-pairs U/V into a
# per-layer C-grid VectorField (fieldset.UV). Layer selection at the
# fine/coarse boundary happens for the (U, V) pair as a unit, so no
# fine.U / coarse.V mixing.
U_nested_field = NestedField("U", [current_fieldset_fine.U, current_fieldset_coarse.U])
V_nested_field = NestedField("V", [current_fieldset_fine.V, current_fieldset_coarse.V])
nested_fieldset = FieldSet(U_nested_field, V_nested_field)
```

# Release locations

```python
gdf_release_area = gpd.read_file(
    data_root / "derived" / "fucus_release_points.geojson"
)

n_release_cells = len(gdf_release_area)
n_total_particles = particles_per_cell * n_release_cells

print(f"release cells: {n_release_cells} | total particles: {n_total_particles}")
```

```python
def relative_position_in_cell(x_rel, y_rel, cell):
    """Map unit-square coordinates to a point inside a quadrilateral cell.

    Uses a bilinear mapping from two edge vectors of the cell. Assumes
    the cell is a parallelogram (exact) or near-parallelogram (approximate).

    Args:
        x_rel: Relative position along the first edge, in [0, 1].
        y_rel: Relative position along the second edge, in [0, 1].
        cell: A `shapely.Polygon` with 5 exterior coordinates (closed quad).

    Returns:
        Tuple of (lon, lat) in the cell's coordinate system.
    """
    (x0, y0), (x1, y1), (x2, y2), (x3, y3), (x4, y4) = cell.exterior.coords
    ex = (x3 - x0, y3 - y0)
    ey = (x1 - x0, y1 - y0)
    return x0 + x_rel * ex[0] + y_rel * ey[0], y0 + x_rel * ex[1] + y_rel * ey[1]
```

```python
# Exactly particles_per_cell particles per cell, random positions within each.
# Loop over cells (~500), broadcast the random map across particles_per_cell
# at once — relative_position_in_cell is pure numpy arithmetic and
# broadcasts over array x_rel / y_rel.
release_lons = np.empty(n_total_particles)
release_lats = np.empty(n_total_particles)
for ci, geom in enumerate(gdf_release_area.geometry):
    sl = slice(ci * particles_per_cell, (ci + 1) * particles_per_cell)
    release_lons[sl], release_lats[sl] = relative_position_in_cell(
        rng.uniform(size=particles_per_cell),
        rng.uniform(size=particles_per_cell),
        geom,
    )

print(f"particles: {len(release_lons)}")
```

# ParticleSet and output

```python
pset = ParticleSet(
    fieldset=nested_fieldset,
    pclass=JITParticle,
    lat=release_lats,
    lon=release_lons,
    time=start_time,
)
```

```python
output_filename = f"Fucus_BSH_{release_date_str}_{regime}_dt{output_dt_mins}min.zarr"
output_path = path_trajectories / regime / str(start_time.year) / output_filename
print(f"Output path: {output_path}")

memory_store = MemoryStore()
output_store = pset.ParticleFile(
    name=memory_store,
    outputdt=output_dt,
    chunks=(chunk_traj, chunk_obs),
)
```

# Execute

```python
pset.execute(
    parcels.AdvectionRK4,
    dt=calc_dt,
    endtime=end_time,
    output_file=output_store,
    verbose_progress=True,
)
```

# Write to disk

```python
ds = xr.open_zarr(memory_store)
# Drop trailing obs from the last zarr chunk that didn't fill completely.
valid_obs = ds.lon.notnull().any(dim="trajectory").compute()
ds = ds.isel(obs=valid_obs)
output_path.parent.mkdir(parents=True, exist_ok=True)
ds.chunk({"trajectory": chunk_traj, "obs": chunk_obs}).to_zarr(output_path, mode="w")
print(f"Written to {output_path} ({ds.sizes['obs']} obs)")
```

# Diagnostics

```python
ds = xr.open_zarr(memory_store)
n_valid = int(ds.lon.isel(obs=0).notnull().sum().values)
print(f"{n_valid} / {ds.sizes['trajectory']} trajectories valid at t=0")
ds
```
