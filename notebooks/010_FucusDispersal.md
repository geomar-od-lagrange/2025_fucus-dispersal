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
RNG_seed = 123

# Release date (ISO format)
start_date = "2019-01-01"

# Experiment type: "surface", "bottom", or "surface_stokes"
experiment_type = "surface"

# Timesteps in minutes
calc_dt_mins = 5
output_dt_mins = 60

# What is the maximum age of the particles
max_age_days = 10

# Fraction of ambient current velocity each particle follows (1.0 = full current)
velocity_factor = 1.0
# How many particles to release per cell
particles_per_cell = 100

# Roots: inputs from data twin, heavy outputs, and (unused here) BSH store.
data_root = "../data"
output_root = "../output"

# Output chunking (trajectory, obs)
chunk_traj = 10000
chunk_obs = 1000
```

# Setup

```python
data_root = Path(data_root)
output_root = Path(output_root)

path_2d_fields = output_root / "2d_fields"
path_trajectories = output_root / "Trajectories"

release_date = datetime.fromisoformat(start_date)
release_date_str = release_date.strftime("%Y%m%d")

calc_dt = timedelta(minutes=calc_dt_mins)
output_dt = timedelta(minutes=output_dt_mins)

last_modeling_date = release_date + timedelta(days=max_age_days)

np.random.seed(RNG_seed)

print(
    "release date:",
    release_date.date(),
    "\nlast modeling date:",
    last_modeling_date.date(),
)
```

```python
output_filename = (
    f"Fucus_BSH_{release_date_str}_{experiment_type}"
    f"_dt{output_dt_mins}min_vf{velocity_factor}"
    f"_seed{RNG_seed}.zarr"
)
output_path = path_trajectories / experiment_type / str(release_date.year) / output_filename
print("Output path:", output_path)
```

# Load 2D velocity fields

```python
file_suffix = f"_{experiment_type}.nc"

current_files_fine = sorted(path_2d_fields.glob(f"c_file_fine_*{file_suffix}"))
current_files_coarse = sorted(path_2d_fields.glob(f"c_file_coarse_*{file_suffix}"))

print(
    "Fine files:",
    len(current_files_fine),
    "\nCoarse files:",
    len(current_files_coarse),
)
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


def stem(f):
    return f.name.split("_")[3]


# Restrict to time intersection of fine and coarse files
fine_stems = {stem(f) for f in current_files_fine}
coarse_stems = {stem(f) for f in current_files_coarse}
common = fine_stems & coarse_stems

if common != fine_stems or common != coarse_stems:
    n_dropped = len(fine_stems) + len(coarse_stems) - 2 * len(common)
    warnings.warn(
        f"Fine/coarse timesteps not aligned: dropping {n_dropped} files "
        f"(fine: {len(fine_stems)} -> {len(common)}, "
        f"coarse: {len(coarse_stems)} -> {len(common)})"
    )

current_files_fine = sorted(f for f in current_files_fine if stem(f) in common)
current_files_coarse = sorted(f for f in current_files_coarse if stem(f) in common)
timestamps_fine = [get_timestamp_from_file(f) for f in current_files_fine]
timestamps_coarse = [get_timestamp_from_file(f) for f in current_files_coarse]
print(f"Timesteps: {len(common)} (fine and coarse aligned)")
```

```python
# Plausibility check: verify derived timestamps against actual file contents
tolerance = np.timedelta64(1, "m")

for label, files, timestamps in [
    ("fine", current_files_fine, timestamps_fine),
    ("coarse", current_files_coarse, timestamps_coarse),
]:
    with xr.open_dataset(files[0]) as ds:
        first_actual = ds.time.values[0]
    with xr.open_dataset(files[-1]) as ds:
        last_actual = ds.time.values[-1]

    first_derived = timestamps[0][0]
    last_derived = timestamps[-1][-1]

    assert abs(first_actual - first_derived) <= tolerance, (
        f"{label}: timestamp mismatch in first file: derived {first_derived}, actual {first_actual}"
    )
    assert abs(last_actual - last_derived) <= tolerance, (
        f"{label}: timestamp mismatch in last file: derived {last_derived}, actual {last_actual}"
    )
    print(f"Timestamps {label} OK: {first_derived} ... {last_derived}")
```

# Fieldset

```python
dimension_dict = dict(lon="lon", lat="lat", time="time")
field_dimensions = [dimension_dict, dimension_dict]

current_variable_ID = ["U", "V"]
current_variable_names = ["uvel", "vvel"]
current_interp_methods = ["cgrid_velocity", "cgrid_velocity"]


def make_fieldset(data_files, variable_ID, variable_names, interp_methods, timestamps):
    data_filenames = dict(zip(variable_ID, [data_files] * len(variable_ID)))
    data_variables = dict(zip(variable_ID, variable_names))
    data_dimensions = dict(zip(variable_ID, field_dimensions))
    interp_method = dict(zip(variable_ID, interp_methods))

    return FieldSet.from_netcdf(
        timestamps=timestamps,
        filenames=data_filenames,
        variables=data_variables,
        dimensions=data_dimensions,
        interp_method=interp_method,
        allow_time_extrapolation=False,
        gridindexingtype="nemo",
    )
```

```python
current_fieldset_fine = make_fieldset(
    data_files=current_files_fine,
    variable_ID=current_variable_ID,
    variable_names=current_variable_names,
    interp_methods=current_interp_methods,
    timestamps=timestamps_fine,
)
current_fieldset_coarse = make_fieldset(
    data_files=current_files_coarse,
    variable_ID=current_variable_ID,
    variable_names=current_variable_names,
    interp_methods=current_interp_methods,
    timestamps=timestamps_coarse,
)
```

```python
U_nested_field = NestedField("U", [current_fieldset_fine.U, current_fieldset_coarse.U])
V_nested_field = NestedField("V", [current_fieldset_fine.V, current_fieldset_coarse.V])
nested_fieldset = FieldSet(U_nested_field, V_nested_field)
```

# Release locations

```python
gdf_release_area = gpd.read_file(data_root / "derived" / "fucus_release_points.geojson")

n_release_cells = len(gdf_release_area)
n_total_particles = particles_per_cell * n_release_cells
release_time = release_date + timedelta(hours=6)

print(
    "total number of particles:",
    n_total_particles,
    "\nnumber of release cells:",
    n_release_cells,
    "\nrelease time:",
    release_time,
)
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
# Exactly particles_per_cell particles per cell, random positions within each
cell_indices = np.repeat(np.arange(n_release_cells), particles_per_cell)
rand_x = np.random.uniform(size=n_total_particles)
rand_y = np.random.uniform(size=n_total_particles)

release_lons, release_lats = zip(
    *[
        relative_position_in_cell(
            rand_x[i], rand_y[i], gdf_release_area.iloc[ci].geometry
        )
        for i, ci in enumerate(cell_indices)
    ]
)

print("particles:", len(release_lons))
```

# Kernels

```python
def AdvectionRK4_2D_BSH(particle, fieldset, time):
    dt = particle.dt
    lat0 = particle.lat
    lon0 = particle.lon
    time0 = time

    u1, v1 = fieldset.UV[time0, 0, lat0, lon0]
    lat1 = lat0 + v1 * 0.5 * dt
    lon1 = lon0 + u1 * 0.5 * dt
    time1 = time0 + 0.5 * dt

    u2, v2 = fieldset.UV[time1, 0, lat1, lon1]
    lat2 = lat0 + v2 * 0.5 * dt
    lon2 = lon0 + u2 * 0.5 * dt
    time2 = time0 + 0.5 * dt

    u3, v3 = fieldset.UV[time2, 0, lat2, lon2]
    lat3 = lat0 + v3 * dt
    lon3 = lon0 + u3 * dt
    time3 = time0 + dt

    u4, v4 = fieldset.UV[time3, 0, lat3, lon3]
    lon4 = lon0 + (u1 + 2 * u2 + 2 * u3 + u4) / 6 * dt
    lat4 = lat0 + (v1 + 2 * v2 + 2 * v3 + v4) / 6 * dt

    particle_dlon += (lon4 - lon0) * particle.velocity_factor
    particle_dlat += (lat4 - lat0) * particle.velocity_factor


def max_age_kernel(particle, fieldset, time):
    particle.age_sec += particle.dt
    if particle.age_sec > particle.max_age_sec:
        particle.delete()
```

# ParticleSet and output

```python
max_age_sec = max_age_days * 24 * 60 * 60

fucus_particle = JITParticle
fucus_particle = fucus_particle.add_variable("age_sec", initial=0)
fucus_particle = fucus_particle.add_variable("max_age_sec", initial=max_age_sec)
fucus_particle = fucus_particle.add_variable("velocity_factor", initial=velocity_factor)

pset = ParticleSet(
    fieldset=nested_fieldset,
    pclass=fucus_particle,
    lat=release_lats,
    lon=release_lons,
    time=release_time,
)
```

```python
memory_store = MemoryStore()

output_particle_file = pset.ParticleFile(
    name=memory_store,
    outputdt=output_dt,
    chunks=(chunk_traj, chunk_obs),
)
```

# Execute

```python
pset.execute(
    [AdvectionRK4_2D_BSH, max_age_kernel],
    dt=calc_dt,
    endtime=last_modeling_date,
    output_file=output_particle_file,
    verbose_progress=True,
)
```

# Write to disk

```python
ds = xr.open_zarr(memory_store)
# Drop trailing obs where all particles have been deleted (NaN-padded)
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
