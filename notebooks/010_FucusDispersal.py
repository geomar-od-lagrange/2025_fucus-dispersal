# ---
# jupyter:
#   jupytext:
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
# # Import

# %%
from parcels import FieldSet, JITParticle, ParticleSet, NestedField, Variable
from zarr.storage import MemoryStore

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from pathlib import Path
import geopandas as gpd

import warnings

warnings.filterwarnings("ignore", category=xr.SerializationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# RNG seed for reproducible release-point sampling within each cell.
RNG_seed = 123

# Particle release time (ISO format). Particles enter the simulation at
# this instant; the fieldset must contain this timestamp. BSH currents are
# on a 15-min grid with the first sample of filename-hour HH at HH:15:00,
# so doy=1 releases use T00:15:00 to land on the first available sample
# (relevant for 2016, where no prior-year file exists).
start_time = "2019-01-01T00:15:00"

# Simulation end time (ISO format). pset.execute terminates here. The
# fieldset must contain timestamps up to end_time unless
# allow_time_extrapolation is True.
end_time = "2019-01-11"

# Forcing regime, "surface_stokes" or "bottom": selects the 2D-field set
# under output_root/2d_fields/ and the wall fields the fieldset reads.
regime = "surface_stokes"

# Uniform factor on the U and V the kernel samples, at every nest level.
# 1.0 is the model current as given. Anything else is a sensitivity variant
# (e.g. 0.8 for an assumed 20 % bottom slowdown) and writes under its own
# regime name, `<regime>_vs<scale>` with "." spelled "p" (`bottom_vs0p8`),
# so it never mixes with the 1x trajectories downstream.
velocity_scale = 1.0

# Wall velocity ε (mm/s): every C-grid face that never carries flow in the
# forcing period and borders a wet cell carries ε pointing away from land,
# added to UV at every RK4 stage (docs/wall_velocity.md). Same ε in m/s on
# both nests; not multiplied by velocity_scale. 0 adds exact zeros.
wall_velocity_mm_s = 2.0

# Drop a step whose end position lies in a land cell (keep the old
# position). False with wall_velocity_mm_s = 0 is plain AdvectionRK4, bit
# for bit.
reject_land_steps = True

# Wall and land fields from 005 (walls_<nest>_<regime>.nc), relative to
# output_root or absolute.
wall_fields = "wall_fields"

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
# one-day demo subset under data/bsh_hbmnoku_demo/.
allow_time_extrapolation = False

# Roots: data twin checkout (inputs) and heavy-output store (outputs).
data_root = "../data"
output_root = "../output"

# Output zarr chunking along (trajectory, obs).
chunk_traj = 10000
chunk_obs = 1000

# %% [markdown]
# # RNG

# %%
rng = np.random.default_rng(RNG_seed)
print(f"RNG seed: {RNG_seed}")

# %% [markdown]
# # Setup

# %%
data_root = Path(data_root)
output_root = Path(output_root)

path_2d_fields = output_root / "2d_fields"
path_trajectories = output_root / "Trajectories"
path_wall_fields = output_root / wall_fields

# The name every downstream stage keys on (Trajectories/<output_regime>/).
# `regime` alone still picks the 2D-field set the fieldset reads.
output_regime = (
    regime if velocity_scale == 1.0
    else f"{regime}_vs{velocity_scale:g}".replace(".", "p")
)

start_time = datetime.fromisoformat(start_time)
end_time = datetime.fromisoformat(end_time)
release_time_str = start_time.strftime("%Y%m%dT%H%M%S")

calc_dt = timedelta(minutes=calc_dt_mins)
output_dt = timedelta(minutes=output_dt_mins)

print(f"start: {start_time} | end: {end_time}")
print(f"forcing: {regime} | velocity_scale: {velocity_scale} | output regime: {output_regime}")
print(f"wall velocity: {wall_velocity_mm_s} mm/s | reject land steps: {reject_land_steps}")

# %% [markdown]
# # Load 2D velocity fields
#
# Layout assumption: 2D-field NetCDFs sit at
# ``output_root/2d_fields/c_file_{fine,coarse}_<YYYYMMDDHH>_<regime>.nc``,
# one per BSH 6-hour window per grid per regime.

# %%
current_files_fine = sorted(path_2d_fields.glob(f"c_file_fine_*_{regime}.nc"))
current_files_coarse = sorted(path_2d_fields.glob(f"c_file_coarse_*_{regime}.nc"))

print(
    f"Fine files: {len(current_files_fine)} | "
    f"Coarse files: {len(current_files_coarse)}"
)

# %%
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


# %%
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

# %%
# The kernel's last field sample is Land, and a NestedField sample sets the
# particle state from its own result, so a time-extrapolation error raised
# by the RK4 stages would be overwritten. Check the time span up front.
last_forcing_time = min(timestamps_fine[-1][-1], timestamps_coarse[-1][-1])
assert allow_time_extrapolation or np.datetime64(end_time) <= last_forcing_time, (
    f"end_time {end_time} beyond the last forcing step {last_forcing_time}"
)


# %% [markdown]
# # Fieldset

# %%
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

# %%
# FieldSet detects the NestedField case and auto-pairs U/V into a
# per-layer C-grid VectorField (fieldset.UV). Layer selection at the
# fine/coarse boundary happens for the (U, V) pair as a unit, so no
# fine.U / coarse.V mixing.
U_nested_field = NestedField("U", [current_fieldset_fine.U, current_fieldset_coarse.U])
V_nested_field = NestedField("V", [current_fieldset_fine.V, current_fieldset_coarse.V])
nested_fieldset = FieldSet(U_nested_field, V_nested_field)

# %% [markdown]
# # Velocity scale
#
# `Field.set_scaling_factor` multiplies the field data as each time chunk is
# loaded (`rescale_and_set_minmax` in `FieldSet.computeTimeChunk`, after the
# NaN-to-zero fill), so the kernel samples scaled U and V. It is set on the
# member fields of each nest level: those are the Field objects the per-level
# C-grid VectorFields wrap and the chunk loader visits. The C-grid
# interpolation is linear in U and V, so the sampled velocity scales by
# exactly the factor. 1.0 multiplies by one, which is exact in floating point.

# %%
for field in [*U_nested_field, *V_nested_field]:
    field.set_scaling_factor(velocity_scale)

# %% [markdown]
# # Wall velocity and land
#
# The 005 fields per nest: `uvel`/`vvel` are a unit wall velocity (±1 on
# wall faces, pointing away from land), `land` is 1 on land T-cells. They
# sit on the same lon/lat as the velocity files, so each nest level gets the
# same grid and NestedField picks the same level for `UV`, `UVwall` and
# `Land`. `UVwall` is built as `FieldSet.add_UVfield` builds `UV`: a
# NestedField of per-level C-grid VectorFields. ε enters as the wall fields'
# scaling factor (time-constant fields are scaled once, at set time); the
# kernel samples `UV + UVwall` at every RK4 stage. `Land` is constant per
# cell (`cgrid_tracer`) and read at each step's end position.

# %%
def make_wall_fieldset(path):
    return FieldSet.from_netcdf(
        filenames={"Uw": path, "Vw": path, "Land": path},
        variables={"Uw": "uvel", "Vw": "vvel", "Land": "land"},
        dimensions={
            "Uw": {"lon": "lon", "lat": "lat"},
            "Vw": {"lon": "lon", "lat": "lat"},
            "Land": {"lon": "lon", "lat": "lat"},
        },
        interp_method={
            "Uw": "cgrid_velocity",
            "Vw": "cgrid_velocity",
            "Land": "cgrid_tracer",
        },
        gridindexingtype="nemo",
    )


wall_files = [path_wall_fields / f"walls_{nest}_{regime}.nc" for nest in ["fine", "coarse"]]
print("wall fields:", *wall_files)
wall_fieldset_fine, wall_fieldset_coarse = (make_wall_fieldset(p) for p in wall_files)

Uw_nested_field = NestedField("Uw", [wall_fieldset_fine.Uw, wall_fieldset_coarse.Uw])
Vw_nested_field = NestedField("Vw", [wall_fieldset_fine.Vw, wall_fieldset_coarse.Vw])
for field in [*Uw_nested_field, *Vw_nested_field]:
    field.set_scaling_factor(wall_velocity_mm_s / 1000.0)
nested_fieldset.add_field(Uw_nested_field)
nested_fieldset.add_field(Vw_nested_field)
nested_fieldset.add_vector_field(NestedField("UVwall", Uw_nested_field, Vw_nested_field))
nested_fieldset.add_field(
    NestedField("Land", [wall_fieldset_fine.Land, wall_fieldset_coarse.Land])
)
nested_fieldset.add_constant("reject_land_steps", int(reject_land_steps))


# %% [markdown]
# # Kernel
#
# `parcels.AdvectionRK4` with `UV + UVwall` at every stage. The step is
# taken unless `reject_land_steps` is set and its end position lies in a
# land cell; then the particle keeps its position and `nrej` counts the
# rejection. With ε = 0 the wall terms are exact zeros and, without the
# reject, every position is the plain RK4 one bit for bit.
#
# Land seeds are not advected at all: they start in a cell that never
# carries flow, where plain RK4 returns exactly zero, so skipping them
# changes no position at ε = 0 and saves their 9 field samples per step
# (~35 % of the particles).

# %%
def AdvectionRK4Wall(particle, fieldset, time):
    if particle.land_seed == 0:
        (u1, v1) = fieldset.UV[particle]
        (uw1, vw1) = fieldset.UVwall[particle]
        u1 = u1 + uw1
        v1 = v1 + vw1
        lon1, lat1 = (particle.lon + u1 * 0.5 * particle.dt, particle.lat + v1 * 0.5 * particle.dt)
        (u2, v2) = fieldset.UV[time + 0.5 * particle.dt, particle.depth, lat1, lon1, particle]
        (uw2, vw2) = fieldset.UVwall[time + 0.5 * particle.dt, particle.depth, lat1, lon1, particle]
        u2 = u2 + uw2
        v2 = v2 + vw2
        lon2, lat2 = (particle.lon + u2 * 0.5 * particle.dt, particle.lat + v2 * 0.5 * particle.dt)
        (u3, v3) = fieldset.UV[time + 0.5 * particle.dt, particle.depth, lat2, lon2, particle]
        (uw3, vw3) = fieldset.UVwall[time + 0.5 * particle.dt, particle.depth, lat2, lon2, particle]
        u3 = u3 + uw3
        v3 = v3 + vw3
        lon3, lat3 = (particle.lon + u3 * particle.dt, particle.lat + v3 * particle.dt)
        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3, particle]
        (uw4, vw4) = fieldset.UVwall[time + particle.dt, particle.depth, lat3, lon3, particle]
        u4 = u4 + uw4
        v4 = v4 + vw4
        dlon = (u1 + 2 * u2 + 2 * u3 + u4) / 6.0 * particle.dt
        dlat = (v1 + 2 * v2 + 2 * v3 + v4) / 6.0 * particle.dt
        land = fieldset.Land[time + particle.dt, particle.depth, particle.lat + dlat, particle.lon + dlon, particle]
        if land > 0.5 and fieldset.reject_land_steps == 1:
            particle.nrej += 1
        else:
            particle_dlon += dlon  # noqa
            particle_dlat += dlat  # noqa


# %% [markdown]
# # Release locations

# %%
gdf_release_area = gpd.read_file(
    data_root / "helcom_fucus_redlist" / "fucus_release_points.geojson"
)

n_release_cells = len(gdf_release_area)
n_total_particles = particles_per_cell * n_release_cells

print(f"release cells: {n_release_cells} | total particles: {n_total_particles}")


# %%
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


# %%
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

# %% [markdown]
# # Land seeds
#
# A particle whose release position lies in a land cell of the 005 land set
# is a land seed. `Land` is sampled through the NestedField, so the nest and
# cell are the ones the kernel uses. Written once per trajectory
# (`land_seed`); every downstream stage keys its real-drifter set on it.

# %%
land_seed = np.array(
    [nested_fieldset.Land[0.0, 0.0, la, lo] for lo, la in zip(release_lons, release_lats)]
) > 0.5
print(f"land seeds: {land_seed.sum()} of {land_seed.size} ({100 * land_seed.mean():.1f} %)")

# %% [markdown]
# # ParticleSet and output
#
# Two variables beyond position and time: `land_seed` (per trajectory) and
# `nrej`, the running count of rejected steps (per output step).

# %%
FucusParticle = JITParticle.add_variables([
    Variable("land_seed", dtype=np.int32, initial=0, to_write="once"),
    Variable("nrej", dtype=np.int32, initial=0),
])

pset = ParticleSet(
    fieldset=nested_fieldset,
    pclass=FucusParticle,
    lat=release_lats,
    lon=release_lons,
    time=start_time,
    land_seed=land_seed.astype(np.int32),
)

# %%
# Output layout: output_root/Trajectories/<output_regime>/<release_year>/<filename>.
# Downstream notebooks (020-025) discover regimes by iterdir on the
# Trajectories directory and load all zarrs underneath each regime.
output_filename = (
    f"Fucus_BSH_{release_time_str}_{output_regime}"
    f"_dt{output_dt_mins}min_seed{RNG_seed}.zarr"
)
output_path = path_trajectories / output_regime / str(start_time.year) / output_filename
print(f"Output path: {output_path}")

memory_store = MemoryStore()
output_store = pset.ParticleFile(
    name=memory_store,
    outputdt=output_dt,
    chunks=(chunk_traj, chunk_obs),
)

# %% [markdown]
# # Execute

# %%
pset.execute(
    AdvectionRK4Wall,
    dt=calc_dt,
    endtime=end_time,
    output_file=output_store,
    verbose_progress=True,
)

# %% [markdown]
# # Write to disk

# %%
ds = xr.open_zarr(memory_store)
# Drop trailing obs from the last zarr chunk that didn't fill completely.
valid_obs = ds.lon.notnull().any(dim="trajectory").compute()
ds = ds.isel(obs=valid_obs)
output_path.parent.mkdir(parents=True, exist_ok=True)
ds.chunk({"trajectory": chunk_traj, "obs": chunk_obs}).to_zarr(output_path, mode="w")
print(f"Written to {output_path} ({ds.sizes['obs']} obs)")

# %% [markdown]
# # Diagnostics

# %%
ds = xr.open_zarr(memory_store)
n_valid = int(ds.lon.isel(obs=0).notnull().sum().values)
print(f"{n_valid} / {ds.sizes['trajectory']} trajectories valid at t=0")
nrej_end = ds.nrej.where(ds.land_seed == 0).max("obs").compute()
print(
    f"land seeds: {int(ds.land_seed.sum())} | rejected steps of water seeds: "
    f"{int(nrej_end.sum())} in {int((nrej_end > 0).sum())} trajectories"
)
ds
