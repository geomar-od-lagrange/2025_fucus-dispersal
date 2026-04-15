# 2D Field Extraction Plan

Extract 2D velocity fields from the 3D BSH sigma-coordinate data. All experiment types use 2D fields — sigma level selection happens in preprocessing, not at runtime.

Reduces I/O from ~3 TB/year (25 sigma levels) to ~120 GB/year. Simplifies the custom kernel — `AdvectionRK4_2D_SIGMABSH` becomes `AdvectionRK4_2D_BSH` (sigma dimension dropped, BSH grid quirks preserved).

## `002_prepare_2d_fields.py`

Two variants:

| Variant | Content | Use case |
|---|---|---|
| `U_surface`, `V_surface` | BSH surface layer (sigma=0) | Surface experiments |
| `U_bottom`, `V_bottom` | BSH deepest valid layer per cell | Bottom experiments |

### Surface extraction

1. Extract BSH surface layer (sigma=0) as 2D
2. Write as netCDF matching BSH file structure (same time chunking, same grid dimensions)

### Bottom extraction

1. For each grid cell, identify the deepest sigma layer with valid (non-NaN) velocity — varies spatially (shallow cells have fewer layers)
2. The deepest-layer index is static (depends on bathymetry, not time) — precompute once, apply to all timesteps
3. Extract that layer's U/V as a 2D field
4. Write as netCDF matching BSH file structure

All outputs are proper C-grid 2D velocity fields.

## Simulation changes

The notebook uses a 2D fieldset with a simplified custom kernel — `AdvectionRK4_2D_BSH`. This is the current `AdvectionRK4_2D_SIGMABSH` with the sigma dimension dropped (hardcoded to 0), preserving the BSH-specific grid indexing and field ordering:

```python
def AdvectionRK4_2D_BSH(particle, fieldset, time):
    dt = particle.dt
    lat0 = particle.lat
    lon0 = particle.lon
    time0 = time

    (u1, v1) = fieldset.UV[time0, 0, lat0, lon0]
    lat1 = lat0 + v1 * 0.5 * dt
    lon1 = lon0 + u1 * 0.5 * dt
    time1 = time0 + 0.5 * dt

    (u2, v2) = fieldset.UV[time1, 0, lat1, lon1]
    lat2 = lat0 + v2 * 0.5 * dt
    lon2 = lon0 + u2 * 0.5 * dt
    time2 = time0 + 0.5 * dt

    (u3, v3) = fieldset.UV[time2, 0, lat2, lon2]
    lat3 = lat0 + v3 * dt
    lon3 = lon0 + u3 * dt
    time3 = time0 + dt

    (u4, v4) = fieldset.UV[time3, 0, lat3, lon3]
    lon4 = lon0 + (u1 + 2 * u2 + 2 * u3 + u4) / 6 * dt
    lat4 = lat0 + (v1 + 2 * v2 + 2 * v3 + v4) / 6 * dt

    particle_dlon += (lon4 - lon0) * fieldset.relative_particle_speed
    particle_dlat += (lat4 - lat0) * fieldset.relative_particle_speed
```

The experiment type (surface vs bottom) is determined by which preprocessed files are pointed at — a data concern, not a code concern.

## Implementation order

1. Create `002_prepare_2d_fields.py` — surface and bottom extraction
2. Simplify 010 to use 2D fieldset + standard `AdvectionRK4`, drop `AdvectionRK4_2D_SIGMABSH`
3. Validate with surface + bottom experiments (25x I/O reduction)
