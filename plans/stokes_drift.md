# Stokes Drift Plan

Depends on: `2d_field_extraction.md` (2D surface fields must exist first)

## Data Source

**BALTICSEA_MULTIYEAR_WAV_003_015** — Baltic Sea Wave Hindcast (FMI, WAM cycle 4.7)

| Property | Value |
|---|---|
| Variables | `VSDX`, `VSDY` (surface Stokes drift, m/s) |
| Resolution | 2 km, hourly |
| Coverage | 53–66°N, 9–30°E (full Baltic + Danish Straits + Kattegat) |
| Period | 1980–2026 (covers our 2016–2025 window) |
| Grid | Regular lat/lon (A-grid) |
| Access | `copernicusmarine` Python toolbox |

The global WAVERYS product (0.2°) is too coarse — can't resolve Baltic straits.

## Preprocessing

Extend `002_prepare_2d_fields.py` to produce a third variant:

| Variant | Content |
|---|---|
| `U_total`, `V_total` | BSH surface + Stokes |

Steps:
1. Download VSDX/VSDY from CMEMS for the required time range (~240 GB compressed for 2016–2025)
2. Fill NaN/land in Stokes fields with zero before regridding — where BSH has water but the wave model doesn't (narrow fjords, sheltered waters), Stokes drift is negligible anyway. Physically defensible, no extrapolation needed.
3. Interpolate VSDX onto BSH **U-points** (eastern cell edges), VSDY onto BSH **V-points** (northern cell edges)
4. Zero out at edges where either adjacent cell is land — ensures no flux through coastline
4. Sum: `U_total = U_surface + U_stokes`, `V_total = V_surface + V_stokes`
5. Write as netCDF matching BSH file structure

The result is a proper C-grid 2D velocity field. Land boundaries respected by construction. No divergence concerns — particles are surface-pinned (buoyant Fucus).

No simulation code changes needed — same 2D fieldset, same `AdvectionRK4_2D_BSH` kernel. Just point at the `U_total`/`V_total` files instead of `U_surface`/`V_surface`.

## Near-shore Tracking and Beaching

### Coast distance field

Create a static 2D `coast_distance` field:
1. Compute from BSH bathymetry (H0 files) — distance transform on the land mask
2. One field per grid (fine + coarse)
3. Static, loaded once — no I/O cost per timestep

### Near-shore tracking kernel

```python
def nearshore_kernel(particle, fieldset, time):
    d = fieldset.coast_distance[time, particle.depth, particle.lat, particle.lon]
    if d < fieldset.nearshore_threshold_km:
        particle.nearshore_time += particle.dt
```

Additional particle variable: `nearshore_time` (cumulative seconds within threshold distance of coast).

### Post-processing: Stochastic beaching model

Applied in visualization notebooks, not in the simulation.

The beaching weight for heatmaps and connectivity can be a function of:
- `nearshore_time` — total near-shore exposure
- Location along coast (sheet pile walls vs sandy beaches)
- Stokes drift magnitude at the coast (wave energy → beaching probability)
- Season (ice cover, storm frequency)

Different beaching models can be tested on the same trajectory dataset without re-running simulations.

## Implementation Order

1. Download Stokes drift data from CMEMS
2. Extend `002` to regrid Stokes onto C-grid edges, land-mask, and sum with BSH surface
3. Create static coast distance field (from H0 bathymetry)
4. Add `nearshore_time` tracking kernel to 010
5. Add beaching weight logic to 020/021 post-processing
