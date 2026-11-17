# Dispersal-distance metric

Mean displacement from release as a function of particle age.
Computed in `notebooks/022_DispersalDistance.md` as the equirectangular
great-circle approximation `111 km · √(dlat² + (dlon · cos lat0)²)`
per trajectory, then `.mean("trajectory")`. The `cos(lat0)` factor
anchors at *release* latitude, so per-trajectory distance is a
deterministic function of `(lon, lat)` at obs and `(lon0, lat0)` at
obs=0; equirectangular vs haversine deviates <0.1 % at Baltic
latitudes.

Land-seeded particles (010's per-trajectory `land_seed`: release position in
a land cell of the regime's 005 land set, see
[wall_velocity.md](wall_velocity.md)) and NaN-padded
obs both fall out via lazy `.where`; the four scope reductions
(global, subbasin, German waters, quarter) are `.where`-masked
variants of the same lazy graph, computed in one shared
`dask.compute(*lazies)`.

## Cross-references

- [seeding.md](seeding.md) · [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) · [visualisations.md](visualisations.md)
