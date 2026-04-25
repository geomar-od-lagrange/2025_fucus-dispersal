# Particle seeding

The release set is the *F. vesiculosus* presence cells in the HELCOM
REDLIST_SIS_Macrophytes shapefile.
`notebooks/000_FucusStartLocations.md` filters to `F_vesiculo != 0`,
reprojects to EPSG:4326, and bakes
`data/helcom_fucus_redlist/fucus_release_points.geojson` (committed
to the data twin; 010 reads it directly).

## Spatial: per-cell uniform sampling

Each release cell is a polygon (REDLIST EEA-grid quad), not a point.
`notebooks/010` samples `particles_per_cell` random points uniformly
inside each polygon by a bilinear map of unit-square coordinates onto
the cell's edge vectors — exact for parallelograms, near-exact for
the REDLIST quads. Uniform-in-polygon (rather than snap-to-BSH-centroid)
spreads particles across intersected BSH cells in proportion to area;
snapping would concentrate releases at the same N points regardless
of cell area, biasing the dispersal kernel.

Production runs use `particles_per_cell = 100`; with ~872 cells, each
run produces ~87 200 trajectories.

## Temporal: 73 releases per year × N years

A single Parcels run releases all particles at one instant; the sweep
is built from many such runs:

```
release_doy = 1 + 5·n   for n ∈ [0, 72]   (doys {1, 6, 11, …, 361})
```

leap-year-agnostic. Every `(year, doy, regime)` triple is one
papermill invocation of 010 producing one zarr; the sweep loops live
in `scripts/010_FucusDispersal_{surface,bottom}_job.sh`. 5-day spacing
resolves seasonal circulation transitions without bloating SLURM
wallclock.

The "single start time" simplification applies *inside* one run (no
age-tracking kernels, no per-particle kill criteria). The study
ensemble aggregates across runs, so `release_quarter` and
`release_year` remain genuine aggregation dims (020/022/023/025 facet
by quarter; 024 partitions by year).

## RNG contract

010 takes one `RNG_seed` papermill parameter and seeds a single
`numpy.random.default_rng(RNG_seed)` at the top of the run; all draws
consume from this generator and the seed is printed. Job scripts pass
a fresh seed per `(start_time, regime)` invocation, and the seed is
embedded in the output filename so reseeded resubmissions of the same
`(date, regime)` produce distinct zarrs side-by-side rather than
colliding. Multiple seeds for the same `(date, regime)` aggregate
additively on `n_obs` in the hex store — they're independent samples
of the same release distribution, useful when sample-size expansion
is wanted without changing the release date or regime.

## Output layout per run

```
output_root/Trajectories/<regime>/<year>/Fucus_BSH_<YYYYMMDD>_<regime>_dt<N>min_seed<S>.zarr
```

Downstream notebooks discover regimes / release_years via
`iterdir()` / `glob("**/*.zarr")`; multiple seed-suffixed zarrs per
`(date, regime)` are read together with no special handling.

## Cross-references

- [2d_field_extraction.md](2d_field_extraction.md) — fieldset variant
  per regime.
- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) —
  per-run trajectories aggregated into the hex store.
- [distance_calculation.md](distance_calculation.md) — distance-vs-age
  metric.
