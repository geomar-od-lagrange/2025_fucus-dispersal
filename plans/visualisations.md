# Visualisations

Plot types as top-level sections. Scopes under each. Types without
scope sections are global / unscoped.

## Cross-cutting rules

These apply to every plot type below.

### Styling

Use the barest xarray / pandas `.plot()` method with default arguments
only. No custom figsize, dpi, line styles, colours, colormaps, norms,
alphas, or axis tweaks. Faceting goes through xarray's built-in
`col=` / `row=` arguments where it helps. Cartopy is fine for maps —
use `ax.coastlines()` with its defaults; no custom coastline overlays
(drop the BSH-shapefile coastlines for now). We'll restyle in a
second pass after seeing the defaults.

### Data pre-processing (shared loader)

- Filter particles seeded on land: drop trajectories where
  `(lon, lat).diff("obs").isel(obs=0) == 0`. Count removed and report
  the count in the global time-stats section as provenance.
- Concat all run zarrs for a given regime into one lazy dataset.

### Histogram computation (compute once, many slices)

- Use xhistogram with `obs` kept as a dim (don't aggregate obs at hist
  time). Build one lazy histogram per regime with dims
  `(lon, lat, obs, release_group)`, then every scoped plot — density,
  mean-age, per-season, per-year, per-age-bin, Baltic/German-waters
  zoom — is a `.sel` / `.sum` / `.groupby` on the lazy array. Dask
  handles the laziness; no need to persist to zarr.
- Mean-age maps reuse the same lazy object with an obs-weighted
  variant (numerator = sum over obs of `count × obs`, denominator =
  density).

### Circulation regimes (3 + 1 + 1)

Five regimes: 3 bottom runs at different velocity factors, 1 surface,
1 surface + Stokes.

- **Maps** (density, mean-age, raw trajectories): one figure per regime.
  Never co-panel bottom with surface.
- **Line / time-series plots**: regime as `hue=` so all 5 appear on
  one axes with default colours.
- **Single-number summaries** (tables, boxplots with regime on an axis):
  mixing is fine.

### Temporal scoping

Per-year and per-season plots are *alternative* slicings, never crossed.
Per-year pools across all seasons; per-season pools across all years.

### Parameters cell (papermill)

Numeric knobs in one parameters cell per viz notebook: lon/lat extents
(Baltic, German), bin sizes (both), season definition, age-bin width,
trajectory subset size, regime selector. No styling knobs.

## Notebooks

Ordered from concrete particles toward aggregated fields — zoom out
slowly. 020 and 022 run per-regime via papermill (5 executions each);
021 and 023 run once and load all regimes.

- `020_RawTrajectories.ipynb` — polyline subsets on maps. Per-regime.
  Scopes: per HELCOM release subbasin, German waters, per season,
  per year.
- `021_TimeStats.ipynb` — global ensemble diagnostics. All regimes
  shown as a regime-keyed categorical where useful. Includes the
  land-seed filter count.
- `022_DispersalDistance.ipynb` — distance-vs-time line plots,
  regimes overlaid via `hue=`. Scopes: global, per HELCOM subbasin,
  per German coast segment, per season, per year.
- `023_Heatmaps.ipynb` — density + mean-age maps, sharing the
  lazy histogram. Per-regime. Scopes: whole Baltic, per HELCOM
  subbasin, German waters, per season, per year.

Existing `030_FucusVisualization.ipynb` and
`031_FucusVisualizationSH.ipynb` move to `notebooks/explore/` as
references.

## Particle-density heatmaps

2D histogram of particle positions, one panel per scope slice, plotted
with `.plot()`. Baseline sanity check: where do particles go.

### Whole Baltic
Single panel, all trajectories and all obs pooled. Overview figure.

### Per HELCOM subbasin of release
Facet by release subbasin. Shows which source populations reach where.

### German waters
Single panel zoomed to SH + MV coast. Primary regional output.

### Per season
Facet by release season (DJF/MAM/JJA/SON), whole Baltic and
German-waters zoom. Captures seasonal circulation differences.

### Per year
Facet by release year (2016–2025), whole Baltic and German-waters zoom.
Reveals interannual variability.

## Mean-age maps

2D histogram weighted by particle age divided by count. Shows how old
particles typically are when they reach a given cell.

### Whole Baltic
Single panel. Complement to the density heatmap.

### Per HELCOM subbasin of release
Facet by release subbasin.

### German waters
Zoomed panel.

### Per season
Facet by release season.

### Per year
Facet by release year.

## Dispersal distance vs. time

Line plot of mean displacement from release point as a function of
particle age. Regime as `hue=` on every panel.

### Global
All trajectories pooled. Benchmark.

### Per HELCOM subbasin of release
One line per release subbasin (subbasin as `hue=`; regime via facet).

### Per German coast segment
One line per SH coast segment.

### Per season
Lines grouped by release season.

### Per year
Lines grouped by release year.

## Raw trajectory plots (subset)

Particle tracks drawn as lines on a map, thinned to a random subset
(subset size in parameters cell). Qualitative check that trajectories
look physical.

### Per HELCOM subbasin of release
One panel per subbasin, subset drawn from that subbasin only.

### German waters
Zoomed panel, subset drawn from German release cells only.

### Per season
Facet by release season, subset per season.

### Per year
Facet by release year, subset per year.

## Time statistics (global)

Single-figure diagnostics of the trajectory ensemble, no regional or
temporal scoping. All plotted with default `.plot()`.

- Count of particles filtered as land-seeded (from the shared loader).
- Distribution of trajectory lifetimes (max valid obs per trajectory).
- Fraction of particles still alive vs. age.
- Count of particles beached / near-coast vs. age (once beaching is
  implemented).
- Distribution of final displacements at end-of-run.
