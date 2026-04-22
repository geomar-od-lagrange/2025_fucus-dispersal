# Viz wrap-up

Simplifications across `notebooks/020_RawTrajectories`, `021_TimeStats`,
`022_DispersalDistance`, `023_Heatmaps`, and shared `helpers.py`. Driven by
(a) only one year of output exists, (b) release-point subbasin assignment is
currently limited to strict `within`, (c) a few executed plots are either
overlapping to the point of uselessness or trivially flat.

## 1. Release-point subbasin assignment (helpers.py)

Goal: every trajectory with a valid release point gets a subbasin — including
release points that fall just outside any HELCOM polygon (coast-hugging
cells). No polygon overlaps. Lazy on the `(trajectory,)` axis, never shipping
`(trajectory, obs)` lon/lat to shapely.

Approach — **raster Voronoi of HELCOM polygons, applied to release points
only**:

1. Eagerly, once on the client:
   - Rasterize `subbasins` onto a regular `(lat, lon)` grid (say 0.01° across
     the Baltic box). Integer IDs 1..N, fill=0. `rasterio.features.rasterize`.
   - Fill zero cells with their nearest non-zero neighbour:
     `idx = scipy.ndimage.distance_transform_edt(raster == 0, return_indices=True)`
     then `filled = raster[tuple(idx)]`. Cap distance so open-ocean / outside-
     Baltic stays 0.
   - This is a raster Voronoi over polygon shapes (not centroids). Overlap-
     free by construction. No geopandas buffer surgery.
2. `lon0 = ds.lon.isel(obs=0, drop=True)` / `lat0 = ds.lat.isel(obs=0, drop=True)`
   **before** any ufunc — only 1-D `(trajectory,)` arrays cross worker bounds.
3. `subbasin_id = xr.apply_ufunc(lookup, lon0, lat0, dask="parallelized",
   output_dtypes=[np.int16])` where
   `lookup(lon, lat) = filled[searchsorted(lat_edges, lat), searchsorted(lon_edges, lon)]`.
   Pure numpy per block; `filled` + edge arrays are small enough to capture by
   closure.
4. Map integer IDs back to subbasin name via a second `apply_ufunc` or a
   xarray `.map_blocks` over a pandas Index. Id 0 → NaN (land-seeded / outside
   all polygons even after fill).

Shapely is only ever touched once at the top, on N polygons — not on
trajectories. The per-trajectory pipeline is pure numpy, chunk-aligned with
`ds.lon`/`ds.lat`.

Also in helpers:
- **Drop `release_year`** from `attach_release_metadata`. Keep
  `release_quarter`. One-year output only.
- Keep `QUARTER_LABELS` / `relabel_quarter` as is.

## 2. `020_RawTrajectories.md` — join experiment types (see `explore/032`)

Current: one regime per run (papermill parameter), plots that regime across
scopes.

New: one notebook run, three regimes overlaid per panel, one color per
regime (bottom / surface / surface_stokes). Pattern from
`explore/032_QuickTrajectoryPlot.md`.

Changes:
- Drop `experiment_type` parameter. Discover regime dirs under
  `output/Trajectories/`.
- Load each regime with `load_trajectories` + `mask_land_seeded` +
  `attach_release_metadata(..., subbasins)` → `regime_dsets: dict`.
- `sample_subsets` is still per-regime; extend the panel loop to call it once
  per regime and render three `LineCollection`s per axes with distinct
  colors. Low alpha (≈0.3) and thin line widths (≈0.3) to cope with three
  overlays.
- Keep scopes: per HELCOM release subbasin, German waters, per release
  quarter. **Drop "per release year".**
- Legend: one entry per regime, drawn once per panel group (not per axes).

Keep all the structural plumbing that matters (lazy dask keys precomputed via
`dask.compute`, NaN-filtered `LineCollection`, warning silencing).

## 3. `023_Heatmaps.md` — keep regime split, drop year scope

Keep papermill `experiment_type`. One regime per run is fine for heatmaps
(three regime panels per scope would crowd the xhistogram facets).

Drop:
- `release_year` precompute, `years = sorted(...)`.
- `h_by_year_lazy`, `h_by_year` from the fused `dask.compute(...)`.
- "Per release year" section at the bottom.

Also rename any lingering "years" references. `attach_release_metadata` no
longer sets `release_year`, so referencing it would fail loudly — good.

## 4. `021_TimeStats.md` — prune dead plots

Executed-notebook inspection shows two plots are informationless under the
current data:

- **Lifetime distribution histogram** (cell 15 in executed nb): every valid
  trajectory runs the full integration, so all three histograms collapse to
  a single bar at obs=5274. `alpha=0.5` bars stack; nothing distinguishable.
  **Drop this cell and its `lifetime=ds.lon.notnull().sum("obs")` entry in
  the lazy dict.**
- **Alive fraction vs age** (cell 17): a flat line at 1.0 for all three
  regimes, for the same reason. **Drop this cell and the `alive_numerator` /
  `alive` entries in the lazy dict.**

Keep:
- Land-seed-count bar plot (one bar per regime, readable).
- Final displacement histogram. Switch from `alpha=0.5` overlapping bars to
  `histtype="step"` so the three regimes are clearly separable (bottom's
  short tail is currently half-obscured by the overlay).

Net effect: one "land-seeded count" bar + one "final displacement" stepped
histogram. Simpler notebook, nothing dropped that carried signal.

## 5. `022_DispersalDistance.md` — drop year scope

Executed-notebook inspection:
- Global / per-subbasin / German waters / per-quarter all read cleanly.
  Keep.
- **"Per release year"** facet renders a single panel labelled
  `release_year = 2.019e+03` — one year of data, one panel, already
  identical to the global panel. **Drop the section and `da_year_lazy` /
  `da_year` from the fused `dask.compute`.**

No other changes in this notebook.

## 6. Order of operations

1. `helpers.py`: rip out `release_year`; add raster-Voronoi subbasin lookup
   keyed on release coords. Verify the lookup reproduces the current
   `sjoin(within)` IDs for points inside original polygons (regression
   check on a handful of known release cells).
2. `023_Heatmaps.md`: drop year scope + cell.
3. `022_DispersalDistance.md`: drop year scope + cell.
4. `021_TimeStats.md`: drop lifetime-hist + alive-fraction cells; switch
   final-displacement to `histtype="step"`.
5. `020_RawTrajectories.md`: restructure to overlay three regimes per panel,
   drop year scope, drop `experiment_type` papermill param.
6. Re-execute all four via the jupytext → nbconvert workflow, re-inspect
   the new PNGs.

## 7. Risks / things to double-check

- Raster Voronoi uses pixel-nearest, so a release point in a narrow strait
  could snap to whichever subbasin happens to own the closer pixel at raster
  resolution. At 0.01° (~1 km) that is fine for HELCOM-scale subbasins but
  worth sanity-checking along the Danish straits.
- Integer-ID → name mapping must survive reindex over NaN (land-seeded)
  release points. Keep a sentinel id=0 → NaN.
- Overlaying three regimes in `020_RawTrajectories` triples the
  `LineCollection` vertex count per panel. Already subsampled to
  `n_traj_subset` per regime, so total ≈ 3·n; keep `n_traj_subset`
  conservative and revisit if rendering gets slow.
- `papermill` configs / scripts invoking `020_RawTrajectories` per regime
  need updating when the regime param goes away.
