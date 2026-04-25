# Pre-prod wrapup

Final pass before the repo and study go public. This plan is a living
punch-list: as more remarks land in `remarks.md` and inline TODOs across
the tree, fold them into the appropriate section below. Once everything
here is closed, write the matching `docs/*.md` files, move this plan to
`plans/done/`, and we're prod-ready.

Sources of remarks folded in so far:

- `remarks.md` (top-level, 5 items as of 2026-04-25)
- Inline TODOs in `notebooks/010_FucusDispersal.ipynb` (16 items)
- Inline TODOs in `notebooks/020_RawTrajectories.ipynb` (13 items)
- Inline TODOs in `notebooks/helpers.py` (2 items, including a policy
  reversal that drives §3 below)
- Inline FB on `docs/bundle_and_layout.md`
- A two-stage agent pass (A: theme rubric from the annotations above;
  B: rubric applied to un-annotated notebooks `000`, `021`–`025`)
  whose findings populate §5 below.

More remarks may still land; extend the sections below in place
rather than starting a new plan.

## 1. Particle-physics simplification (010)

The single-start-time regime makes several historical contraptions
unnecessary. Strip them, then verify a fresh sweep still produces
sensible trajectories.

- [x] Remove `velocity_factor` end-to-end: kernel multiplication,
      particle attribute, papermill parameter, sweep loops in
      `scripts/010_FucusDispersal_*_job.sh`, downstream filename
      conventions, and any aggregation code keyed on it. The
      bottom-focused experiment will get its own kernel handling
      bottom slowdown — no need to keep a generic per-particle scale.
- [x] **New zarr filename schema**. Drop `vf{velocity_factor}` and
      rename `{experiment_type}` → `{regime}`. Keep only the fields
      `parse_zarr_stem` needs downstream: release date and regime.
      Landed: `Fucus_BSH_{release_date_str}_{regime}_dt{output_dt_mins}min.zarr`
      — three fields, one parsed by date, one by regime token, one
      retained for run-config traceability. `parse_zarr_stem` in
      `helpers.py` updated to match (will be inlined in §3).
- [x] Remove `max_age_kernel` and the `age_sec` / `max_age_sec`
      particle attributes. With a common `release_date`, end-time on
      `pset.execute` is the only kill criterion needed.
- [x] Replace `last_modeling_date = release_date + timedelta(days=max_age_days)`
      with explicit `start_time` / `end_time` parameters using names
      consistent with `pset.execute`'s kwargs. Drop the legacy "+6h"
      shift; let the fieldset-bounds exception raise if a user picks a
      start time that collides with available data. (Also added
      `allow_time_extrapolation` papermill parameter, default `False`,
      for short verification runs against the demo subset.)
- [x] Revisit whether the custom `AdvectionRK4_2D_BSH` can be replaced
      by Parcels' standard 2D RK4 once `velocity_factor` is gone. If
      yes, delete the custom kernel. (Resolved: yes; replaced with
      `parcels.AdvectionRK4`.)
- [x] Once both kernels are gone, audit particle attribute schema for
      anything else now unused (`age_sec`, etc.).
- [x] **Scope of "single start time"**: the simplification is
      *per-execution*, not *per-study*. Production design is 73
      releases/year × N years of papermill-injected 010 runs at
      `release_doy = 1 + 5*n` for `n ∈ [0, 72]` (doys
      `{1, 6, 11, …, 361}`, leap-year-agnostic), each producing one
      zarr with one `release_date`, one `release_quarter`, one
      `release_year`. Downstream notebooks aggregate across runs, so
      `release_quarter` and `release_year` remain genuine aggregation
      dimensions. Conclusion: `max_age_kernel`, `age_sec`,
      `last_modeling_date`, and the custom RK4 still go (per-run
      simplifications); per-quarter facets in `022_DispersalDistance`
      / `023_Heatmaps` / `025_HexHeatmaps` and the `release_year`
      partition in `024_BuildHexAggregates` stay (cross-run
      aggregation).
- [x] Update `scripts/010_FucusDispersal_*_job.sh` sweep loops to
      `n ∈ [0, 72]` (73 doys ≤ 361). Drop any `doy == 366` branch
      from sweep generators.

## 2. Notebook 010 readability cleanup

Smaller hygiene items flagged inline. Do these in the same pass as §1
so the rerun covers both.

- [x] Comment every parameter in the parameters cell (010 done; §5
      propagates the convention to the other notebooks).
- [x] Move `np.random.seed(RNG_seed)` to the first cell after params,
      and unify the seeding idiom across notebooks (single
      `np.random.default_rng(seed)` per run, printed for
      reproducibility).
- [x] Collapse the release-date / last-modelling-date `print(...)` to
      a single line, placed immediately after the dates are defined.
- [x] Move construction of `output_filename` to *after* fieldset
      creation; it belongs adjacent to the `ParticleFile` /
      particleset block, not next to date arithmetic.
- [x] Rename `output_particle_file` → `output_store` (it writes to a
      `MemoryStore`, not a file).
- [x] Inline `file_suffix`.
- [x] Inline `def stem(f)`.
- [x] Split the timestamp/stem cell: first cell computes the common
      stems across coarse / fine groupings and *warns* if they don't
      fully overlap; second cell builds the `np.timedelta64` array.
      (Static-grid alignment not added — current 010 doesn't read
      static files.) Also dropped the redundant timestamp
      plausibility check that re-read file metadata to verify the
      precomputed timestamps; trust the derivation and let the
      FieldSet raise if it's wrong.
- [x] Expand `make_fieldset` body — type the `data_filenames`,
      `data_variables`, `data_dimensions` dicts out explicitly instead
      of three layered `dict(zip(...))` calls. Verified
      `data_filenames = dict(zip(variable_ID, [data_files] * len(variable_ID)))`
      is correct: it produces `{"U": data_files, "V": data_files}`.
      `allow_time_extrapolation` plumbed through as a parameter.
- [x] Rewrite the `release_lons, release_lats = zip(*[...])`
      comprehension explicitly — replaced with a cell-major loop that
      broadcasts `relative_position_in_cell` across `particles_per_cell`
      per cell. ~872 iterations instead of ~87K, shapely accessed once
      per cell. `rng.uniform` calls inlined in the loop body.
      Smoke-tested at production particle count: release-points cell
      runs in <50 ms.
- [x] Fix the NaN-trim comment: particles aren't "deleted" — the
      trailing zarr chunk just didn't fill up completely.

Additionally landed during Phase A but not in the original plan:

- [x] Added `parcels>=3,<4` and pinned `zarr>=2.18,<3` in
      `pixi.toml` (parcels was missing from the project env; the
      former `min_data` kernelspec pointed at a deleted submodule
      env). Both `pixi.toml` and `pixi.lock` updated.
- [x] Added a one-line comment to the `FieldSet(U_nf, V_nf)` block
      explaining that `_add_UVfield` re-pairs the two scalar
      `NestedField`s into a per-layer C-grid `VectorField`, so layer
      selection at the fine/coarse boundary happens for the (U, V)
      pair as a unit.
- [x] Updated `AGENTS.md` pipeline-stages section: two job scripts
      (`surface`, `bottom`), drops `velocity_factor` from the sweep
      tuple, names the `<regime>/<year>/` output layout.
- [x] Smoke-test verified at production particle count
      (87,200 trajectories, 3 days): notebook runs end-to-end against
      the demo subset (`data/bsh_hbmnoku_demo/`) with
      `allow_time_extrapolation=True`. Filename, layout, and kernel
      attribute (`JITParticleAdvectionRK4`) all match the new schema.

## 3. Eliminate `notebooks/helpers.py`

**Policy reversal.** The previous rule ("extract to `helpers.py` only
when the logic is non-trivial and reused verbatim") is dropped. The new
rule: notebooks define their own helpers locally or inline the logic.
Nothing currently in `helpers.py` is non-trivial enough to justify a
shared module once the sub-bullets below are addressed.

- [x] Simplify `attach_release_metadata`'s subbasin lookup. **Function
      deleted entirely** — split into two inline pieces per consumer
      notebook (020/022/023): a one-line `release_quarter` assignment
      from `ds.time.dt.quarter`, and a notebook-local
      `assign_release_subbasin` using shapely `STRtree.nearest` kept
      lazy via `xr.apply_ufunc(..., dask="parallelized")`. Eager
      `gpd.sjoin_nearest` was tried first and rejected — would OOM at
      production scale (60M+ trajectories).
- [x] Audit every callsite of every `helpers` import. Resolved as
      follows: `load_trajectories` and `mask_land_seeded` → trivial
      inline (no def) at each callsite; `attach_release_metadata` →
      split per above; `relabel_quarter` → deleted (panel titles
      render as `release_quarter = 1..4`; section headers carry the
      JFM/AMJ/JAS/OND legend); `parse_zarr_stem` → local `def` in
      024 (single caller, regex earns the name); `QUARTER_LABELS` →
      4-entry dict literal at the loop site in 020 and 025.
- [x] `relabel_quarter` deleted; `QUARTER_LABELS` inlined as
      per-notebook dict literal `quarter_labels = {1:"JFM", …}` in
      020:283 and 025:290.
- [x] Deleted `notebooks/helpers.py` and `notebooks/__pycache__/`.
- [x] Updated `AGENTS.md` "Notebook-local utilities" subsection
      (replaces "Shared helpers"). Now four explicit gating tests:
      (1) shared module not forbidden but bar is really high — three
      copies of a 12-line function isn't enough; (2) idiom test
      before def-ing; (3) one concern per function; (4) don't mutate
      coords for presentation.

Phase B execution scope vs. Phase G full-sweep verification: Phase B
landed the helpers replacement and a smoke-test against existing
demo zarrs is left for Phase G (mixed-format zarrs from pre-Phase A
runs would confuse `parse_zarr_stem`; cleanest verify is against a
freshly-swept output tree).

## 4. Notebook 020 readability cleanup

The biggest single cleanup target after 010. Several items here close
naturally once §3 is done; do §3 first.

- [x] Rename extents with consistent prefixes:
      `(lon_min, lon_max, lat_min, lat_max)` → `baltic_*`;
      `(de_lon_min, …)` keeps `de_*`. Both prefix groups, no bare
      versions.
- [x] Drop the hardcoded `regime_colors = {"bottom": "tab:orange",
      "surface": "tab:blue"}`. Take colours from the default cycle by
      regime index, or use a colourblind-safe palette mapped by index
      (no `"tab:..."` literals in the notebook).
- [x] Move all imports to the top of the notebook. The late
      `import os, time, dask.distributed, ...` cell collapses into the
      header import block.
- [x] Configure warning filters once at the top of the notebook, not
      with per-cell `with warnings.catch_warnings()` context managers.
- [x] Add a one-line comment to the
      `regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())`
      line stating the assumed
      `output_root/Trajectories/<regime>/...` layout — so a failure
      points the reader at the layout assumption rather than at the
      autodetect mechanics.
- [x] Drop the `regime_keys` precompute loop. Filter the dataset where
      the plot is built; xarray + dask will fuse the access.
- [x] Inline `lonlat_aspect` (one-liner) at its sole callsite.
- [x] Inline `plot_lines` into the per-regime figure loop. The
      `with warnings.catch_warnings():` block disappears with §4 above.
- [x] Drop the `subbasins_list` precompute and the empty-subbasin
      branch (`if avail.size == 0: ax.set_visible(False); continue`).
      Plot every subbasin every regime; empty panels are fine and keep
      the layout stable across regimes. `rng.choice(avail,
      size=min(n_traj_subset, avail.size), replace=False)` already
      handles `avail.size == 0` (returns empty array, plots nothing).
- [x] Drop the legend; use `fig.suptitle(regime)` for context (and
      `set_title(subbasin)` per panel as already done).

## 5. Per-notebook walk (`000`, `021`–`025`)

Per-notebook punch-lists from the rubric pass. Each item is concrete
and cell-anchored. Items that overlap §1–§4 / §7–§8 conventions are
not duplicated here unless the notebook has a distinctive instance.

### 5a. `notebooks/000_FucusStartLocations.md`

- [x] Comment the lone `data_root` parameter (role: read root of the
      data twin checkout).
- [x] Add a layout-assumption comment naming the
      `helcom_fucus_redlist/REDLIST_SIS_Macrophytes.shp` path the
      notebook reads.
- [x] Revisit the `data_root / "derived"` write target once §7's
      `data/` rename lands — `derived/` is being dropped as a
      category. (Resolved: writes to
      `helcom_fucus_redlist/fucus_release_points.geojson`,
      co-located with its source shapefile.)

### 5b. `notebooks/021_TimeStats.md`

- [x] Comment the `output_root` parameter.
- [x] Hoist the late `# Dask cluster` cell (`import os`, `import time`,
      `from dask.distributed import Client`) into the top import
      block.
- [x] Replace `from helpers import load_trajectories,
      mask_land_seeded` with notebook-local definitions (per §3).
      (No-op: Phase B already removed; verified no `helpers` import.)
- [x] Rewrite the
      `lazy = dict(...); results = dict(zip(lazy.keys(), dask.compute(*lazy.values())))`
      block with five named lazies and an explicit unpack — the
      dict-of-lazies + zip-rehydrate hides the fan-in.
- [x] Add the `output_root/Trajectories/<regime>/...` layout comment
      next to `regimes = sorted(p.name for p in
      trajectory_root.iterdir() if p.is_dir())`.
- [x] Collapse the multi-line per-regime f-string `print(...)` to one
      line, or move the reporting into a small dedicated summary
      cell.

### 5c. `notebooks/022_DispersalDistance.md`

- [x] Comment every parameter in the parameters cell.
- [x] Rename Baltic-side extents to `baltic_*` so the prefix
      convention matches `de_*` (currently the Baltic side is bare
      and the DE side is prefixed — same lie as 020).
      (No-op: 022 has no Baltic-side extent parameters; only `de_*`.
      Nothing to rename.)
- [x] Hoist the late `# Dask cluster` cell.
- [x] Replace `from helpers import attach_release_metadata,
      load_trajectories, mask_land_seeded, relabel_quarter` per §3.
      (No-op: Phase B already removed.)
- [x] Inline `def in_de_mask(ds)` (one expression, one callsite).
      Keep `distance_km` (called twice, modestly substantive).
- [x] Inline `def load_regime(regime)` into the dict-comp or rewrite
      as an explicit loop.
- [x] Rewrite `def _grouped_mean(group_key)` without the closure
      over `regimes` / `regime_dsets` / `regime_distance`; either
      pass the dependencies explicitly or write the two `xr.concat`
      blocks out. (Wrote both out as explicit `xr.concat` blocks.)
- [x] Inside `_grouped_mean`, drop the
      `d.assign_coords({group_key: ds[group_key].compute()})` —
      forces the group coord eager when the rest of the graph is
      lazy. Pre-compute once outside or keep it lazy. (Pre-computed
      once per regime in `regime_group_coords` before the concat.)
- [x] Add the layout-assumption comment next to the regimes
      autodetect.

### 5d. `notebooks/023_Heatmaps.md`

- [x] Comment every parameter.
- [x] Rename `experiment_type = "surface"` → `regime = "surface"`.
      Downstream use is exclusively as a regime
      (`trajectory_root / experiment_type`); the name is a leftover
      from an earlier vocabulary. Match `024` / `025`.
- [x] Rename Baltic-side extents to `baltic_*` (`lon_min/lat_min`,
      `n_lon_baltic` already mixes prefixes — pick one).
- [x] Hoist the late `# Dask cluster` cell.
- [x] Replace `from helpers import attach_release_metadata,
      load_trajectories, mask_land_seeded, relabel_quarter` per §3.
      (No-op: Phase B already removed.)
- [x] Inline `def lonlat_aspect(extent)` (duplicate of the 025
      definition; trivial). (Inlined the two top-level callsites and
      the body inside `single_map`.)
- [x] Inline `def count_hist(ds_, lon_bins, lat_bins)` (one-expression
      wrapper). Keep `density` and `mean_age_hours` (two callsites
      each over different histograms).
- [x] Drop the precompute block
      `sb_np, quarter_np = dask.compute(ds.subbasin, ds.release_quarter)`
      + the `subbasins_list` / `quarters` distinct-set construction.
      Build histograms with lazy `.where(key_da == v)` over a known
      list (subbasins from the `subbasins` GeoDataFrame; quarters as
      `[1, 2, 3, 4]`).
- [x] Replace the
      `sorted({s for s in sb_np if isinstance(s, str)})` /
      `sorted({int(q) for q in quarter_np if not np.isnan(q)})`
      idiom with a direct dropna+`np.unique`. (Subsumed by the
      previous item: precompute gone, lists derived from the static
      inputs directly — `subbasins["subbasin"].tolist()` and
      `[1, 2, 3, 4]`.)
- [x] Add the layout-assumption comment next to
      `trajectory_path = output_root / "Trajectories" / experiment_type`.

### 5e. `notebooks/024_BuildHexAggregates.md`

- [x] Comment every parameter (`data_root`, `output_root`, `bsh_root`,
      `release_year`; `age_bin_days` and `output_dt_mins` already
      have inline comments).
- [x] Replace `hex_configs = [...]` with a single scalar
      `hex_radius` parameter. Multiple radii come from papermill
      sweeps over 024, not from an in-notebook loop. Hard-code
      projection and centering — equal-area projection centred on the
      BSH domain centroid — in a code cell below params. The
      per-config loop disappears entirely; `_h0_hex_frame` and
      `build_counts` become straight-line code without a closure
      (next bullet then drops to a one-line "no rewrite needed").
      (Domain centroid computed from the coarse-grid H0 bbox.)
- [x] Output layout change: `output_root/HexAggregates/<config_name>/`
      → `output_root/HexAggregates/r{hex_radius}{unit}/` (pick the
      unit when implementing — m, km, or whatever the hex library
      expects). Update the partition-layout comment block accordingly.
      (Picked metres — matches the HexProj `hex_size_meters` API
      exactly; `r{hex_radius}m` keeps param value and path in lockstep.)
- [x] Hoist the late `from dask.distributed import Client` import.
- [x] Replace `from helpers import load_trajectories, mask_land_seeded,
      parse_zarr_stem` per §3. Inline `parse_zarr_stem` against the
      simplified filename schema from §1 (date + regime; no `vf{}`).
      (No-op: Phase B already inlined `parse_zarr_stem` and dropped
      the other two helper imports.)
- [x] Inline `def _zarr_for(regime, ...)` (three-line glob wrapper,
      one callsite).
- [x] `_h0_hex_frame` and `build_counts`: with the per-config loop
      gone, both either become inline code in the main cell or stay
      as plain functions with explicit args (no closure trick, no
      `arg_=value` defaults). The current comment acknowledging the
      closure trap goes with the closure. (Both kept as plain
      functions with explicit args; renamed `_h0_hex_frame` →
      `h0_hex_frame` since it's no longer a private closure helper.)
- [x] Simplify the
      `subbasin_id_to_name = {-1: "_outside"}` /
      `subbasin_name_to_id = {v: k for k, v in subbasin_id_to_name.items() if k >= 0}`
      round-trip: build name→id from `enumerate` directly, derive the
      reverse from it.
- [x] Add a single comment block at the top of the path-construction
      cell listing the four layout assumptions it encodes
      (`bsh_root/static_file_<grid>/H0_file_<grid>.nc`,
      `output_root/HexAggregates/<config>`,
      `output_root/Trajectories/<regime>/<release_year>`,
      and the `counts/regime=.../release_year=...` partition layout).
- [x] Drop the `release_doy == 366` skip outright. With the sweep
      bounded to `n ∈ [0, 72]` (max doy 361, see §1), the branch is
      dead code, not leap-year handling.
- [x] Tighten the per-config / output-sizes summary prints; consider a
      small DataFrame display instead of multi-line f-strings.
      (Per-regime summary now a small `pd.DataFrame` display; the
      output-sizes block kept as f-string table — concise enough.)

### 5f. `notebooks/025_HexHeatmaps.md`

- [x] Comment every parameter (only `cmap` and the panel-height pair
      currently are).
- [x] Rename Baltic-side extents to `baltic_*` (same prefix lie as
      022/023).
- [x] Replace `from helpers import QUARTER_LABELS` with a local
      `QUARTER_LABELS = {1: "JFM", 2: "AMJ", 3: "JAS", 4: "OND"}`.
      (No-op: Phase B already inlined as `quarter_labels` literal.)
- [x] Justify the `cmap = "viridis"` parameter and the per-overlay
      style overrides (`color="black"`, `color="magenta"`,
      `linewidth=...`, `edgecolor="face"`) inside `log_density_plot`
      in `docs/visualisations.md` (per §6) rather than as inline
      comments. Where defaults read fine, drop the override.
      (Per-locked-decision: cmap stays; deferred per-override
      justify/strip to Phase F via `# TODO(phase-f)` comments —
      Phase C scope fence forbids §6 docs work.)
- [x] Inline `def lonlat_aspect(extent)` at its two callsites.
- [x] Drop the `subbasins_ordered` precompute in Panel B and the
      empty-panel `set_visible(False)` skips in Panels B and D. Empty
      panels render fine and keep layout stable. (Iterates every
      named subbasin in stable id order; removed the
      `set_visible(False)` skip in Panel D.)
- [x] **Panel D**: per-quarter facet is load-bearing across the
      cross-run aggregation (73 releases/year × N years span all
      quarters). Keep, but clean per the §4-style simplifications
      (drop empty-panel skips, lazy `.where`).
- [x] Replace the doy→quarter Timestamp arithmetic with the one-liner
      `((doy - 1) // 90 + 1).clip(1, 4)` (or month-from-doy →
      `(month - 1) // 3 + 1`). Add a comment naming the conversion.
      (See judgment-call note in the Phase C report — `// 90` and
      `// 91` are both off by ≥1 quarter for some sweep doys, so
      went with calendar-correct `pd.to_datetime(format="%Y%j").dt.quarter`,
      a single timestamp call instead of the original two-step
      `pd.to_datetime + pd.to_timedelta`.)
- [x] Decide on the four `groupby...rename...reset_index().merge(...)
      .pipe(gpd.GeoDataFrame, ...)` chains: extract one
      notebook-local `to_hex_gdf(counts_df, key_df)` (one helper, four
      callsites — justified) or accept the duplication. (Extracted
      per locked decision.)
- [x] Add the layout-assumption comment naming
      `output_root/HexAggregates/<baltic_config_name>` and the
      `counts/regime=…/release_year=…` partition layout.

## 6. Methodology docs (replace `docs/bundle_and_layout.md`)

`bundle_and_layout.md` documents plumbing that decays as the layout
moves. Drop it. Replace with focused method docs — one per scientific
or methodological decision a reader of the published study needs.

The criterion for "this deserves a doc": a future reader of the paper
who reads only the code can't reconstruct the *why*.

**Source-material audit of `plans/` (open and `plans/done/`)** — most
candidate docs have an existing plan as draft. Each bullet names the
target `docs/` file and the plan(s) to mine; the plans then move to
`plans/done/` (or stay there) with a one-liner pointing at the new
doc.

- [x] `docs/h0_semantics.md` — BSH H0 sign convention (z-up,
      MSL-zero, tidal flats can be `H0 < 0`); always-wet mask.
      **Source**: `plans/seafloor-location-H0-semantics.md` is
      already in docs-format — practically a move-and-rename. Trim
      the duplicate H0 paragraph in CLAUDE.md down to a one-line
      pointer at the doc.
- [x] `docs/2d_field_extraction.md` — surface and bottom 2D-field
      extraction from the 3D BSH sigma grid; I/O reduction (~3 TB/yr
      → ~120 GB/yr); why sigma-layer selection happens at preprocess
      rather than at runtime; the resulting file structure consumed
      by 010. Add one sentence noting that the 2D-field strategy
      sidesteps the live-3D-interpolation failure modes (corner
      trapping, sigma-coord registration ambiguity) that the
      archived `plans/done/corner_*.md` and
      `plans/done/grid_registration*.md` plans investigated; point
      readers there for the historical record. **Source**:
      `plans/done/2d_field_extraction.md` plus the bottom-cell
      zero-velocity finding from `plans/bottom_stationary_audit.md`
      (the audit *finding*, not the audit *process*).
- [x] **Not** a separate doc: the T-point/F-point grid-registration
      story and the C-grid corner-trapping story are *historical*.
      Both were investigated against the 3D-sigma pipeline retired
      by the "Go 2d" refactor (commit `5ebd46b`, 2026-04-15); their
      fixes were never wired into the current 2D-field code path.
      Leave `plans/done/grid_registration*.md` and
      `plans/done/corner_*.md` archived as research record. The
      one-sentence motivation pointer in
      `docs/2d_field_extraction.md` is enough; documenting the
      fixes themselves would mislead future readers into thinking
      they ship.
- [x] `docs/stokes_drift.md` — CMEMS Baltic Wave Hindcast
      (`BALTICSEA_MULTIYEAR_WAV_003_015`, variables `VSDX`/`VSDY`,
      2 km hourly), why not the global WAVERYS, regridding onto the
      BSH grid, summation with Eulerian currents. **Source**:
      `plans/done/stokes_drift.md` is most of the way there;
      modernise the access recipe to match
      `notebooks/002_download_stokes.py`.
- [x] `docs/seeding.md` — release-point sourcing from the Fucus
      shapefile; per-cell uniform sampling; n-particles bookkeeping;
      RNG contract; the 73-releases-per-year × N-years sweep design
      (doys `1 + 5n` for `n ∈ [0, 72]`, see §1). **Source**: derive
      from 000 + 010 source after §1/§5a are done — no plan covers
      this end-to-end yet.
- [x] `docs/hexbinning_and_connectivity.md` — unified hex grid for
      source and target, equal-area projection centred on Baltic,
      `key.parquet` + `counts/` schema, the symmetric source↔sink
      self-join that makes `ostrea`-style queries work.
      **Source**: `plans/hex_aggregate_store.md` is already
      methodology-shaped — practically a rewrite-as-doc. Cross-link
      with §5e's `hex_radius` parameterisation.
- [x] `docs/distance_calculation.md` — distance-vs-time metric
      definition (great-circle vs cumulative path vs from-release),
      edge cases (NaN-padded obs, land-seeded particles excluded).
      **Source**: derive from 022 source after §5c is done — no
      plan covers this.
- [x] `docs/visualisations.md` — per-plot-type rationale: what each
      notebook (020–025) shows, scope decisions, where overrides
      are justified. Styling *rules* stay in CLAUDE.md.
      **Source**: trim `plans/visualisations.md` to the
      per-plot-type sections; promote those to docs.
- [x] After the new docs land: delete `docs/bundle_and_layout.md`. Do
      *not* leave a stub or redirect; git history is the changelog.

### Plans that stay in `plans/done/` as historical record

These are process / cleanup / governance notes, not methodology — no
`docs/` promotion:

- `plans/data_licensing_public_bundle.md` (compliance audit, lands
  here once green per §9).
- `plans/experiment_tracking.md` (superseded by §1's filename
  simplification — archive once the new schema ships).
- `plans/viz_wrap_up.md` (subsumed by §4/§5 — archive once §4/§5
  close).
- `plans/done/010_notebook_cleanup.md`, `bundle_and_layout.md`,
  `job_script_review.md`, `next_steps.md`, `notebook_review.md`,
  `portable_data_paths.md`, `rollout.md`, `simplification.md`,
  `hextraj_hex_counts_oom.md` (already archived; leave in place).
- `plans/done/grid_registration.md`,
  `plans/done/grid_registration_bug.md`,
  `plans/done/corner_rounding.md`, `plans/done/corner_theory.md`,
  `plans/done/corner_vortex.md`, `plans/done/corner_*.png` —
  research record for the 3D-sigma pipeline retired by the
  "Go 2d" refactor; the investigated fixes were never shipped.
  Pointed at from `docs/2d_field_extraction.md` motivation
  paragraph; otherwise leave in place.

## 7. Repo surface presented to public users

These are the first things a visitor sees; tidy before announcing.

- [x] Rename `data/` subdirs to `<source>_<dataset>` form (source
      first, then specific dataset). Walk every existing dir and
      decide on its merits — don't preserve a name just because
      something already uses it. Drop `data/derived/` as a category;
      either co-locate the derived blob with its source dir or give
      it its own clearly-named dir.
- [x] Update consumers in lock-step: every script / notebook reading
      `data/...` paths, the `obtain/*.sh` recipes, `fetch_data.sh`,
      and `ATTRIBUTION.md`. Grep before, grep after. (`helpers.py`
      is gone after §3, so no path constants live there any more.)
- [x] Mirror the rename in the data twin repo
      (`git.geomar.de/od-lagrange/2025_fucus_dispersal_data`) and
      bump the submodule pointer in the same PR.
- [x] `README.md`: introduce the term "twin" explicitly when first
      mentioning the data submodule, so newcomers see the project
      vocabulary before they encounter it elsewhere.
- [x] `scripts/000_FucusStartLocations_job.sh`: drop. Stage 000
      doesn't need a job script — it's quick local prep. Delete the
      `.sh`, leave the `.md` / `.ipynb` notebook in place.
- [ ] Walk every other `scripts/0??_*_job.sh`: confirm each is still
      load-bearing, names match the current notebook stage numbering,
      `output_root` plumbing is consistent, no stale `--region` /
      `--year` flags from the old filename-encoded era. Rewrite where
      cleaner than patching.

## 8. Conventions to enforce across the whole pipeline

These are the rules §5's per-notebook items make concrete. Restate them
here as a final cross-notebook audit so nothing slips through.

- [x] **Parameters cell**: every parameter commented, primitives only,
      RNG seeding (`np.random.default_rng(seed)`) as the first
      post-params cell where the notebook uses randomness. Notebooks
      that don't use randomness skip the seed cell rather than adding
      a no-op one. (Verified in 000/010/020/021/022/023/024/025.
      Fix: 020's `output_root` comment said "write root for derived
      stores" — corrected to "Read root of trajectory zarrs".)
- [x] **Imports at the top**: no late `import …` cells. The repeated
      `# Dask cluster` block in 021/022/023/024 is the recurring
      offender — hoist its imports. (No-op: §5 hoists already landed.
      Re-grep confirmed every `import`/`from` in 0??_*.md sits in the
      header import block.)
- [x] **Warning filters at the top**: configure once near the imports;
      no `with warnings.catch_warnings()` blocks scattered through
      cells. If a notebook currently has none, leave it alone.
      (010 and 020 set filters at the top; rest have none. No
      `with warnings.catch_warnings()` anywhere.)
- [x] **No `helpers` import** anywhere (§3 removes the file; this is
      the cross-cutting check). (Re-grepped: zero hits.)
- [x] **Consistent extent prefixes**: pairs of related extent
      variables both prefixed (`baltic_*` / `de_*`); no bare
      `lon_min` / `lat_min` paired with prefixed `de_lon_min`. Hits
      currently in 020, 022, 023, 025. (No-op: §5 renames already
      landed. Function-local trailing-underscore unpacks
      `lon_min_, lon_max_, lat_min_, lat_max_ = extent` in
      `single_map`/`log_density_plot` are scope-local and don't
      pair with module-scope vars.)
- [x] **Variable names match current behaviour**: notable rename to
      apply: `experiment_type` → `regime` in `023_Heatmaps`. Look for
      similar lies elsewhere (e.g. `*_file` for in-memory stores,
      `last_modeling_date` for end-time semantics). (No-op: §5d
      already renamed. Re-grep finds no further lies; the remaining
      `*_file` matches are all real file paths or parcels API kwargs.)
- [x] **No styling overrides without a doc-justified reason**: no
      `color="tab:..."` literals, no `cmap=` / `figsize=` / `linewidth=`
      overrides except where `docs/visualisations.md` explicitly
      argues for them. (Phase F gates the doc-justify pass; per the
      Phase E fence, new unjustified overrides get `# TODO(phase-f)`
      markers, not deletions. Marked: 020 line 211 trajectory-line
      `linewidth=0.5`. Existing 025 markers untouched. Aspect-driven
      `figsize=(panel_height_in*aspect, panel_height_in)` in
      020/023/025 is treated as accepted layout sizing — 025's
      figsize calls also carry no TODO.)
- [x] **Layout assumptions are commented**: any `iterdir()` /
      glob walk / hard-coded sub-path encodes a layout — name the
      assumption in a one-liner so failures point at the contract,
      not the autodetect mechanics. Hits in 020, 021, 022, 023, 024,
      025. (Existing markdown blocks above each `iterdir` cover
      020/021/022/023/024. Added: 010 `# Load 2D velocity fields`
      glob and 010's output-path layout block now name the contracts
      explicitly.)
- [x] **Notebooks read standalone**: parameter cells, path
      construction, and the Dask-cluster boilerplate are *expected* to
      be duplicated per notebook. Do not re-factor the
      `SCHEDULER_FILE`-or-local-`Client` snippet into a helper — that
      is the test case for §3's policy reversal. If it needs
      documenting, add a snippet to `docs/` rather than a Python
      module. (Verified: cluster boilerplate is duplicated verbatim
      in 020/021/022/023/024 and helpers.py is gone.)

## 9. Final verification before flipping to prod

- [ ] Re-execute all `notebooks/0??_*.md` end-to-end via
      `pixi run jupytext --sync --execute …` (papermill where
      parameters need injecting). Commit the freshly-rendered
      `.ipynb` alongside the `.md` so GitHub renders figures.
- [ ] Smoke-test from a fresh clone:
      `git clone --recurse-submodules` → `pixi install` →
      `scripts/fetch_data.sh` no-op → run stages 000 through 025
      against the BSH demo subset.
- [ ] `ATTRIBUTION.md` walkthrough: every dataset currently shipped
      in the twin still listed; no licence-incompatible additions
      since the last audit (cross-check against
      `plans/data_licensing_public_bundle.md`).
- [ ] Review `AGENTS.md` shape after F + G have landed. Look for
      content trimmable to pointers (rules-vs-reference split),
      `Conventions › Notebooks` bullets that duplicate the jupytext
      skill, and stale references after the rename + docs extraction.
      Decide whether a structural pass earns its keep or the file
      reads fine post-trim. No new sub-plan unless the answer is yes.
- [ ] Once green: move this file to `plans/done/wrapup.md` with a
      one-liner pointing at the new method docs as the durable
      record.

## 10. Implementation gaps discovered during prod-prep

- [x] **Stokes drift WAVERYS fallback for German Bight.** The
      Baltic high-res CMEMS Stokes product (`cmems_mod_bal_wav_my_PT1H-i`)
      starts at ~9 °E; the BSH fine grid extends west to ~6.2 °E. Pre-
      fix, `interpolate_stokes` zeroed the strip and `surface_stokes`
      particles in the German Bight silently degraded to BSH-only.
      Implementation: `notebooks/002_download_stokes.py` now pulls both
      Baltic high-res and WAVERYS (`cmems_mod_glo_wav_my_0.2deg_PT3H-i`)
      side by side under `<output-root>/stokes/<product>/<year>/`;
      `notebooks/003_prepare_2d_fields.py` `interpolate_stokes` layers
      Baltic-where-defined / WAVERYS-as-fallback. Each wave-model
      field is spread by N=5 iterations of a 3×3 rolling mean before
      interp (covers 100 % of BSH potentially-wet cells in the
      wave-model file extent for both grids; fixed empirically against
      the BSH fine + coarse domains). Stokes is shut off per-timestep
      where BSH says the face is blocked (`u_surf == 0` or NaN) so
      tidal flats receive Stokes when wet and zero when dry — likely
      source of the "near-shore slowdown" beaching artifact in earlier
      runs. Verified end-to-end against the demo (BSH 2020-01-01 +
      Baltic high-res + WAVERYS): 003 produces no-NaN 2D fields with
      Stokes contribution in every lon band 5–15 °E; 010 advects 872
      particles over 4 h with sensible displacements (max ~13 km).

## Open questions (collect here as they arise)

- **2026-04-25 — Stokes spread bridges thin land barriers.** The N=5
  3×3 rolling-mean spread in `interpolate_stokes` extends ~10 km on
  the Baltic 2 km wave-model grid, wide enough to push open-ocean
  Stokes across barriers like the Curonian Spit (~1–3 km) into the
  sheltered Curonian Lagoon — unphysical (real fetch resets at the
  spit). Three candidate mitigations laid out in
  [docs/stokes_drift.md](../docs/stokes_drift.md) §"Open concern".
  Kept at N=5 + per-timestep face mask for now; needs further thought
  before settling on a fix.

### Resolved

- **2026-04-25 — "Single start time" scope.** Per-execution, not
  per-study. Production design: 73 releases/year × N years at
  `release_doy = 1 + 5n` for `n ∈ [0, 72]` (doys `{1, 6, …, 361}`,
  leap-year-agnostic). Cross-run aggregation keeps `release_quarter`
  and `release_year` as load-bearing dimensions; quarter facets and
  year partitions stay. The `release_doy == 366` skip becomes dead
  code under the bounded sweep and is dropped. Folded into §1, §5e.
- **2026-04-25 — `hex_configs` shape.** Don't move the list-of-dicts
  into the parameters cell. Replace with a *scalar* `hex_radius`
  parameter; multiple radii come from papermill sweeps over 024, not
  an in-notebook loop. Hard-code equal-area projection centred on the
  BSH domain centroid in a code cell below params. The per-config
  loop disappears, `_h0_hex_frame` and `build_counts` lose their
  closure tricks, and the output layout becomes
  `HexAggregates/r{hex_radius}{unit}/`. Folded into §5e.
- **2026-04-25 — Output filename schema.** Simplify. Drop
  `vf{velocity_factor}`, rename `{experiment_type}` → `{regime}`,
  retain only what `parse_zarr_stem` actually needs (release date +
  regime + a stable `dt{output_dt_mins}min` for traceability).
  Folded into §1, §5e.
