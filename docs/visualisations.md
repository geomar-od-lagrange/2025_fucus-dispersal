# Visualisations: per-plot rationale

Why each notebook in `020`–`031` shows what it shows, and where
styling overrides earn their deviation from
[../AGENTS.md](../AGENTS.md)'s plotting rules. AGENTS.md states the
**rules** (no `cmap=`, no `figsize=`, no `color=` literals); this doc
states the **rationale** for each surviving rule-break. Anything not
listed here should be dropped.

## Plot types and scopes at a glance

| Notebook | Plot type | Per-regime? | Default scopes |
|----------|-----------|-------------|----------------|
| 020 | Raw trajectory lines | Overlaid | Per HELCOM release subbasin · German waters · Per release quarter |
| 021 | Time-series stats | Overlaid | Global only |
| 022 | Distance vs. time | Overlaid (`hue="regime"`) | Global · Per subbasin · German waters · Per quarter |
| 023 | Density + mean-age maps | Per-run (papermill regime) | Whole Baltic · Per subbasin · German waters · Per quarter |
| 025 | Hex relative-density maps | Per-run (regime + season + radius) | Baltic percentage+dilution · German-waters zoom · per origin subbasin |
| 026 | Hex time-horizon relative-density maps | Per-run (regime + season + radius) | Cumulative elapsed-time horizons (10/20/50/100 d), all release years pooled |
| 026a | 026 per origin subbasin | Per-run (regime + season + radius) | One percentage and one dilution figure per origin (HELCOM) subbasin |
| 026b | 026a per origin subbasin and year | Per-run (regime + season + radius) | One percentage and one dilution figure per (origin subbasin, release year) |
| 027 | Hex distance-quantile maps | Per-run (regime + radius + season) | One 1xN panel row, one panel per quantile (0.1/0.5/0.9), the season's releases pooled across years |
| 028 | Subbasin connectivity matrices | Per-run (regime + radius + season + member) | Raw residence · emission fraction per age horizon · interannual range |
| 029 | Beaching maps | Per-run (regime + radius + member + season) | Where-stranded · beached fraction per source hex · age horizons (10/20/50 d) · beached-fraction curve with interannual band · travel distance at stranding |
| 030 | Survival heatmaps | Per-run (regime + radius + member + season) | Per horizon (20/50/100 d): occupancy · survival-weighted; surviving fraction as a separate 1xN figure; drifting-fraction curve with interannual band |
| 031 | Beaching sweep | Per-run (regime + radius + season) | Five statistics across members with interannual bars · where-stranded maps across members |

## Cross-cutting choices

**Per-regime: overlay vs per-run.** Overlay regimes when the comparison
*is* the point (lines, distance-vs-time). Render per-regime when mixing
would muddle the spatial signal — density and mean-age maps fall here,
because cross-panel reading invites artefacts from regime physics rather
than from the spatial field.

**Release season is the primary temporal axis (025–031).** Every store
consumer takes one `season` parameter and pools all four release years.
The month set is a dict literal in a code cell (papermill injects only
primitives), and the season tag goes in every figure and export filename.

| season | DJF | MAM | JJA | SON | ALL |
|---|---|---|---|---|---|
| months | 12, 1, 2 | 3, 4, 5 | 6, 7, 8 | 9, 10, 11 | 1–12 |

The month comes from `release_doy` through
`pd.to_datetime(f"{year}{doy}", format="%Y%j").dt.month` with the
partition's own year, so it is leap-correct; DJF pairs a year's December
with its own January and February, not the next year's. `020`–`023` still
scope by **quarter**, derived the same way via `.dt.quarter` (the
integer-floor `(doy-1)//90+1` is wrong by ≥1 quarter near boundaries
because Q1–Q4 month lengths differ). Quarter pools across years and year
pools across quarters; crossing them everywhere multiplies panel counts
beyond useful inspection density.

**Interannual spread is the error bar, not a product.** Where a figure
shows one scalar per horizon or per member, the 2016–2019 min–max of the
same reduction — grouped by the release year parsed from the partition
filename — is drawn as a shaded band (029, 030), a lo–hi bar (031) or a
CSV column (028). Drawn as an explicit range rather than a symmetric
`yerr`, because a pooled value need not lie inside it: `beach_hexes` pools
as a union of hexes and exceeds every single year's count. Maps stay
pooled.

**Relative density, not counts (025/026 family).** `percentage` =
100 · n_obs(hex) / Σ n_obs over the selected release set (regime, season,
years, and the age horizon or facet the panel covers) — the share of
particle-time. `dilution` = that share per km² of hex water area
(`water_area_m2` from the 024a key), so a hex that is mostly land is not
credited with a whole hex of dilution volume. Percentage is linear,
dilution log. Absolute counts are dropped from the headline maps — a count
is only readable against a release-set size the reader does not have.
Hexes with `water_area_m2 == 0` carry a percentage but no dilution. The
percentage scale saturates at the 99th percentile with an over-range
colorbar arrow: occupancy is concentrated enough that a few hexes sit two
decades above the bulk, and scaling to the maximum paints every other hex
the same dark colour; the dilution floor at the 1st percentile is the same
argument at the other end.

**Sentinel rows.** `release_hex` or `target_hex` of -1 (land-seeded
particles, obs that left the BSH domain) carry no geometry and are dropped
from numerator and denominator. Each notebook prints the excluded share
(~34.9 % of raw `n_obs` in production — land-seeded releases dominate it).
The land-seeded mask itself is the same `(diff lon at obs=0 == 0) & (diff
lat at obs=0 == 0)` test in every viz; rationale and edge cases in
[distance_calculation.md](distance_calculation.md).

**3 m isobath overlay (025/026 family).** Fucus no longer grows below
~3 m, so every map panel carries the `H0 = 3 m` contour as the viable-band
marker, drawn from the BSH static H0 grids
(`data/bsh_hbmnoku_static/static_file_{fine,coarse}/H0_file_*.nc`), not
from the key's `mean_depth_m` (a 6 km hex mean cannot resolve a 3 m
contour). H0 is NaN over land (see [h0_semantics.md](h0_semantics.md)), so
`contour(levels=[3.0])` masks land for free. Fine grid takes precedence
where it covers; the coarse grid is blanked inside the fine bbox.
Override: magenta at `linewidths=0.3` — anti-viridis, hairline so it does
not read as a second coastline.

**Print-ready figures (025–031).** Every saved figure is one full
text-width (180 mm) figure at 300 dpi: `fig.set_size_inches(180/25.4, h)`
+ `savefig(dpi=300)` ⇒ 2125 px wide. Deliberate break from the AGENTS.md
"no figsize/dpi" rule: these panels are the manuscript figures.
`bbox_inches="tight"` is deliberately not used — it would crop the width
away from 180 mm. The page geometry is module constants
(`FIGURE_WIDTH_IN` / `FIGURE_DPI`) in every consumer, not papermill
parameters: the manuscript's text width is not a per-run knob, and the
inline display size in an executed notebook is irrelevant to the PNG.
Nothing mutates `mpl.rcParams`.

**Panel layout and colorbar placement.** Grids are `layout="constrained"`,
each map's aspect fixed by its extent, and the panel height derived from
page width ÷ ncols ÷ extent aspect plus a fixed decoration allowance — no
draw-measure-rescale loop. Per-panel colorbars are inset axes whose
axes-fraction extent equals the map's, so they track the map at any size
or DPI; they stay per-panel even under a shared norm so columns keep
identical geometry. The two families place them differently, for the same
reason read on different axes:

- **025/026 family** — `ax.inset_axes([0, -0.09, 1, 0.035])`, a horizontal
  bar *under* the map at width 1.0. These grids are up to four columns
  wide, and a side colorbar plus its tick labels would take a quarter of
  the page width per column.
- **027/029/030/031** — `ax.inset_axes([1.02, 0, 0.035, 1.0])`, a vertical
  bar *right* of the map at height 1.0. These rows are short and their
  extent is cropped to the Baltic proper (9–30.7 °E), so the maps are tall
  and narrow: a bar under the map would be far wider than the values need
  and would eat the row height the tall maps want.

Facet grids split across figures past ~12 panels rather than shrinking the
maps; 027/029/030/031 set `hex_seam_lw = 1.1 * 72 / FIGURE_DPI`, and 031's
sweep panels put member tags on the y axis so the long tags read
horizontally.

**Histogram dim layout: keep `obs` as a dim** in 023's `xhist` outputs.
Every scoped plot is then a `.sel`/`.sum`/`.groupby` on the lazy array,
computed once in a shared dask pass. Mean-age reuses the same lazy object
via an obs-weighted variant.

**Coastline overlay** in the hex-map notebooks is Natural Earth 10 m via
cartopy's shapereader, clipped once per extent, drawn black at
`linewidth=0.5` (black reads against the dark low-density band of the
default colormap, which is where registration matters most). Natural Earth
is preferred to the BSH wet-cell coastline because the audience reads
geographic context, not model-grid registration. The cartopy notebooks
(020/022/023) use `add_feature(cfeature.COASTLINE)` per the AGENTS.md
default.

## Notebook 025 — Hex relative-density maps (special case)

Plots on a plain (non-cartopy) lon/lat axis: hex polygons, Natural Earth
coastline, HELCOM subbasin overlay and the H0 isobath are all already in
EPSG:4326; cartopy reprojection would introduce sub-pixel mismatches that
double-paint hex edges. Three products over one (regime, season) with all
four release years pooled: Baltic percentage + dilution, the same clipped
to the German Bight viewport (normalisation stays basin-wide, so the
percentages still read as shares of the whole basin's particle-time), and
a per-origin-subbasin percentage grid where each panel is normalised over
its own releases (colour not comparable between panels). Surviving styling
overrides: hex `edgecolor="face"` / `linewidth=0.4` (the default black
hairline stroke dominates at Baltic zoom, and matching edge to fill makes
the grid read as a continuous field), black coastline `linewidth=0.5`, and
the magenta isobath. The `cmap` parameter is gone — the default colormap
carries the linear percentage and log dilution scales.

## Notebook 026 — Time-horizon relative-density maps (special case)

025's hex rendering re-keyed by the counts store's `age_bin` axis over the
full hex-key domain (the BSH model bbox, North Sea included, so the
geometry rather than a hardcoded box sets the range). A horizon `T` is
**cumulative**: every `age_bin` whose window starts before `T`, i.e.
elapsed time in `[0, T)`, so panels nest and `T` reads as "by 100 days".
Each horizon is normalised to its own 100 %, so the panels show how the
cloud's shape fills in; the concentration trend is read off the colorbars
and the closing quantile-function plot (each hex's share, log y).
Percentage and dilution go to separate figures so each map stays readable
at 180 mm. Horizons must be whole multiples of the store's `age_bin_days`
— a tighter snapshot means rebuilding 024, not a different plot.

## Notebook 026a — per origin subbasin (special case)

026 emitted once per origin (HELCOM) subbasin, filtered on the **release**
side via the key's `helcom_subbasin`, while target hexes still range over
the whole domain — it answers "where do propagules *from this subbasin*
go?". `_outside` (-1) and land-seeded releases are reported, not mapped.
Each origin is normalised and colour-scaled over its own releases
(cross-origin totals differ by orders of magnitude, so a shared scale
would wash out the smaller sources); colour is not comparable between
figures.

## Notebook 026b — per origin subbasin and year (special case)

026a split by release year (parsed from each counts-partition filename),
exposing the interannual variability 026a hides by pooling. Every panel of
an origin divides by that origin's all-year particle-time at the same
horizon, and one norm spans the origin's years and horizons, so a faint
year is a genuinely weaker year and the four years sum to 100 % at a given
horizon (in production each carries ~25 % — release counts are equal by
construction, so the interannual signal is spatial). GeoJSON exported at
the longest horizon only, where the release year rather than the horizon
is the product's axis.

## Exports (025 / 026 / 026a / 026b)

GeoJSON (EPSG:4326, hex geometry) under `output_root/Exports/<nb>/`, named
`<nb>_<regime>_<season>[_T<horizon>d][_<origin>][_<year>].geojson`,
columns `hex_id, n_obs, percentage, dilution` plus `origin_subbasin` /
`release_year` where faceted. GeoJSON over GeoTIFF because the hex
tessellation is vector.

## Notebook 027 — Hex distance-quantile maps (special case)

A **distance-store consumer**: reads the per-source-hex distance histogram
built by 024b (`(release_hex, release_doy, distance_bin)→n_traj`) plus the
024a key for geometry. Pools the run's season across years by summing
histograms, then derives each quantile from the pooled cumulative count —
one panel per quantile level in a 1xN row. Like 025 it draws hex polygons
+ Natural Earth coastline on a plain EPSG:4326 axis, with the same hex
`edgecolor="face"` / black coastline `linewidth=0.5` registration
overrides.

It does **not** inherit 025's density colour scaling: distance is a linear
physical quantity, so each quantile map uses the default colormap with a
matplotlib auto-ranged norm and a colorbar labelled in km. The one
explicit label override is `legend=True` +
`legend_kwds={"label": "final displacement (km)"}`: geopandas derives no
unit-aware axis label from a column (unlike xarray reading
`long_name`/`units`), so the bare column name `distance_km` would reach
the reader without units. Source hexes below `min_traj_per_hex` are
dropped — a quantile from a handful of particles is noise. The quantile
GeoJSON is exported as `Exports/027/027_<regime>_<season>.geojson`.

## Notebook 028 — Subbasin connectivity matrices (special case)

A connectivity-store consumer reading **either** store, selected by the
`member` parameter: `""` reads the 024c partitions
(`(origin_subbasin, target_subbasin, release_doy, age_bin) → n_obs`) and
plots `n_obs`; a rate-model tag reads the 024g survconn partitions and
plots the survival-weighted `w_obs`. It pools the release years, restricts
to one release `season`, drops `-1` (unnamed/outside) rows and reports the
dropped fraction. No hex geometry, no cartopy — every panel is `imshow` on
a plain subbasin × subbasin grid with names as tick labels set at plot
time.

Three views, each its own figure: the raw residence matrix (log, pooled
over every age bin), the emission fraction `x / x.sum(axis=1)` per
cumulative age horizon (half-open, ages `< T` days), and the interannual
mean of that fraction over the four release years (min/max go to CSV).
Row-normalising is the point of view 2: raw rows differ by orders of
magnitude with release count, so only the normalised rows compare origins.

Overrides: 180 mm × 300 dpi as above; one figure per horizon, not a panel
grid (17 HELCOM names are illegible once four panels share a page); run
context in `fig.suptitle`, horizon in `ax.set_title` (long y tick labels
push the axes box right and an axes title overruns the page);
`norm=LogNorm()` (the self-retention diagonal dominates by orders of
magnitude, so a linear scale washes out every off-diagonal; zeros masked
to NaN so empty cells render blank; the emission-fraction norm is floored
four decades below the peak, as in 029). Default colormap — the readable
object is a small labelled grid, not a wide-range map. Cell annotation
with mean and min–max range is gated at `ANNOTATE_MAX_N = 12`; the axis is
17 wide, so it is normally off and the range is read from
`Exports/028/connectivity_<regime>_<season>_<member>_T<h>d_emission_fraction_by_year.csv`.

## Notebook 029 — Beaching maps (special case)

A **beaching-store consumer**: reads the store built by 024e
(`(release_hex, release_doy, beach_hex, beach_age_bin, disp_bin, …) →
weight`) plus the 024a key. The `member` parameter is the opaque
rate-model tag in the store filename
(`HexAgg_beaching_r<radius>m_<regime>_<year>_mMM_<member>.parquet`); it
selects the partitions, tags every figure filename, and is never parsed.
`season` selects which monthly partitions are pooled. Views on one plain
EPSG:4326 axis, reusing 025's hex registration but computing
`hex_seam_lw` from the figure DPI rather than pinning a point value: the
seam is a fixed-pixel anti-aliasing artefact, not a fraction of hex width,
so pinning the stroke to ~1 px keeps the seam closed while letting the
resulting hex dilation shrink as resolution rises.

- **Where-stranded density** — beached weight (expected particles) per
  stranding hex, summed over the store's `disp_bin` (and `shore_type`
  where the reducer emits it) axes, on a `LogNorm` floored four decades
  below the peak (weighted deposition's sparse tail carries fractional
  weight). `shore_type` is deliberately **not shown at all** — neither
  faceted nor broken out in the summary. `shore_type` is a threshold label
  on the HELCOM BRISK + CLMS flat fraction `ff` at the stranding hour
  (`flat` if `ff >= 50%` else `wall`); `trap_flat`/`trap_wall` are reducer
  parameters, named in the member tag only when they differ from 1 — see
  [beaching.md](beaching.md).
- **Beached fraction per source hex** — a ratio in [0, 1], so the default
  linear norm, no colour override to defend.
- **Age-horizon where-stranded maps** — the same log-density map at
  cumulative age cut-offs (`beach_age_bin < T // age_bin_days`), a shared
  `LogNorm` across horizons so the fill-in over time is legible; this
  mirrors 026's horizon logic on the beaching axis.
- **Beached-fraction curve** — cumulative beached fraction against age,
  with the 2016–2019 min–max as a shaded band.
- **Travel distance at stranding** — stranded weight over the store's
  `disp_bin` axis (crow-flies displacement from release at the deposit
  step), pooled over sources, with the weight-weighted median printed.
  Drawn with `drawstyle="steps-mid"`: the axis is a binned quantity, so a
  step reads as the histogram it is, where a categorical bar plot would
  label hundreds of bins. Short-travel strandings stay in — under a
  threshold rate model a strong-wave stranding of freshly released
  material near home is signal, so `disp_bin` is a diagnostic axis and
  never a mask.

## Notebook 030 — Survival heatmaps (special case)

A **survival-occupancy consumer**: reads the store built by 024f
(`(release_doy, age_bin, target_hex) → occ, surv`) plus the 024a key,
selected by the same opaque `member` tag and the same `season` pooling as
029. A horizon `T` names the snapshot bin that **ends** at `T` — ages in
`[T − age_bin_days, T)`, index `T/age_bin_days − 1` — so `T` reads as "at
`T` days old", and the panel titles carry the bin's age span. (026's
horizons are cumulative `[0, T)` instead; this store is a per-age
snapshot, so there is nothing to accumulate.) Per horizon, a 3-column row — plain occupancy, survival-weighted
occupancy (beaching removed), and the surviving fraction `surv/occ` —
reusing 025's hex registration. The two density columns **share one
`LogNorm`** (029's four-decades-floored norm, since survival weights have
a long sub-1 tail) so the age-thinning and the beaching-removal read on
the same scale; the fraction column uses a **fixed linear 0–1** scale so
it is comparable across horizons. A drifting-fraction curve against age,
with the 2016–2019 min–max band, closes the notebook.

## Notebook 031 — Beaching parameter sweep (special case)

Compares the **members** of the 024e store — one rate-model setting each,
named by the reducer's opaque tag (`members_csv`), drawn in the given
order with a display label per tag (`labels_csv`, default = the tags) and
an optional `baseline_member`, for one `season`. Each member is a full
`(year, month)` set. It reports a **range** rather than a number, because
in the Baltic the beaching scheme can dominate the result. A statistics
grid plus a map row:

- **Five statistics across members** — beached fraction, Gini
  concentration of the stranded weight, number of stranding hexes,
  weight-weighted median stranding age, weight-weighted median travel
  distance — one panel each, members on the **y axis** in the given order
  so the long tags read horizontally off one shared column, each pooled
  marker carrying its interannual lo–hi bar. Members differ in several
  parameters at once, so a numeric axis would invite reading a slope where
  there is only a list. Separate panels rather than twin-y axes: the
  quantities share no units and a crossing point would mean nothing.
- **Where-stranded maps across members**, one column per member titled by
  its label, on a **shared `LogNorm`** so the difference read off the row
  is pattern, not scale — the same reasoning as 026's horizon row.

The Gini and travel-distance panels are the load-bearing ones: beached
totals can coincide between members whose stranding patterns differ, so
concentration and how far the stranded material got are what distinguish a
wave-selective member from a residence-driven one.

## Cross-references

- [beaching.md](beaching.md) — 029/031's store, model, and limitations.
- [survival_occupancy.md](survival_occupancy.md) — 030's store and the
  survival-weighting behind it.
- [seeding.md](seeding.md) — release-set semantics every viz reads.
- [distance_calculation.md](distance_calculation.md) — 022's metric.
- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — the
  stores 025/026/027 read, and 028's 024c/024g connectivity section.
- [h0_semantics.md](h0_semantics.md) — the BSH H0 grids behind the 3 m
  isobath overlay.
- [job_scripts.md](job_scripts.md) — how each of these is submitted.
- [../AGENTS.md](../AGENTS.md) — the styling rules this doc defends
  exceptions to.
