# Visualisations: per-plot rationale

Why each notebook in `020`–`036` shows what it shows, and where
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
| 025 | Hex relative-density maps | Per-run (regime + season + radius + member) | Baltic percentage+dilution · German-waters zoom · per origin subbasin |
| 026 | Hex time-horizon relative-density maps | Per-run (regime + season + radius + member) | Cumulative elapsed-time horizons (30/60/90/180 d), `release_years` pooled · mean particle age per hex over the window · area holding 90 % of the cloud at every age-bin edge (CSV) |
| 026a | 026 per origin subbasin | Per-run (regime + season + radius + member) | One percentage and one dilution figure per origin (HELCOM) subbasin |
| 026b | 026a per origin subbasin and year | Per-run (regime + season + radius + member) | One percentage and one dilution figure per (origin subbasin, release year) |
| 027 | Hex distance-quantile maps | Per-run (regime + radius + season) | One 1xN panel row, one panel per quantile (0.1/0.5/0.9), the season's releases pooled across `release_years` |
| 028 | Subbasin connectivity matrices | Per-run (regime + radius + season + member) | Raw residence over the 180 d window · emission fraction per cumulative horizon (30/60/90/180 d) · interannual range |
| 029 | Beaching maps | Per-run (regime + radius + member + season) | Where-stranded with the release sites overlaid · source/sink substrate contrast · stranded share per source hex · cumulative horizons (30/60/90/180 d) · afloat by origin subbasin (table + CSV) · stranded-share curve with interannual band · displacement at stranding |
| 030 | Survival heatmaps | Per-run (regime + radius + member + season) | Per snapshot bin ending at 30/60/90/180 d: occupancy · survival-weighted (% of released); surviving fraction as a separate 1xN figure; drifting-fraction curve with interannual band |
| 031 | Beaching sweep | Per-run (regime + radius + season) | Five statistics across members with interannual bars · where-stranded maps across members |
| 032 | None — selects, numbers and downscales other notebooks' PNGs | Selection (primary + comparison regime) | Tables A–G · one curated reading order over 025–031 |
| 033 | Hex shoreline-type map + stranding-vs-substrate panels | Per-run (regime + radius + season + member pair) | Flat fraction per hex · coastline length per hex · stranded weight vs. flat fraction for the production and the substrate-inert member |
| 034 | Bottom-slowdown comparison | Two regimes (`bottom` vs `bottom_vs0p8`) + one year + radius | Area / √area holding 90 % of particle-time at 30/60/90/180 d · map pair at 90 d |
| 035 | Area time series | Overlaid (one line per regime, each at its 026 weighting tag) + season | Area / √area holding 90 % of particle-time vs. cumulative horizon, every age-bin edge, pooled line + per-year min–max band |
| 036 | Seeding map | Regime-free (one radius) | The 872 release cells by HELCOM subbasin · land-seeded share per cell |

## Cross-cutting choices

**Per-regime: overlay vs per-run.** Overlay regimes when the comparison
*is* the point (lines, distance-vs-time). Render per-regime when mixing
would muddle the spatial signal — density and mean-age maps fall here,
because cross-panel reading invites artefacts from regime physics rather
than from the spatial field.

**Release season is the primary temporal axis (025–031).** Every store
consumer takes one `season` parameter and pools the release years of its
`release_years` parameter. The month set is a dict literal in a code cell
(papermill injects only primitives), and the season tag goes in every figure
and export filename.

**Release years 2016–2024.** Every consumer that finds partitions by
filename (025–033) takes `release_years = "2016-2024"`, an inclusive
`"first-last"` range or a single year (`"2019"`), parsed in a code cell and
applied to the filenames before any partition is opened; it prints the range
and the years on disk it left out. Pooled values and per-year bands therefore
cover exactly those years. The partial 2025 (January–May releases only) stays
in the stores and is not read: its missing seasons would pull every
per-year band down, and a full 2025 needs 2026 forcing. The driver and job
scripts pass no year argument, so the notebook default is the study's year
set. 026b keeps a separate `mapped_years_csv` for which years get panels,
since each origin is normalised over all of `release_years`; 034 compares one
fixed `year`; 035 reads 026's CSVs, whose per-year rows already end at the
range. 032 lists every year on disk in its coverage table and pools only
`release_years` in tables A–C (D–G read the consumers' exports).

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
shows one scalar per horizon or per member, the per-release-year min–max of the
same reduction — grouped by the release year parsed from the partition
filename — is drawn as a shaded band (029, 030, 035), a lo–hi bar (031), a
CSV column (028) or per-year CSV rows (026 areas). Drawn as an explicit
range rather than a symmetric
`yerr`, because a pooled value need not lie inside it: `beach_hexes` pools
as a union of hexes and exceeds every single year's count. Maps stay
pooled.

**Weighted or unweighted (025/026 family).** One `member` parameter picks
the store, and the notebook calls the value it reads `weight` from then on:
`""` reads the 024 counts store (`n_obs`, every hourly sample counts 1), a
rate-model tag reads that member's 024f survocc store (`surv`, each sample
weighted by its un-beached fraction `exp(−A)`, see
[survival_occupancy.md](survival_occupancy.md)). Both are cut to the same
age window `[0, max_age_days)` (180 d): the counts store is uncapped (220 d),
so without the cut a weighted and an unweighted map would cover different
ages. A survocc read asserts that every pooled partition's
`df.attrs["max_age_days"]` agrees and covers the consumer's window. Origin
subbasin is the key's `helcom_subbasin` of `release_hex` for counts and the
store's own `origin_subbasin` for survocc — the same assignment, `-1` for
a release hex outside every named polygon. Every figure and export stem
carries the member tag or the literal `unweighted`. Considered a separate
weighted notebook family; rejected, the maps, norms and exports are
identical code over a different value column.

**One window, recorded in the store.** The 024e/024f/024g partitions carry
their build parameters as `df.attrs` (window `max_age_days`, `age_bin_days`,
the rate-model point; [beaching.md](beaching.md) § Store schema). Every
consumer of them reads the attrs per file *before* concatenating — pandas
drops attrs that differ between the frames it concatenates, and
`release_year`/`release_month`/`git_sha` always do — and asserts that every
pooled partition agrees on `max_age_days`, `age_bin_days` (`disp_bin_km` for
024e) and the rate-model keys, that the store's binning equals the
notebook's parameters, and that `max_age_days` covers the largest horizon.
Provenance keys are never compared. A partition with empty attrs (written
before the window was recorded) raises with a rebuild request. Where members
are compared (029's contrast member, 031, 032, 033) the rate model differs by
design and only the window and binning must agree. Considered trusting the
filenames; rejected, a partial rerun at another window would mix silently.

**Units: % of released drifters (029, 030, 031).** The beaching and survival
maps and curves are percentages of the released water-seeded drifters of the
pooled partitions, so a map, its curve and the 029 afloat table share one
denominator and read without knowing the release-set size. In 029 and 031 the
denominator is the pooled `Σweight` of the 024e store (deposits + residual per
release hex = its released drifters). In 030 it is the exact 024f denominator:
with the `target_hex = -1` sentinel rows kept, a snapshot bin's `Σocc` over all
origins is `N · age_bin_days · 24`, so `100 · surv / Σocc` is the mean share of
released drifters afloat in that hex during the bin
([survival_occupancy.md](survival_occupancy.md)). 029's GeoJSON keeps the
absolute `stranded_weight` beside `pct_released`.

**Log-scale floors.** Every log map of a share — 028's emission fraction,
029/031's stranded share, 030's densities — is floored a fixed **four decades
below its peak**: the `exp(−A)` tail of a weighted quantity runs arbitrarily
far down, and how far depends on the rate-model member, so a percentile floor
would give every member its own dynamic range. The 025/026 family's
**dilution** maps, weighted or not, keep the 1st-percentile floor (below).
Considered the four-decade floor there too; rejected, dilution divides by
water area, so its peak sits in near-land hexes with a sliver of water: on
a survival-weighted surface_stokes JJA test pool (2019–2020) at 180 d the
median hex sits 5.3 decades below the peak (the unweighted maps span ~8), and
a four-decade floor painted more than half the basin the floor colour.

**Relative density, not counts (025/026 family).** `percentage` =
100 · weight(hex) / Σ weight over the selected release set (regime, season,
years, and the age horizon or facet the panel covers) — the share of
particle-time (survival-weighted particle-time for a member). `dilution` =
that share per km² of hex water area
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
(~34.9 % of raw `n_obs` in production — land-seeded releases dominate it;
024f holds water-seeded drifters only, so a survocc read drops just the
`target_hex = -1` rows). Origin `-1` releases stay in the basin-wide maps
and are left out of the per-origin panels.
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

**Print-ready figures (025–031, 033–036).** Every saved figure is one full
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
- **027/029/030/031/033** — `ax.inset_axes([1.02, 0, 0.035, 1.0])`, a vertical
  bar *right* of the map at height 1.0. These rows are short and their
  extent is cropped to the Baltic proper (9–30.7 °E), so the maps are tall
  and narrow: a bar under the map would be far wider than the values need
  and would eat the row height the tall maps want.

Facet grids split across figures past ~12 panels rather than shrinking the
maps; 027/029/030/031/033 set `hex_seam_lw = 1.1 * 72 / FIGURE_DPI`, and 031's
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
double-paint hex edges. Three products over one (regime, season) with the
`release_years` pooled: Baltic percentage + dilution, the same clipped
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
and ≤ `max_age_days` — a tighter snapshot means rebuilding 024, not a
different plot.

**Area holding a share of the cloud.** Per horizon `T` in
`area_horizons_days_csv` — default `""`, every age-bin edge up to
`max_age_days` (10, 20, …, 180 d), so 035 can draw it as a time series;
the maps keep `time_horizons_days_csv` — and share
`q` (`area_quantiles_csv`, default 0.9 only): hexes ranked by density
`weight / water_area_m2` descending, weight accumulated until it reaches
`q`, and their water area summed — the smallest water area holding a share
`q` of the cumulative weight. The hex that crosses `q` counts only the
fraction of its area needed (a 6 km hex is 93.5 km² of water, so whole
hexes would quantise small-T areas); hexes with weight but no water area
rank first at zero area. `length_km = sqrt(area_km2)` is the side of the
equal-area square, a length scale rather than a distance travelled.
`pooled` ranks the summed weights of every year in `release_years`; each year is the
same reduction on its own releases, so the per-year rows give the
interannual range (a pooled value need not lie inside it). A year whose
partitions hold no release in the season (the partial last release year
has no JJA or SON) is dropped and reported. The ranking is over cumulative
particle-time, so the area measures concentration, not spread: where
drifters accumulate in few hexes (coastal trapping), later horizons add
weight there and the area can shrink as `T` grows. Only `q = 0.9` is
reported: under advection-only physics the 0.5 area is dominated by particles
stalling in a few coastal hexes and is hard to interpret. Each scope's
per-`(target_hex, age_bin)` weights are accumulated over `age_bin` once, so a
dense horizon set costs one ranking per horizon, not one regrouping. No
figure here: the CSV feeds 035 (the time series) and 032 (table F, which
filters the CSV to its own `area_quantile`).

**Mean particle age.** One map per run over the whole window: per hex,
`Σ age_centre · weight / Σ weight` with `age_centre = (age_bin + 0.5) ·
age_bin_days`, ages in `[0, max_age_days)`, from the same per-`(age_bin,
target_hex)` reduction as the horizon maps, so a member run is
survival-weighted and an unweighted run is plain residence. Default
colormap on a linear scale over the mapped range; colorbar label `mean
particle age (d)`. **Mask:** a hex is drawn only if its dilution at
`T = max_age_days` reaches the 1st percentile that floors the window's
dilution map (`dilution_norm`) — the hexes that map shows above its floor,
zero-water hexes excluded. At production this drops ~1 % of hexes and
0.1–0.3 % of the weight. Considered an absolute sample-count cut; rejected,
a survival weight has no sample count and the two weightings would need two
rules. The run prints the weight-weighted median of the hex mean ages, the
unweighted median, the range and the basin-wide mean age. Unweighted, the
basin-wide mean is the window midpoint (90 d at 180 d) up to the share that
leaves the domain — a check on the reduction, not a result. Figure
`MeanAge_<regime>_<season>_r<radius>m_<member|unweighted>.png`, export
`026_meanage_<regime>_<season>_<member|unweighted>.geojson` (mapped hexes
only, the window's `weight, percentage, dilution` plus `mean_age_days`).

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
year is a genuinely weaker year and the release years sum to 100 % at a given
horizon. Unweighted, each year carries an equal share (release counts are
equal by construction, so the interannual signal is spatial); weighted, a
year's share also carries its survival. GeoJSON exported at
the longest horizon only, where the release year rather than the horizon
is the product's axis.

## Exports (025 / 026 / 026a / 026b)

GeoJSON (EPSG:4326, hex geometry) under `output_root/Exports/<nb>/`, named
`<nb>_<regime>_<season>_<member|unweighted>[_T<horizon>d][_<origin>][_<year>].geojson`,
columns `hex_id, weight, percentage, dilution` plus `origin_subbasin` /
`release_year` where faceted. GeoJSON over GeoTIFF because the hex
tessellation is vector. 026's mean-age export
`026_meanage_<regime>_<season>_<member|unweighted>.geojson` adds
`mean_age_days` (window-wide, no horizon facet). 026 also writes
`026_areas_<regime>_<season>_<member|unweighted>.csv` with columns
`horizon_days, q, scope, area_km2, length_km`, `scope` = `pooled` or a
release year. Figures follow the same rule —
`<Name>_<regime>_<season>_r<radius>m_<member|unweighted>[_<facet>…].png`:
the member tag follows the run identity and precedes the per-figure facets
(column, origin, year, part).

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
dropped fraction. With the shoreward subbasin rule
([hexbinning_and_connectivity.md](hexbinning_and_connectivity.md)) that
fraction is the genuinely un-Baltic part of the domain plus, for the
unweighted store only, the land-seeded particles counts keeps as
`release_hex == -1`; 024g's population excludes those by construction, so
the weighted matrix drops less than the unweighted one. No hex geometry, no cartopy — every panel is `imshow` on
a plain subbasin × subbasin grid with names as tick labels set at plot
time.

Three views, each its own figure: the raw residence matrix (log, pooled
over the window, ages `< max_age_days`, stem `_T<max_age_days>d`), the
emission fraction `x / x.sum(axis=1)` per cumulative age horizon (half-open,
ages `< T` days, default 30/60/90/180 d), and the interannual mean of that
fraction over the release years (min/max go to CSV). Both stores are cut to
the same window: the 024c store is uncapped (220 d), 024g records its own
window in the attrs, so without the cut the weighted and the unweighted raw
matrix would cover different ages.
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
`season` selects which monthly partitions are pooled. Two files come from the
data twin rather than the store: the release-site geojson behind the overlay
and the shoreline-class table behind the substrate readout, both below. Views
on one plain EPSG:4326 axis, reusing 025's hex registration but computing
`hex_seam_lw` from the figure DPI rather than pinning a point value: the
seam is a fixed-pixel anti-aliasing artefact, not a fraction of hex width,
so pinning the stroke to ~1 px keeps the seam closed while letting the
resulting hex dilation shrink as resolution rises.

- **Where-stranded density** — stranded weight per stranding hex as % of the
  released drifters, summed over the store's `disp_bin` (and `shore_type`
  where the reducer emits it) axes, on a `LogNorm` floored four decades
  below the peak (weighted deposition's sparse tail carries fractional
  weight). `shore_type` is deliberately **not shown at all** — neither
  faceted nor broken out in the summary. It is a threshold label on the
  HELCOM BRISK + CLMS flat fraction `ff` at the stranding hour (`flat` if
  `ff >= 50%` else `wall`), a coarser cut than the rate, which scales
  continuously with `ff` through `trap_wall + (trap_flat − trap_wall)·ff/100`
  — at the production member (`trap_flat = 1`, `trap_wall = 0`) the substrate
  signal is already inside the plotted weight, so a flat/wall split of it
  would re-report the rate model as if it were a result. See
  [beaching.md](beaching.md).
- **Release-site overlay** (this panel only, `show_release_points`) — the 872
  REDLIST cells with `F_vesiculo != 0`
  (`data/helcom_fucus_redlist/fucus_release_points.geojson`, baked by `000`),
  one hollow marker per cell at its representative point. This is the file
  `010` opens and seeds particles inside, so the overlay is the seeding set,
  not a proxy. Without it the panel is a deposition field and the reader has
  to supply the sources from memory; with it the figure is the source→sink
  statement the study is about. **Overrides: `markerfacecolor="none"`,
  `markeredgecolor="magenta"`, `markersize=2.5`, `markeredgewidth=0.4`.**
  Hollow because sources and sinks overlap across the western Baltic and a
  filled marker would hide the hex value the reader is comparing against;
  magenta because the field underneath is a four-decade `LogNorm` of the
  default map and the marker has to read against both a near-black low end
  and a yellow high end (the same anti-viridis argument as the 3 m isobath in
  025/026); 2.5 pt because that is roughly the 10 km REDLIST cell's own size
  at 180 mm page width, so the marker stands for the cell rather than
  inventing a point source. Considered and rejected: the cell polygons drawn
  as outlines (faithful, but 10 km squares are larger than the 6 km hexes and
  mesh together in the Danish straits and the Åland archipelago); the key's
  `fucus_area_m2 > 0` hexes drawn as outlines (1 827 hexes against 872 cells —
  a 10 km cell touches several 6 km hexes, so the footprint is overstated, and
  a continuous magenta band along the whole Swedish and Finnish coast
  outweighs the field it annotates); filled dots (clean, but they read as
  scattered noise rather than as sites, and they occlude). The overlay is
  deliberately **not** repeated on the beached-fraction or age-horizon panels:
  the beached-fraction panel's coloured hexes already *are* the release set,
  and the horizon row is about time, not about where the material came from.
- **Source and sink substrate** — printed, not plotted: the shoreline flat
  fraction `ff` at the release sites against the flat fraction where the
  weight lands, aggregated per hex the way 033 does it (nest seam dropped,
  unattributed coast counted 0.5) and reported for `member` and for a
  substrate-inert `contrast_member`, plus a release-cells-against-stranded-
  weight table per HELCOM subbasin. What it shows: *Fucus* sits on coast that
  is **wall-like relative to the domain** (`fucus_area_m2`-weighted mean `ff`
  0.427 against a domain mean of 0.527; 57 % of release cells on hexes with
  `ff < 0.5`), while under the production member the deposit sits on flat
  coast (`ff` 0.801, 80.6 % of weight on flat hexes) — and the substrate-inert
  member deposits at `ff` 0.389, i.e. **the flat bias is the rate model's trap
  term, not the drift**, which is exactly why both members are printed. The
  subbasin table shows where the two diverge: Western Gotland Basin, Bothnian
  Sea and Åland Sea hold 33.3 % of the release cells and receive 5.6 % of the
  production stranded weight (29.5 % inert), while under production Great
  Belt and Gulf of Finland receive more than they release (10.7 → 16.7 %,
  10.6 → 19.5 %). Numbers from the `ALL` run, release years 2016–2024, 180 d
  window. Reading the per-hex correlation behind this is 033's job; 029 states the two means so the overlaid map has numbers
  attached.
- **Stranded share per source hex** — % of each hex's own released drifters
  stranded within the window, in [0, 100], so the default linear norm, no
  colour override to defend.
- **Age-horizon where-stranded maps** — the same log-density map at
  cumulative age cut-offs (`beach_age_bin < T // age_bin_days`, default
  30/60/90/180 d), a shared `LogNorm` across horizons so the fill-in over time
  is legible; this mirrors 026's horizon logic on the beaching axis. Exported
  per horizon as `Exports/029/029_<regime>_<season>_<member>_T<h>d.geojson`
  with `stranded_weight` (expected particles) and `pct_released`.
- **Afloat by origin subbasin** — printed and exported, not plotted: per
  HELCOM origin subbasin (the key's `helcom_subbasin` of `release_hex`),
  `afloat(T) = Σweight − Σweight[beach_hex ≥ 0 & beach_age_bin < T/age_bin_days]`
  as % of the origin's `Σweight` (its released drifters), at the same
  horizons; pooled from the summed weights plus the per-release-year min–max.
  Origin `-1` (release hex outside every named polygon) is its own row, and
  `all` is the basin-wide value. "Afloat" is everything not stranded: it
  includes drifters that left the BSH domain, and no lifetime or sinking is
  applied — the CSV
  (`Exports/029/afloat_by_origin_<regime>_<season>_<member>.csv`, columns
  `origin_subbasin, origin_name, horizon_days, released_weight, afloat_pct,
  afloat_pct_min, afloat_pct_max, n_years`) says so in its `#` header.
- **Stranded-share curve** — cumulative stranded % of released against age,
  with the per-release-year min–max as a shaded band.
- **Displacement at stranding** — stranded share (% of released) over the store's
  `disp_bin` axis (crow-flies displacement from release at the deposit
  step), pooled over sources, with the weight-weighted median printed.
  Drawn with `drawstyle="steps-mid"`: the axis is a binned quantity, so a
  step reads as the histogram it is, where a categorical bar plot would
  label hundreds of bins. Short-displacement strandings stay in — under a
  threshold rate model a strong-wave stranding of freshly released
  material near home is signal, so `disp_bin` is a diagnostic axis and
  never a mask.

## Notebook 030 — Survival heatmaps (special case)

A **survival-occupancy consumer**: reads the store built by 024f
(`(release_doy, age_bin, origin_subbasin, target_hex) → occ, surv`) plus the
024a key, selected by the same opaque `member` tag and the same `season`
pooling as 029. A horizon `T` (default 30/60/90/180 d) names the snapshot bin
that **ends** at `T` — ages in `[T − age_bin_days, T)`, index
`T/age_bin_days − 1` — so `T` reads as "at `T` days old", and the panel
titles carry the bin's age span. A snapshot, not 026's cumulative `[0, T)`:
being afloat at `T` is an instantaneous state, and a survival map
accumulated over `[0, T)` would be a weighted residence map, which 026
already draws. The two density panels — plain occupancy and
survival-weighted occupancy (beaching removed) — are in % of the released
drifters (`100 · occ / Σocc`, `100 · surv / Σocc` with the bin's exact `Σocc`,
see Units above) and **share one `LogNorm`** floored four decades below the
peak, so the age-thinning and the beaching-removal read on the same scale;
the surviving fraction `surv/occ` is a separate 1×N figure on a **fixed
linear 0–1** scale so it is comparable across horizons. Every panel sums over
`origin_subbasin`; the `target_hex = -1` sentinel drops at the key join.

A drifting-fraction curve `Σsurv / Σocc` against age, with the
per-release-year min–max band, closes the notebook. It **keeps** the
sentinel: the curve is then the mean surviving fraction of every released
drifter, including those that left the BSH domain (they keep the survival
they had on leaving), which is the same population 029's afloat table
counts. Dropping the sentinel would condition on the drifters still in the
domain, a population that shrinks with age. The per-bin table prints `N` and
the share of samples at `target_hex = -1`.

## Notebook 031 — Beaching parameter sweep (special case)

Compares the **members** of the 024e store — one rate-model setting each,
named by the reducer's opaque tag (`members_csv`; the production driver passes
every member on disk), drawn in the given
order with a display label per tag (`labels_csv`, default = the tags) and
an optional `baseline_member`, for one `season`. Each member is a full
`(year, month)` set, and every member must share one window
(`max_age_days`, `age_bin_days`, `disp_bin_km`, asserted from the attrs) —
a member reduced over a longer window has had longer to strand. Stranded
share and maps are in % of released, the maps on a shared `LogNorm` floored
four decades below the peak. It reports a **range** rather than a number, because
in the Baltic the beaching scheme can dominate the result. A statistics
grid plus a map row. Each partition is reduced as it is read — its total
weight (the denominator) and its stranded weight per `(beach_hex,
beach_age_bin, disp_bin)` with the release year — so memory holds the
reduced frames, not the ~10⁸ raw rows per member a year-pooled season has
(measured peak 6.2 GiB for 14 members at season ALL, 2016–2024; the raw
rows would be ~112 GB):

- **Five statistics across members** — stranded % by `max_age_days`, Gini
  concentration of the stranded weight, number of stranding hexes,
  weight-weighted median stranding age, weight-weighted median displacement — one panel each, members on the **y axis** in the given order
  so the long tags read horizontally off one shared column, each pooled
  marker carrying its interannual lo–hi bar. Members differ in several
  parameters at once, so a numeric axis would invite reading a slope where
  there is only a list. Separate panels rather than twin-y axes: the
  quantities share no units and a crossing point would mean nothing.
- **Where-stranded maps across members**, one column per member titled by
  its label, on a **shared `LogNorm`** so the difference read off the row
  is pattern, not scale — the same reasoning as 026's horizon row.

The Gini and displacement panels are the load-bearing ones: beached
totals can coincide between members whose stranding patterns differ, so
concentration and how far the stranded material got are what distinguish a
wave-selective member from a residence-driven one.

## Notebook 032 — Headline report (special case)

[`032_HeadlineReport`](../notebooks/032_HeadlineReport.py) is the only stage
here that draws nothing. It regenerates `output_root/Headline/` — a
`HEADLINE.md` carrying the study's seven headline tables plus a numbered,
downscaled copy of the figures those tables are read against. Submitted by
[`scripts/headline_report_job.sh`](../scripts/headline_report_job.sh), the
last stage of the downstream chain, so the report regenerates whenever new
release years land.

Every number is **recomputed**, never transcribed:

| table | what | read from |
|---|---|---|
| A | stranded % within the window per release season and year, production member | `HexAgg_beaching_*` |
| B | the same pooled value across every rate-model member on disk | `HexAgg_beaching_*` |
| C | `Σw_obs/Σn_obs` per season and cumulative horizon (`connectivity_horizons_csv`, 30/60/90/180 d), plus the top-target comparison at `top_horizon_days` (60 d) | `HexAgg_survconn_*`, `Exports/028/*.csv` |
| D | median over source hexes of each displacement quantile | `Exports/027/*.geojson` |
| E | concentration and viable-band share at `horizon_days` (180 d) | `Exports/026/026_<regime>_<season>_<weighting>_T<h>d.geojson` + the 024a key |
| F | water area (and its √) holding `area_quantile` (0.9) of the cumulative particle-time at the `area_horizons_csv` (30/60/90/180 d) map horizons, pooled [per-year min–max; the pooled value can lie outside it]; CSV rows at other `q` or horizons are filtered out | `Exports/026/026_areas_<regime>_<season>_<weighting>.csv` |
| G | afloat % per origin subbasin per horizon, pooled [per-year min–max], plus the basin-wide row per season | `Exports/029/afloat_by_origin_<regime>_<season>_<member>.csv` |

`<weighting>` is the regime's 025/026 tag: the production member where the
regime has a 024f store (surface_stokes), `unweighted` otherwise (bottom) —
the rule the consumer driver applies, discovered from the partitions on
disk. Table F reads its CSV generically — `scope` splits pooled from
per-year rows, every float column but `q` is a value, every other column
but `q` a row key — so a key column the export gains later becomes rows
without an edit; a requested horizon the CSV lacks is a gap. The window every number refers to comes from the 024e/024g attrs
(asserted equal across members and between the two stores), not from a
parameter.

Tables A–C go to the parquet partitions rather than to `Exports/`, because
the quantities they need are not exported: a pooled multi-year beached
fraction is a ratio of summed weights (not a mean of yearly ratios), and
`Σw/Σn` needs both weights on the same rows, which only the 024g store
carries — 028 writes one value column per CSV. D and E read `Exports/`
instead, so the table and the figure beside it are the same reduction; 027's
`min_traj_per_hex` cut and 026's cumulative horizon would otherwise be
reimplemented and drift.

**Everything that grows is discovered.** Release years, regimes and
rate-model member tags are parsed out of the partition filenames, so a store
growing from few release years to many, or from one member to a full trap
sweep, changes the report's content and never this notebook. A coverage
table states the years each store holds per regime — the regimes are not in
step (`surface` stops at the aggregates; beaching and survival exist only
where a 024d sidecar was built), and a headline that pools silently over an
uneven year set is worse than no headline.

**Figure selection is a rule, not a list.** Each entry names the producing
notebook and a filename *stem* built from the run's own parameters
(`{regime}`, `{other_regime}`, `{season}`, `{pooled}`, `{radius}`,
`{member}`, `{tag}`/`{other_tag}` — the regime's 025/026 weighting —
`{horizon}` = `top_horizon_days`); resolution takes the exact `<stem>.png`
only. Considered also accepting the first `<stem>_*.png`; rejected, a member
tag can extend another (`step_wc0p1_ts3_tcinf` →
`step_wc0p1_ts3_tcinf_tf1_tw0`), so a missing figure would silently resolve
to a different member's. The map entries use the `figure_season` parameter
and the curve/sweep entries the year-pooled season, so the manifest is one
reading order rather than twelve hardcoded names. Captions take their
horizons and windows from the parameters, the store attrs and the 026
exports on disk. An entry that resolves to nothing is recorded in a "Gaps"
section of `HEADLINE.md` rather than raising: a partially-run consumer
matrix should still yield a report that says what is missing. An entry read
against another (`pairs_with`) names its figure number, or says it is
absent.

**Stale is missing.** Every export a table reads and every figure entry
carries its sources — `(store kind, regime, member)` — and is a gap if its
mtime is older than the newest `release_years` partition of those stores or
than the 024a key. Otherwise a store rebuilt at a new window would leave the
old consumer products in place under the same names, captioned with the new
window. Sources: 024f at the production member for a weighted 025/026
product, else the 024 counts; 024b for 027; 024c / 024g for unweighted /
weighted 028; 024e for 029 and (every member) 031; 024f for 030. 035
draws 026's area CSVs, so its sources are those CSVs, each checked against
its own stores — a 035 figure drawn from a stale CSV is stale too. 036 has
no store source and is checked against the 024a key alone (its subbasin
assignment; the trajectory zarrs behind its land share are not tracked —
the seeding does not change between runs). Considered recording input build ids in
every product; rejected, mtime needs no consumer change (copying outputs
without preserving mtimes defeats it).

Copies are downscaled to `figure_width_px` (1400) with PIL; the 180 mm /
300 dpi originals stay under `Figures/` as the manuscript versions, so the
print-ready rule above is untouched — 032 adds no styling override because
it opens no figure. Stale `NN_*.png` are deleted before the copy, since the
numbers are positions in the reading order, not identities.

The narrative lines between the tables are computed from the tables
themselves (seasonal spread in points, the member range, the regime
distance ratio), never written as prose that assumes a result.

## Notebook 033 — Shoreline type map (special case)

The substrate axis of the beaching rate model, drawn once as a map and then
tested against the stranding it drives. Reads the data twin's
`shoreclass_bsh_coastline/bsh_coastline_k2_flatfraction.parquet` (65 876
BSH-coastline sub-segment midpoints in EPSG:3035, each with a `flat_fraction`,
a `seg_len_m` and an `on_seam` flag), the 024a key, and the 024e beaching
store — parquet only, no Dask.

**The nest seam is dropped first.** That table traces the outline of the BSH
model's wet cells, and the outline closes along the edge of the fine nest as
well as along real land: the model's water simply stops there. `on_seam` marks
those sub-segments — 1 100 of them, 472 km, 1.6 % of the modelled coastline,
all but 9 km being the nest's northern edge at 56.4500 N in two dead-straight
open-water runs (North Sea 6.2–8.2 E, Kattegat 10.9–12.7 E south of Anholt).
Left in, they drew a row of `flat_fraction = 0.5` hexes across open water on
both maps, inflated the coastline-length map, and put 11 phantom hexes in the
fit set. 024d keeps the same rows and should: a particle reaching the nest edge
does meet a wall of land-flagged cells there. The beaching products are
unaffected either way — the hexes that hold seam coastline and nothing else
carry exactly 0 stranded weight under the production member and 0.0003 % of it
under the substrate-inert contrast.

**Per hex: the `seg_len_m`-weighted mean flat fraction**, plus the total
coastline length. The midpoints are reprojected 3035 → 4326 (the same `pyproj`
transform 024d uses for the reverse direction) and labelled through the key's
`HexProj`, so hex assignment is the one every `024x` stage uses. A sub-segment
with no class attribution counts **0.5**, the convention 024d applies
([beaching.md](beaching.md)); the length share that is unattributed is carried
as its own column, annotated in the map title and printed. **A hex with no
coastline at all gets no row** — it renders absent rather than as a mid-scale
"mixed" value, which a fill would be indistinguishable from.

Two members are reported, one relationship figure each, and every figure and
printed block names which is which: the
production `step_wc0p1_ts3_tcinf_tf1_tw0` (`trap_flat = 1`, `trap_wall = 0`, so
the realised factor is `trap = ff/100` and **substrate drives the trapping
term**) against `contrast_member = step_wc0p1_ts3_tcinf` (`trap_flat =
trap_wall = 1`, **substrate inert**). Read together they separate "stranding is
on flat coast because the model puts it there" from "stranding is on flat coast
because that is where the drift goes" — the latter is the only thing the inert
member can show.

- **Flat fraction per hex** — the headline map, 029's hex registration and
  Baltic extent verbatim. **One styling override: `cmap="coolwarm"` with a
  fixed `Normalize(0, 1)`.** Justified because this is not a magnitude but a
  bounded two-ended axis (0 = fully walled, 1 = fully flat) with a
  *distinguished midpoint*: 0.5 is exactly what the pipeline assigns coastline
  it has no class for. The default sequential map has no distinguished
  midpoint, so "half flat" and "unattributed" would land in an arbitrary
  mid-ramp colour and the two ends would not read as opposites. The norm is
  pinned rather than auto-ranged so the white stays at 0.5 across runs.
- **Coastline length per hex** — the denominator of the normalised stranded
  weight, and the map that says how much of the flat-fraction map rests on a
  sliver of coast. Default colormap, default linear norm.
- **Stranded weight vs. shore type**, one figure per member, laid out as a
  `subplot_mosaic` — two half-width panels over one full-width panel, because
  the subbasin panel carries 17 long HELCOM names on its y axis and is
  illegible at a third of the page (028's reasoning), while the other two read
  fine at half width. Top left: the scatter of stranded weight per km of coast
  against hex flat fraction with its Spearman ρ (log y; the per-km value is the
  fair comparison, since a hex with 40 km of coast catches more than one with
  4 km whatever it is made of, and the raw-weight ρ is printed beside it; the
  count of zero-weight hexes that a log axis cannot show goes in the panel
  title). Top right: the share of stranded weight and the share of coastline
  length over flat-fraction bins, as grouped bars — **if stranding merely
  followed the coastline's own composition the two series would coincide**, so
  the gap between them is the whole answer. Bottom: the per-subbasin ρ, names
  on the y axis (031's reasoning: long labels read horizontally off one shared
  column).

  The fit set is the coastal hexes holding at least `min_coast_len_m` (1 km)
  of coastline whose coastline is at least `min_baltic_fraction` (half) inside
  the Baltic. The length cut is because a length-weighted mean over a sliver of
  coast is noise and its per-km weight divides by almost nothing. The Baltic
  cut is the sentinel rule adapted: the key tiles the whole BSH domain, so
  North Sea and Skagerrak coast the Baltic releases never reach would let "no
  Fucus released within reach" masquerade as a substrate signal. It uses the
  shoreline table's own `in_baltic` flag rather than the key's
  `helcom_subbasin >= 0`, because the HELCOM polygons leave 859 genuine Baltic
  coastal hexes unassigned — a coastal hex centred on land or inside an
  archipelago falls outside every polygon — and dropping those would discard
  about half the stranded weight. The per-subbasin panel can only show named
  subbasins, so the share of fit-set weight it leaves out is printed alongside
  the share the fit set itself drops. The maps are **not** cut: they show every
  hex that carries coastline.

The `shore_type` flat/wall split of stranded weight is printed, not plotted,
per subbasin and pooled over the Swedish-coast subbasins
(`swedish_subbasins_csv`: Western Gotland Basin, Eastern Gotland Basin,
Northern Baltic Proper, Bothnian Sea). `shore_type` is a *threshold* on `ff`
while the trap factor is *linear* in it, so under the production member
wall-labelled sites still strand at `ff/100` of the flat rate and only `ff = 0`
is fully non-trapping; under the contrast member the label moves no weight at
all and the split is the coastline's own composition as the drift samples it.

Both members' partitions are checked as 029/031 do (attrs per file, one
window and rate-model point per member, one window across the two).

Export: `Exports/033/033_<regime>_<season>_<member>.geojson`, columns `hex_id,
flat_fraction, coast_len_m, stranded_weight` plus the 024a hex geometry
(EPSG:4326), matching the 025–029 export convention.

## Notebook 034 — Bottom slowdown (special case)

The one sensitivity illustration for bottom drift: the `bottom` regime against
`bottom_vs0p8`, the same 010 runs with U and V scaled by 0.8
([job_scripts.md](job_scripts.md), velocity-scaled variant regimes). The
variant exists for one release year and one release per month, so **both
sides are cut to those release doys** (`release_doys_csv`); a doy missing from
either counts partition raises rather than shrinking one side. Parquet only.

- **Footprint table and CSV** — the area holding `q` (default 0.9 only, as
  in 026: the 0.5 area is dominated by particles stalling in a few coastal
  hexes) of the
  particle-time accumulated by `T` (cumulative `[0, T)`, as in 026): hexes
  ranked by `n_obs / water_area_m2` descending, zero-water hexes first at zero
  area, the crossing hex taken fractionally; `length_km = √area_km2`.
  `Exports/034/034_areas_<year>.csv` carries `regime, horizon_days, q,
  area_km2, length_km` and the variant-over-reference area and length ratios.
  Implemented notebook-locally rather than read from 026's areas export,
  because 026 pools a whole season over every year, not a doy subset.
- **Length vs. horizon** — one panel per `q`, one line per regime (`hue`), so
  no line-style override is needed to tell four series apart.
- **Map pair** at `map_horizon_days` — 026's percentage and dilution rendering
  (horizontal colorbar under each map), one figure per view, the two regimes
  side by side under **one shared norm** so the same colour is the same share
  in both; unlike 026's per-panel norms, the comparison is the point here.

Figures at 180 mm × 300 dpi as above.

## Notebook 035 — Area time series (special case)

[`035_AreaTimeSeries`](../notebooks/035_AreaTimeSeries.py) draws 026's area
exports `Exports/026/026_areas_<regime>_<season>_<tag>.csv` for several
regimes on shared axes — R1's "how much water holds the cloud by T". No
parquet; it runs after 026 for every regime it draws.

| parameter | default | |
|---|---|---|
| `regimes_csv`, `tags_csv` | `surface_stokes,bottom`, `<production member>,unweighted` | one 026 weighting tag per regime, position by position (D3) |
| `season` | `ALL` | one release season per figure |
| `area_quantile` | 0.9 | must be a `q` the CSV holds, else it raises |
| `hex_radius` | 6000 | labels the figure only; the CSV names carry no radius |

- **Two panels**, area (km²) and √area (km), against the cumulative horizon
  `T` (ages in `[0, T)`, every age-bin edge 026 exported), y from zero so the
  regimes compare as ratios.
- **One colour per regime** from the default cycle, indexed by its position
  in `regimes_csv` so a regime keeps its colour when another is missing; the
  pooled value is a line with markers (one per exported horizon), the
  per-release-year min–max a band in the same colour at `alpha=0.3`. The
  pooled value can sit outside the band — it ranks the summed weights, which
  spread over more hexes than any one year's.
- **Missing export** — a regime whose CSV is absent is printed as `MISSING:`
  and left out (the figure draws what exists, so a one-regime interim works);
  no CSV at all raises.
- A computed table of the plotted values (pooled, per-year min and max, both
  quantities) is printed; no export of its own — the numbers are 026's.

`output_root/Figures/035/AreaTimeSeries_<season>_r<radius>m.png`, 180 mm ×
300 dpi; 032 picks the year-pooled season's figure. Considered drawing the
series inside 026; rejected, 026 runs one regime at a time and the
comparison is the point.

## Notebook 036 — Seeding map (special case)

[`036_SeedingMap`](../notebooks/036_SeedingMap.py) draws the release set for
the paper's Data section: the 872 REDLIST cells `010` seeds from
(`data/helcom_fucus_redlist/fucus_release_points.geojson`), as filled cell
polygons on 029's Baltic extent, in two side-by-side panels.

| panel | value | source |
|---|---|---|
| left | HELCOM level-2 subbasin of the cell | 024a key `helcom_subbasin` of the hex at the cell's representative point (029's mapping) |
| right | land-seeded share of the cell's particles (%) | first two time steps of every trajectory zarr of `land_regime` / `land_release_year` (default `surface_stokes` 2019, 73 releases) |

- **Land share from the zarrs.** 010 seeds cell-major, so trajectory `i`
  belongs to cell `i // particles_per_cell`; a particle is land-seeded under
  the pipeline's test (zero lon and lat displacement over the first step,
  024's `release_hex = -1`). The 024 counts store holds the land-seeded total
  but no cell id, so the zarrs are the only per-cell source; reading two
  time steps of 73 zarrs takes seconds. One regime-year (7 300 particles per
  cell) suffices: the BSH land mask is the same for every regime and release,
  and the first three 2019 releases of all three regimes agree within
  0.4 points (34.7–35.1 %). Pooled
  34.9 %, matching the sentinel share of raw `n_obs`.
- **Subbasin borders, not outlines.** The HELCOM polygons trace every island
  of the coast they were digitised against, which in the Åland and Stockholm
  archipelagos smears over the cells being shown; 036 draws only the edges
  two polygons share. The Kattegat's open edge to the Skagerrak, shared with
  no polygon, is therefore not drawn.
- **Override: categorical palette beyond the default cycle.** 16 subbasins
  hold cells; the default 10-colour cycle would give pairs of subbasins one
  colour and the legend could not separate them. The palette is the default
  cycle followed by `tab20`'s light tints of the same hues (the cycle's own
  paired extension), with grey (index 7 of both) held back for the 3
  `_outside` cells. Subbasins take colours in order of their cells' mean
  latitude, so a hue and its tint go to a southern and a northern basin,
  never to neighbours. Considered `tab20` in its native order; rejected, it
  puts each hue next to its tint and adjacent basins would read as one.
- **Override: subbasin borders** in mid-grey (`"0.5"`) at `linewidth=0.4`,
  so a border reads neither as coastline (black, as 025–033) nor as a
  category fill.
- **Land-share panel**: default colormap, fixed `Normalize(0, 100)`, colorbar
  right of the map (027/029/033 placement for the tall Baltic crop).
- **Legend** under both panels, three columns, entries in colour order with
  the cell count per subbasin.
- A table of cells, share of cells and mean land-seeded share per subbasin is
  printed; no export.

`output_root/Figures/036/SeedingMap_r<radius>m.png`, 180 mm × 300 dpi; first
entry of 032's figure list.

## Cross-references

- [beaching.md](beaching.md) — 029/031/033's store, model, and limitations.
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
