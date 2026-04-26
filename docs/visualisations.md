# Visualisations: per-plot rationale

Why each notebook in `020`–`025` shows what it shows, and where
styling overrides earn their deviation from
[../AGENTS.md](../AGENTS.md)'s plotting rules. AGENTS.md states the
**rules** (no `cmap=`, no `figsize=`, no `color=` literals); this doc
states the **rationale** for each surviving rule-break. Anything not
listed here should be dropped.

## Plot types and scopes at a glance

| Notebook | Plot type           | Per-regime?                   | Default scopes                                                            |
|----------|---------------------|-------------------------------|----------------------------------------------------------------------------|
| 020      | Raw trajectory lines| Overlaid                      | Per HELCOM release subbasin · German waters · Per release quarter          |
| 021      | Time-series stats   | Overlaid                      | Global only                                                                |
| 022      | Distance vs. time   | Overlaid (`hue="regime"`)     | Global · Per subbasin · German waters · Per quarter                        |
| 023      | Density + mean-age maps | Per-run (papermill regime) | Whole Baltic · Per subbasin · German waters · Per quarter                  |
| 025      | Hex density maps    | Per-run (regime + year + radius) | Panel A: Baltic · Panel B: per subbasin · Panel C: DE-zoom · Panel D: per quarter |

## Cross-cutting choices

**Per-regime: overlay vs per-run.** Overlay regimes when the
comparison *is* the point (lines, distance-vs-time). Render per-regime
when mixing would muddle the spatial signal — density and mean-age
maps fall here, because cross-panel reading invites artefacts from
regime physics rather than from the spatial field.

**Temporal scoping: per-quarter and per-year are alternatives, not
crossed.** Quarter pools across years; year pools across quarters.
Crossing them at every plot site multiplies panel counts beyond useful
inspection density. Quarter is derived from `release_doy` via
`pd.to_datetime(format="%Y%j").dt.quarter` (the integer-floor
`(doy-1)//90+1` is wrong by ≥1 quarter near boundaries because Q1–Q4
month lengths differ).

**Histogram dim layout: keep `obs` as a dim** in 023's `xhist` outputs.
Every scoped plot is then a `.sel`/`.sum`/`.groupby` on the lazy array,
computed once in a shared dask pass. Mean-age reuses the same lazy
object via an obs-weighted variant.

**Land-seeded filter** is the same `(diff lon at obs=0 == 0) &
(diff lat at obs=0 == 0)` mask in every viz; rationale and edge cases
in [distance_calculation.md](distance_calculation.md).

**Coastline overlay** in 025 is Natural Earth 10 m via cartopy's
shapereader, clipped once per extent. Natural Earth is preferred to
the BSH wet-cell coastline because the audience reads geographic
context, not model-grid registration. The cartopy notebooks
(020/022/023) use `add_feature(cfeature.COASTLINE)` per the AGENTS.md
default.

## Notebook 025 — Hex heatmaps (special case)

025 plots on a plain (non-cartopy) lon/lat axis. Hex polygons,
Natural Earth coastline, and HELCOM subbasin overlay are all already
in EPSG:4326; cartopy reprojection would introduce sub-pixel
mismatches that double-paint hex edges and create coastline-vs-overlay
drift. The cost is no scale-bar / gridlines from cartopy machinery,
but the audience reads density not bearing, so the plain axis is
cleaner.

### Justified styling overrides

The hex maps carry four overrides defended below. Aspect-driven
`figsize=(panel_height_in*aspect, panel_height_in)` is treated as
parametric layout (cartographic correctness — 1° lon at mid-latitude
visually equals 1° lat), not styling — same role cartopy plays for
020/022/023.

| Override                         | Where                              | Rationale                                                                                  |
|----------------------------------|------------------------------------|--------------------------------------------------------------------------------------------|
| `cmap = "viridis"`               | parameters cell, all panels        | Log-density spans 4–5 decades; perceptually uniform colormap is load-bearing. Named so a future reader doesn't swap in a non-uniform map. |
| `edgecolor="face"`, `linewidth=0.4` | hex polygons in `log_density_plot` | Default polygon stroke is a black hairline that dominates at Baltic-wide zoom; matching the edge to the fill makes the grid read as a continuous density field. |
| `color="black"`, `linewidth=0.5` | coastline in `log_density_plot`    | Black reads against the dark band of viridis (low density), which is the band where coastline registration matters most. Regime-aware colour was rejected as visual noise. |
| `color="magenta"`, `linewidth=1.05` | subbasin overlay in `log_density_plot`, Panel B only | Magenta is the canonical anti-viridis; sits at the colour-wheel position viridis avoids, so it stands clear of every fill colour and the black coastline. |

## Cross-references

- [seeding.md](seeding.md) — release-set semantics every viz reads.
- [distance_calculation.md](distance_calculation.md) — 022's metric.
- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) —
  the store 025 reads.
- [../AGENTS.md](../AGENTS.md) — the styling rules this doc defends
  exceptions to.
