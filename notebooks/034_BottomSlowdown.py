# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: py:percent,md,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Bottom slowdown: spread at 0.8× bottom velocity
#
# How much does the bottom-regime dispersal footprint shrink if the bottom-cell
# current is assumed 20 % slower? Compares the 1× `bottom` regime against the
# `bottom_vs0p8` variant (010 with `velocity_scale = 0.8`, plain
# `AdvectionRK4`, no drag model) on the **same release doys of one year** —
# the variant is only run for one release per month, so the 1× side is
# restricted to those doys too. Reads the 024 counts partitions and the 024a
# key; parquet only, no Dask.
#
# **Footprint metric** (the R1 area of plans/done/horizons_and_weighting.md D6):
# for each cumulative elapsed-time horizon `T` (every `age_bin` whose window
# starts before `T`), rank target hexes by density `n_obs / water_area_m2`,
# descending, with zero-water hexes first at zero area; take the cumulative
# share of particle-time; the area holding a share `q` is the water area
# summed up to the crossing hex, the crossing hex taken fractionally (a 6 km
# hex is 93.5 km², so whole hexes would quantise small-`T` areas).
# `length_km = sqrt(area_km2)` is the side of a square of that area.

# %%
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, Normalize
from matplotlib.ticker import FuncFormatter
from shapely.geometry import box
from shapely.ops import unary_union

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the data twin (BSH static H0 for the isobath overlay).
data_root = "../data"
# Read root of the hex-aggregate store; write root for figures and exports.
output_root = "../output"
# Hex radius of the store (built by 024a/024). Must match the files on disk.
hex_radius = 6000
# Release year compared; both regimes' counts partitions for it must exist.
year = 2019
# Release doys compared, restricted on both sides. Default: the first
# doy = 1 + 5n of every month, the set the 0.8x variant was run for.
release_doys_csv = "1,36,61,91,121,156,186,216,246,276,306,336"
# Reference (1x) and variant regime names, as in Trajectories/<regime>/.
regime_ref = "bottom"
regime_var = "bottom_vs0p8"
# Cumulative elapsed-time horizons, in days; multiples of age_bin_days.
time_horizons_days_csv = "30,60,90,180"
# Age-bin granularity of the counts store (must match 024).
age_bin_days = 10
# Shares of particle-time whose holding area is reported. 0.9 only: the 0.5
# area is set by particles stalling in a few coastal hexes
# (docs/visualisations.md).
quantiles_csv = "0.9"
# Horizon of the map pair; a multiple of age_bin_days.
map_horizon_days = 90
# Isobath to overlay on the maps, in metres of BSH H0 (floor position, see
# docs/h0_semantics.md). Fucus vesiculosus no longer grows below ~3 m.
isobath_m = 3.0
# Floor and top of the dilution colour scale, fraction of the released
# particle-time per km², whole decades. The same in every dilution figure
# (026, 026a, 034): the authors' fixed scale, two classes per decade
# (docs/visualisations.md). Every hex with weight is drawn: values below the
# floor take the under colour, values above the top the over colour; white is
# only "no particles".
dilution_floor = 1e-7
dilution_top = 0.1
# Map extent (degrees). The 024a key tiles the whole BSH domain including the
# North Sea, which the bottom regime does not reach; cropping to the Baltic
# (western Baltic and Belt Sea to the Gulf of Finland and Bothnian Bay) spends
# the page on the hexes that carry particle-time. Set all four
# to 0 to fall back to the key's full bounds.
extent_lon_min = 9.0
extent_lon_max = 30.7
extent_lat_min = 53.0
extent_lat_max = 66.0

# %% [markdown]
# # Parse parameters

# %%
output_root = Path(output_root)
data_root = Path(data_root)
store_root = output_root / "HexAggregates"
release_doys = sorted(int(x) for x in release_doys_csv.split(",") if x)
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
quantiles = [float(x) for x in quantiles_csv.split(",") if x]
regimes = [regime_ref, regime_var]
# Plain-language figure titles; the run keys stay in filenames and prints.
REGIME_LABELS = {
    "bottom": "bottom",
    "bottom_vs0p8": "bottom, velocity × 0.8",
}
for h in [*time_horizons_days, map_horizon_days]:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d; "
        f"it would cut a bin in half"
    )

# Print-ready output: every saved figure is one full text-width (180 mm)
# figure at 300 dpi, so it lands in the manuscript at its final size
# (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# The map pair is two maps stacked; the Baltic is taller than wide, so two
# maps at full page width would be ~2.4x the page width tall. The figure is
# capped at a height that fits one manuscript page with its caption, and the
# maps take the width that height leaves them.
MAP_PAIR_HEIGHT_IN = 220 / 25.4

# Shared dilution scale: two classes per decade from `dilution_floor` up to
# `dilution_top` (at most 1 per km², all of the released particle-time in
# one km² of water).
DILUTION_CLASSES_PER_DECADE = 2
for bound in (dilution_floor, dilution_top):
    assert np.isclose(np.log10(bound), round(np.log10(bound))), f"{bound} is not a whole decade"
dilution_bounds = np.logspace(
    np.log10(dilution_floor), np.log10(dilution_top),
    round(DILUTION_CLASSES_PER_DECADE * np.log10(dilution_top / dilution_floor)) + 1,
)
# turbo, not the default colormap (docs/visualisations.md): dilution is
# multiplied by a per-site released biomass downstream, so it has to read in
# absolute terms over many decades — one high-contrast colour per half
# decade. Every hex with weight is drawn (masking sub-floor hexes made
# material look as if it vanished): below the floor in the under colour, above
# the top in the over colour — the lowest and highest class colours, shown as
# the colorbar's extend="both" triangles. White is only "no particles".
dilution_cmap = plt.get_cmap("turbo").resampled(len(dilution_bounds) - 1)
dilution_norm = BoundaryNorm(dilution_bounds, dilution_cmap.N)
# Tick labels on the whole decades only; the half-decade edges stay unlabelled.
decade_format = FuncFormatter(
    lambda x, _: f"$10^{{{round(np.log10(x))}}}$"
    if np.isclose(np.log10(x), round(np.log10(x))) else ""
)

figure_dir = output_root / "Figures" / "034"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "034"
export_dir.mkdir(parents=True, exist_ok=True)

print(f"year {year}, {len(release_doys)} release doys: {release_doys}")
print(f"reference {regime_ref!r} vs variant {regime_var!r}")

# %% [markdown]
# # Read key + counts
#
# One partition per regime,
# ``HexAgg_counts_r<radius>m_<regime>_<year>.parquet``, read with the release
# doys pushed into the parquet filter and reduced at once to per-`(age_bin,
# target_hex)` sums. Every requested doy must be present in both partitions —
# the comparison is only fair on the same release set, so a missing doy
# raises rather than shrinking one side.
#
# The `-1` sentinel in either hex column marks land-seeded particles and obs
# that left the BSH domain; those rows carry no geometry and are dropped from
# numerator and denominator alike.

# %%
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")

# dict.fromkeys de-duplicates while keeping order: regime_var == regime_ref
# (a dry run of this notebook before the variant exists) reads one partition.
counts = {}
# Released particle-time per age bin: every released water-seeded drifter
# contributes one sample per output step to age bin 0 (`release_hex ≥ 0`,
# sentinel targets included), the same denominator as 026.
released_per_bin = {}
for regime in dict.fromkeys(regimes):
    path = store_root / f"HexAgg_counts_r{hex_radius}m_{regime}_{year}.parquet"
    df = pd.read_parquet(
        path,
        columns=["release_hex", "release_doy", "target_hex", "age_bin", "n_obs"],
        filters=[("release_doy", "in", release_doys)],
    )
    missing = sorted(set(release_doys) - set(df["release_doy"].unique()))
    if missing:
        raise ValueError(f"{path.name}: release doys {missing} absent")
    valid = (df["release_hex"] >= 0) & (df["target_hex"] >= 0)
    released_per_bin[regime] = float(
        df.loc[(df["age_bin"] == 0) & (df["release_hex"] >= 0), "n_obs"].sum()
    )
    counts[regime] = (
        df[valid].groupby(["age_bin", "target_hex"], as_index=False)["n_obs"].sum()
    )
    print(
        f"{path.name}: {len(df):,} rows, sentinel (-1) share of n_obs "
        f"{100 * df.loc[~valid, 'n_obs'].sum() / df['n_obs'].sum():.2f} %, "
        f"max age_bin {int(df['age_bin'].max())}"
    )
    del df

# %% [markdown]
# # Area holding a share of particle-time
#
# `area_to_share` is D6's ranking with the crossing hex taken fractionally.
# Zero-water hexes (all land in the BSH coastline, yet a sub-grid trajectory
# position fell inside) carry particle-time but no area, so they rank first:
# any share they hold costs nothing, which is the only ordering under which
# the density ranking is still "smallest area first".

# %%
def area_to_share(n_obs, water_area_m2, q):
    """Smallest water area (m²) holding a share ``q`` of ``n_obs``, hexes
    ranked by ``n_obs / water_area_m2`` descending (zero-water hexes first)
    and the crossing hex counted by the fraction of it the share needs."""
    density = np.full_like(n_obs, np.inf)
    np.divide(n_obs, water_area_m2, out=density, where=water_area_m2 > 0)
    order = np.argsort(-density, kind="stable")
    w = n_obs[order] / n_obs.sum()
    a = water_area_m2[order]
    share_after = np.cumsum(w)
    # Clamp: at q = 1 rounding can leave the last cumulative share just short.
    i = min(int(np.searchsorted(share_after, q)), len(w) - 1)
    share_before = share_after[i] - w[i]
    return a[:i].sum() + a[i] * (q - share_before) / w[i]


water_area = key.set_index("hex_id")["water_area_m2"]
rows = []
for regime, c in counts.items():
    for h in time_horizons_days:
        per_hex = c[c["age_bin"] * age_bin_days < h].groupby("target_hex")["n_obs"].sum()
        area = water_area.reindex(per_hex.index).to_numpy(dtype=float)
        assert np.isfinite(area).all(), "target hexes missing from the key"
        for q in quantiles:
            area_km2 = area_to_share(per_hex.to_numpy(dtype=float), area, q) / 1e6
            rows.append(
                {"regime": regime, "horizon_days": h, "q": q,
                 "area_km2": area_km2, "length_km": np.sqrt(area_km2)}
            )
areas = pd.DataFrame(rows)
ref = areas.loc[areas["regime"] == regime_ref, ["horizon_days", "q", "area_km2", "length_km"]]
areas = areas.merge(ref, on=["horizon_days", "q"], suffixes=("", "_ref"))
areas["area_ratio_to_ref"] = areas["area_km2"] / areas.pop("area_km2_ref")
areas["length_ratio_to_ref"] = areas["length_km"] / areas.pop("length_km_ref")

csv_path = export_dir / f"034_areas_{year}.csv"
areas.to_csv(csv_path, index=False, float_format="%.6g")
print(f"wrote {csv_path}")
print(
    areas.pivot_table(
        index=["horizon_days", "q"], columns="regime",
        values=["area_km2", "length_km"],
    ).round(1).to_string()
)
print()
print(
    areas[areas["regime"] == regime_var]
    .set_index(["horizon_days", "q"])[["area_ratio_to_ref", "length_ratio_to_ref"]]
    .round(3).to_string()
)

# %% [markdown]
# # Footprint length ratio vs. horizon
#
# One panel per share `q`: the variant's footprint length over the
# reference's, `length_km(regime_var) / length_km(regime_ref)`, against the
# horizon. The ratio is the quantity the comparison asks for; two absolute
# lines a few per cent apart hide it in the gap between them. `length_km`
# rather than area because it reads as a distance (the area ratio is its
# square, in the CSV). One release set on each side, so there is one ratio per
# horizon and no per-release spread.

# %%
length_ratio = (
    areas[areas["regime"] == regime_var]
    .set_index(["q", "horizon_days"])["length_ratio_to_ref"].to_xarray()
)
length_ratio.attrs = {"long_name": "footprint length ratio"}
length_ratio["horizon_days"].attrs = {"long_name": "elapsed time T", "units": "d"}

fig, axes = plt.subplots(1, len(quantiles), squeeze=False, layout="constrained")
# Print-ready width; one line needs no map aspect, so 0.45 of the page width
# per panel row reads as a normal wide chart.
fig.set_size_inches(FIGURE_WIDTH_IN, 0.45 * FIGURE_WIDTH_IN)
for ax, q in zip(axes.flat, quantiles):
    length_ratio.sel(q=q).plot.line(x="horizon_days", ax=ax)
    # Ratio 1 = no slowdown effect. Drawn in the axes-frame colour so it reads
    # as a reference, not as a second data series in the cycle's first colour.
    ax.axhline(1.0, color=plt.rcParams["axes.edgecolor"], zorder=0)
    ax.set_title(f"{REGIME_LABELS.get(regime_var, regime_var)} / {REGIME_LABELS.get(regime_ref, regime_ref)}, area holding {100 * q:g} % of particle-time")
    ax.set_xticks(time_horizons_days)
fig_path = figure_dir / f"BottomSlowdown_length_{year}_r{hex_radius}m.png"
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Map pair at one horizon
#
# 026's rendering, notebook-local: `percentage` is the share of the regime's
# particle-time present in `[0, map_horizon_days)`, `dilution` the hex's
# particle-time as a fraction of the regime's released particle-time over
# the same ages (`N_released · T`) per km² of hex water area. The two regimes
# are stacked (reference on top) under one colour norm and one colourbar
# under the pair, so the same colour is the same value in both. `dilution`
# is the fixed log scale shared with 026 and 026a, below-floor hexes in the
# under colour; `percentage` is
# linear to the pooled 99th percentile. Extent is the Baltic crop of the
# parameters (the key's bounds if all four are 0). Coastline is the outline
# of the BSH HBMnoku wet cells, merged fine-first as in 024a (fine polygons
# plus the coarse ones outside the fine bbox, minus a 0.001° strip along the
# bbox where the grids meet in open water); the isobath is the BSH H0 3 m
# contour, fine grid taking precedence where it covers.

# %%
def to_hex_gdf(counts_df, key_df, released):
    """Per-target-hex relative density as a GeoDataFrame; ``released`` is
    the released particle-time (samples) over the same ages."""
    gdf = (
        counts_df.groupby("target_hex", as_index=False)["n_obs"].sum()
        .merge(
            key_df[["hex_id", "geometry", "water_area_m2"]],
            left_on="target_hex", right_on="hex_id",
        )
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )
    gdf["percentage"] = 100 * gdf["n_obs"] / gdf["n_obs"].sum()
    gdf["dilution"] = gdf["n_obs"] / released / (
        gdf["water_area_m2"].where(gdf["water_area_m2"] > 0) / 1e6
    )
    return gdf


COLUMN_LABEL = {
    "percentage": "share of particle-time (%)",
    "dilution": "fraction of released particle-time per km²",
}

if any((extent_lon_min, extent_lon_max, extent_lat_min, extent_lat_max)):
    domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max = (
        extent_lon_min, extent_lat_min, extent_lon_max, extent_lat_max
    )
else:
    domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max = key.total_bounds

h0_fine = xr.open_dataset(
    data_root / "bsh_hbmnoku_static/static_file_fine/H0_file_fine.nc"
).H0
h0_coarse = xr.open_dataset(
    data_root / "bsh_hbmnoku_static/static_file_coarse/H0_file_coarse.nc"
).H0
h0_coarse_outside_fine = h0_coarse.where(
    ~(
        (h0_coarse.lon >= h0_fine.lon.min()) & (h0_coarse.lon <= h0_fine.lon.max())
        & (h0_coarse.lat >= h0_fine.lat.min()) & (h0_coarse.lat <= h0_fine.lat.max())
    )
)

wet = gpd.read_file(data_root / "bsh_hbmnoku_static/coastline.geojson")
# Fine bbox padded by half a fine cell, so every fine cell lies inside it.
fine_bbox = box(
    float(h0_fine.lon.min()) - 0.5 * abs(float(h0_fine.lon.diff("lon").mean())),
    float(h0_fine.lat.min()) - 0.5 * abs(float(h0_fine.lat.diff("lat").mean())),
    float(h0_fine.lon.max()) + 0.5 * abs(float(h0_fine.lon.diff("lon").mean())),
    float(h0_fine.lat.max()) + 0.5 * abs(float(h0_fine.lat.diff("lat").mean())),
)
wet_fine_first = unary_union([
    unary_union(wet.loc[wet["grid"] == "fine", "geometry"].tolist()),
    unary_union(wet.loc[wet["grid"] == "coarse", "geometry"].tolist()).difference(fine_bbox),
])
coast = (
    gpd.GeoSeries([wet_fine_first], crs="EPSG:4326").boundary
    .difference(fine_bbox.exterior.buffer(1e-3))
    .clip(box(domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max))
)


def hex_map(gdf, column, ax, norm, cmap, title):
    """One relative-density panel on a plain lon/lat axis: hexes through
    ``norm`` and ``cmap`` (None: the default colormap), coastline, H0
    isobath."""
    # Only hexes touching the extent: off-extent patches still count in
    # constrained_layout's bounding box and push the title off the page.
    gdf = gdf.cx[domain_lon_min:domain_lon_max, domain_lat_min:domain_lat_max]
    gdf.plot(
        ax=ax, column=column, norm=norm, cmap=cmap,
        missing_kwds={"color": "none"},
        # Default polygon stroke is a black hairline that dominates at Baltic
        # zoom; matching edge to fill makes the hexes read as a field.
        edgecolor="face", linewidth=0.4, zorder=1,
    )
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    # Magenta hairline isobath: anti-viridis, legible over every fill colour.
    # Drawn last at zorder 10, above every other artist: where it runs along
    # the coast it must sit on the black line, not under it.
    for h0 in (h0_fine, h0_coarse_outside_fine):
        ax.contour(
            h0.lon, h0.lat, h0, levels=[isobath_m],
            colors="magenta", linewidths=0.3, zorder=10,
        )
    ax.set_xlim(domain_lon_min, domain_lon_max)
    ax.set_ylim(domain_lat_min, domain_lat_max)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (domain_lat_min + domain_lat_max))))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)


gdfs = {
    regime: to_hex_gdf(
        c[c["age_bin"] * age_bin_days < map_horizon_days], key,
        released_per_bin[regime] * map_horizon_days // age_bin_days,
    )
    for regime, c in counts.items()
}
for regime, g in gdfs.items():
    print(
        f"{regime}: remaining / released particle-time "
        f"{g['n_obs'].sum() / (released_per_bin[regime] * map_horizon_days // age_bin_days):.4f}"
    )
both = pd.concat(gdfs.values())
for regime, g in gdfs.items():
    print(
        f"{regime}: {g.loc[g['dilution'] < dilution_floor, 'percentage'].sum():.3f} % of the "
        f"weight below the dilution floor (under colour), "
        f"{g.loc[g['dilution'].isna(), 'percentage'].sum():.3f} % on zero-water hexes, "
        f"{g.loc[g['dilution'] > dilution_top, 'percentage'].sum():.3f} % above the top"
    )
# Shared across the pair: percentage saturates at the 99th percentile of the
# pooled panels (a few hexes sit decades above the bulk; the arrow marks
# it), dilution is the fixed log scale, below-floor hexes in the under colour.
VIEWS = {
    "percentage": (
        Normalize(vmin=0, vmax=both["percentage"].quantile(0.99)), None, "max", {},
    ),
    "dilution": (dilution_norm, dilution_cmap, "both", {"format": decade_format}),
}
for column, (norm, cmap, extend, colorbar_kwargs) in VIEWS.items():
    fig, axes = plt.subplots(2, 1, squeeze=False, layout="constrained")
    fig.set_size_inches(FIGURE_WIDTH_IN, MAP_PAIR_HEIGHT_IN)
    for ax, regime in zip(axes.flat, regimes):
        hex_map(gdfs[regime], column, ax, norm, cmap, title=f"{REGIME_LABELS.get(regime, regime)}, {map_horizon_days} d")
    # One colourbar for both maps, inset under the lower one at axes-fraction
    # width 1.0: the stacked maps are height-limited, so a bar under the
    # whole figure would run past them to the page edges.
    fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        cax=axes[-1, 0].inset_axes([0.0, -0.07, 1.0, 0.03]),
        orientation="horizontal", label=COLUMN_LABEL[column], extend=extend,
        **colorbar_kwargs,
    )
    fig_path = (
        figure_dir
        / f"BottomSlowdown_map_{year}_T{map_horizon_days}d_r{hex_radius}m_{column}.png"
    )
    fig.savefig(fig_path, dpi=FIGURE_DPI)
    print(f"wrote {fig_path}")
    plt.show()

# %% [markdown]
# # Validation prints

# %%
for regime, g in gdfs.items():
    # Particle-time the map extent cuts off, by a point inside each hex.
    c = g.geometry.representative_point()
    outside = ~(
        c.x.between(domain_lon_min, domain_lon_max)
        & c.y.between(domain_lat_min, domain_lat_max)
    )
    print(
        f"{regime}: {len(g):,} occupied hexes at {map_horizon_days} d, "
        f"sum(n_obs)={int(g['n_obs'].sum()):,}, "
        f"top-10 hexes hold {g['percentage'].nlargest(10).sum():.2f} %, "
        f"outside the map extent {g.loc[outside, 'percentage'].sum():.3f} %"
    )
