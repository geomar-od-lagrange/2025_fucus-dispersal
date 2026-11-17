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
# # Bottom friction: spread with and without bed friction
#
# How much does bed friction shrink the bottom-regime dispersal footprint?
# Compares the frictionless `bottom` regime against `bottom_ub49`, the same
# 010 runs with the lab bed-friction model (drift `max(|U| − U_b, 0)` along
# the bottom current, at rest below the break-off speed U_b = 49 mm/s;
# docs/bottom_friction.md), pooled over one release season of the
# `release_years`, as the other consumers pool. Only the releases present in
# **both** regimes enter (per year, the common release doys), so the two
# sides always describe the same release set. Reads the 024 counts
# partitions (the common release set, and the occupancy when unweighted), the
# 024f survocc partitions when `weighting` names a rate-model member, the
# 024a key and the 024e beaching partitions of one rate-model member; parquet
# only, no Dask.
#
# **Weighting**, as 026: `"unweighted"` maps the counts store's `n_obs`, a
# 024f member tag maps its survival-weighted `surv` (each sample weighted by
# its un-beached fraction `exp(−A)`, docs/survival_occupancy.md). Either is
# called `weight` below. Unweighted isolates friction's effect on transport;
# weighted adds the stranding each regime sees.
#
# **Footprint metric** (the R1 area of plans/done/horizons_and_weighting.md D6):
# for each cumulative elapsed-time horizon `T` (every `age_bin` whose window
# starts before `T`), rank target hexes by density `weight / water_area_m2`,
# descending, with zero-water hexes first at zero area; take the cumulative
# share of particle-time; the area holding a share `q` is the water area
# summed up to the crossing hex, the crossing hex taken fractionally (a 6 km
# hex is 93.5 km², so whole hexes would quantise small-`T` areas).
# `length_km = sqrt(area_km2)` is the side of a square of that area. Pooled
# over the release years and per release year, for the interannual range.
#
# **Stranded share** (029's stranding curve at the same horizons): the
# percentage of the released real drifters stranded at ages below `T`, from
# the 024e store at `member`, on the same release set.

# %%
import itertools
import re
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
# Read root of the data twin (BSH static H0 for the fine-grid bbox, the BSH
# wet-cell coastline, the HELCOM subbasins).
data_root = "../data"
# Read root of the hex-aggregate store; write root for figures and exports.
output_root = "../output"
# Hex radius of the store (built by 024a/024). Must match the files on disk.
hex_radius = 6000
# Reference (frictionless) and bed-friction regime names, as in
# Trajectories/<regime>/.
regime_ref = "bottom"
regime_var = "bottom_ub49"
# 024e rate-model member of the stranded-share table, the same for both
# regimes (docs/beaching.md).
member = "step_wc0p1_ts3_tcinf_tf1_tw0_b750"
# Weighting of the areas and maps: "unweighted" reads the 024 counts store
# (n_obs); a 024f rate-model member tag reads that member's survocc store
# (survival-weighted surv). Opaque — never parsed; tags every output name.
weighting = "unweighted"
# Release season: DJF, MAM, JJA, SON or ALL (every month).
season = "ALL"
# Release years pooled: "first-last" (inclusive) or one year ("2019").
# Partitions of other years stay on disk unopened (2025 is partial).
release_years = "2016-2024"
# Age-bin granularity of the counts store (must match 024).
age_bin_days = 10
# Reduction window (days): only ages in [0, max_age_days) are read, as in
# 026. A multiple of age_bin_days.
max_age_days = 180
# Cumulative elapsed-time horizons, in days; multiples of age_bin_days and
# ≤ max_age_days.
time_horizons_days_csv = "30,60,90,180"
# Shares of particle-time whose holding area is reported. 0.9 only: the 0.5
# area is set by particles stalling in a few coastal hexes
# (docs/visualisations.md).
quantiles_csv = "0.9"
# Horizon of the map pair; a multiple of age_bin_days, ≤ max_age_days.
map_horizon_days = 90
# Dilution maps leave a hex blank (NaN dilution) when less than this
# fraction of its area is water: a sliver's tiny water area turns a few
# samples into a dilution decades above its neighbours
# (docs/visualisations.md). Its percentage, and every non-dilution product,
# keep it.
dilution_min_water_fraction = 0.1
# Floor and top of the dilution colour scale, fraction of the released
# particle-time per km², whole decades. The same in every dilution figure
# (026, 026a, 034): the authors' fixed scale, two classes per decade
# (docs/visualisations.md). Every hex with weight is drawn: values below the
# floor take the under colour, values above the top the over colour; white is
# "no particles", or a hex below `dilution_min_water_fraction` water.
dilution_floor = 1e-9
dilution_top = 1e-3
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
#
# `SEASON_MONTHS` lives here rather than in the parameters cell because
# papermill only injects primitives; `season` selects one entry.

# %%
SEASON_MONTHS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ALL": list(range(1, 13)),
}
if season not in SEASON_MONTHS:
    raise ValueError(f"season {season!r} not one of {sorted(SEASON_MONTHS)}")
season_months = SEASON_MONTHS[season]
first_year, _, last_year = release_years.partition("-")
analysis_years = range(int(first_year), int(last_year or first_year) + 1)
assert len(analysis_years), f"release_years {release_years!r} is an empty range"

output_root = Path(output_root)
data_root = Path(data_root)
store_root = output_root / "HexAggregates"
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
quantiles = [float(x) for x in quantiles_csv.split(",") if x]
regimes = [regime_ref, regime_var]
# Plain-language figure titles; the run keys stay in filenames and prints.
REGIME_LABELS = {
    "bottom": "bottom",
    "bottom_ub49": "bottom, bed friction",
}
# The weighting in words, appended to the map titles as on every other map
# figure: an 024f rate-model member keeps the stranding loss.
weighting_label = f"{'without' if weighting == 'unweighted' else 'with'} stranding"
assert max_age_days % age_bin_days == 0, (
    f"max_age_days {max_age_days} d is not a multiple of age_bin_days {age_bin_days} d"
)
for h in [*time_horizons_days, map_horizon_days]:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d; "
        f"it would cut a bin in half"
    )
    assert 0 < h <= max_age_days, f"horizon {h} d outside (0, {max_age_days}] d"
n_age_bins = max_age_days // age_bin_days

# Print-ready output: every saved figure is one full text-width (180 mm)
# figure at 300 dpi, so it lands in the manuscript at its final size
# (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Subbasin border stroke, as 029 (docs/visualisations.md).
basin_border_lw = 0.52

# Shared dilution scale: two classes per decade from `dilution_floor` up to
# `dilution_top`. The physical bound is 1 per km² (all of the released
# particle-time in one km² of water); the scale stops where the data do
# (docs/visualisations.md).
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
# the colorbar's extend="both" triangles. White is "no particles", or a
# zero-water or sliver hex (below `dilution_min_water_fraction`).
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

print(f"reference {regime_ref!r} vs bed friction {regime_var!r}")
print(f"season {season} (months {season_months}), release years "
      f"{analysis_years[0]}–{analysis_years[-1]}, window [0, {max_age_days}) d")

# %% [markdown]
# # Common release set
#
# One partition per regime and release year,
# ``HexAgg_counts_r<radius>m_<regime>_<year>.parquet``. The comparison is
# only fair on the same releases, so it runs on the years both regimes have
# and, within a year, on the in-season release doys both partitions hold; a
# release one side lacks (a lost 010 run) is dropped from both and reported.
# The regex pins four year digits right after the regime, so `bottom` does
# not also take the `bottom_ub49` partitions.

# %%
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")


def counts_partitions(regime):
    """{release_year: path} of the regime's counts partitions in `release_years`."""
    part_re = re.compile(rf"HexAgg_counts_r{hex_radius}m_{regime}_(\d{{4}})\.parquet$")
    out = {}
    years_not_read = []
    for f in sorted(store_root.glob(f"HexAgg_counts_r{hex_radius}m_{regime}_*.parquet")):
        m = part_re.search(f.name)
        if m is None:
            continue
        if int(m[1]) in analysis_years:
            out[int(m[1])] = f
        else:
            years_not_read.append(int(m[1]))
    if years_not_read:
        print(f"{regime}: on disk but outside release_years, not read: {years_not_read}")
    return out


def season_doys(year):
    """The release doys of `year` in the season (leap-correct)."""
    days = pd.date_range(f"{year}-01-01", f"{year}-12-31")
    return set(days[days.month.isin(season_months)].dayofyear.tolist())


# dict.fromkeys de-duplicates while keeping order: regime_var == regime_ref
# (a dry run of this notebook before the variant exists) reads one regime.
partitions = {regime: counts_partitions(regime) for regime in dict.fromkeys(regimes)}
for regime, parts in partitions.items():
    print(f"{regime}: counts partitions for {sorted(parts)}")
common_years = sorted(set.intersection(*(set(p) for p in partitions.values())))
if not common_years:
    raise FileNotFoundError(
        f"no release year in {release_years} with counts partitions of both "
        f"{regimes} at {store_root} — run 024"
    )

common_doys = {}
for year in common_years:
    doys = {
        regime: set(
            pd.read_parquet(
                parts[year], columns=["release_doy"],
                filters=[("release_doy", "in", season_doys(year))],
            )["release_doy"].unique().tolist()
        )
        for regime, parts in partitions.items()
    }
    common_doys[year] = set.intersection(*doys.values())
    for regime, d in doys.items():
        if d - common_doys[year]:
            print(f"  {year}: {regime} releases {sorted(d - common_doys[year])} "
                  f"absent from the other regime, dropped")
    print(f"{year}: {len(common_doys[year])} common {season} releases")
years = [y for y in common_years if common_doys[y]]
if not years:
    raise ValueError(f"no common {season} releases in {common_years}")

# %% [markdown]
# # Read counts
#
# Each partition is read with the common doys and the age window pushed into
# the parquet filter and reduced at once to per-`(age_bin, target_hex)`
# sums, kept per release year for the interannual range.
#
# The `-1` sentinel in either hex column marks land-seeded particles and obs
# that left the BSH domain; those rows carry no geometry and are dropped from
# numerator and denominator alike. Released particle-time per age bin: every
# released water-seeded drifter contributes one sample per output step to
# age bin 0 (`release_hex ≥ 0`, sentinel targets included), the same
# denominator as 026.

# %%
counts = {}
released_per_bin = {}
for regime, parts in partitions.items():
    reduced = []
    released = {}
    for year in years:
        filters = [("release_doy", "in", common_doys[year]), ("age_bin", "<", n_age_bins)]
        if weighting == "unweighted":
            df = pd.read_parquet(
                parts[year], columns=["release_hex", "target_hex", "age_bin", "n_obs"],
                filters=filters,
            )
            valid = (df["release_hex"] >= 0) & (df["target_hex"] >= 0)
            released[year] = float(
                df.loc[(df["age_bin"] == 0) & (df["release_hex"] >= 0), "n_obs"].sum()
            )
            df = df.rename(columns={"n_obs": "weight"})
            source = parts[year].name
        else:
            # Monthly 024f partitions of the year; the season's months only.
            months = sorted({d.month for d in pd.to_datetime(
                [f"{year}-{doy:03d}" for doy in common_doys[year]], format="%Y-%j")})
            files = [
                store_root
                / f"HexAgg_survocc_r{hex_radius}m_{regime}_{year}_m{m:02d}_{weighting}.parquet"
                for m in months
            ]
            df = pd.concat(
                [pd.read_parquet(f, columns=["target_hex", "age_bin", "surv", "occ"],
                                 filters=filters) for f in files],
                ignore_index=True,
            )
            # 024f holds water-seeded drifters only; age-bin-0 occupancy is the
            # released particle-time, sentinel targets included (as 026).
            released[year] = float(df.loc[df["age_bin"] == 0, "occ"].sum())
            df = df.drop(columns="occ").rename(columns={"surv": "weight"})
            valid = df["target_hex"] >= 0
            source = f"{len(files)} survocc partitions of {regime} {year}"
        reduced.append(
            df[valid].groupby(["age_bin", "target_hex"], as_index=False)["weight"].sum()
            .assign(release_year=year)
        )
        print(
            f"{source}: {len(df):,} rows, sentinel (-1) share of weight "
            f"{100 * df.loc[~valid, 'weight'].sum() / df['weight'].sum():.2f} %"
        )
        del df, valid
    counts[regime] = pd.concat(reduced, ignore_index=True)
    released_per_bin[regime] = pd.Series(released)
    print(
        f"{regime}: {released_per_bin[regime].sum() / (age_bin_days * 24):,.0f} "
        f"water-seeded drifters released over {years}, weighting {weighting}"
    )

# %% [markdown]
# # Area holding a share of particle-time
#
# `area_to_share` is D6's ranking with the crossing hex taken fractionally.
# Zero-water hexes (all land in the BSH coastline, yet a sub-grid trajectory
# position fell inside) carry particle-time but no area, so they rank first:
# any share they hold costs nothing, which is the only ordering under which
# the density ranking is still "smallest area first". `scope` is "pooled"
# (every year) or one release year.

# %%
def area_to_share(weight, water_area_m2, q):
    """Smallest water area (m²) holding a share ``q`` of ``weight``, hexes
    ranked by ``weight / water_area_m2`` descending (zero-water hexes first)
    and the crossing hex counted by the fraction of it the share needs."""
    density = np.full_like(weight, np.inf)
    np.divide(weight, water_area_m2, out=density, where=water_area_m2 > 0)
    order = np.argsort(-density, kind="stable")
    w = weight[order] / weight.sum()
    a = water_area_m2[order]
    share_after = np.cumsum(w)
    # Clamp: at q = 1 rounding can leave the last cumulative share just short.
    i = min(int(np.searchsorted(share_after, q)), len(w) - 1)
    share_before = share_after[i] - w[i]
    return a[:i].sum() + a[i] * (q - share_before) / w[i]


water_area = key.set_index("hex_id")["water_area_m2"]
rows = []
for regime, c in counts.items():
    for scope in ["pooled", *years]:
        cs = c if scope == "pooled" else c[c["release_year"] == scope]
        for h in time_horizons_days:
            per_hex = cs[cs["age_bin"] * age_bin_days < h].groupby("target_hex")["weight"].sum()
            area = water_area.reindex(per_hex.index).to_numpy(dtype=float)
            assert np.isfinite(area).all(), "target hexes missing from the key"
            for q in quantiles:
                area_km2 = area_to_share(per_hex.to_numpy(dtype=float), area, q) / 1e6
                rows.append(
                    {"regime": regime, "scope": str(scope), "horizon_days": h, "q": q,
                     "area_km2": area_km2, "length_km": np.sqrt(area_km2)}
                )
areas = pd.DataFrame(rows)
ref = areas.loc[
    areas["regime"] == regime_ref, ["scope", "horizon_days", "q", "area_km2", "length_km"]
]
areas = areas.merge(ref, on=["scope", "horizon_days", "q"], suffixes=("", "_ref"))
areas["area_ratio_to_ref"] = areas["area_km2"] / areas.pop("area_km2_ref")
areas["length_ratio_to_ref"] = areas["length_km"] / areas.pop("length_km_ref")

csv_path = export_dir / f"034_areas_{season}_r{hex_radius}m_{weighting}.csv"
areas.to_csv(csv_path, index=False, float_format="%.6g")
print(f"wrote {csv_path}")
pooled = areas[areas["scope"] == "pooled"]
print(
    pooled.pivot_table(
        index=["horizon_days", "q"], columns="regime",
        values=["area_km2", "length_km"],
    ).round(1).to_string()
)
print()
var_ratio = areas[areas["regime"] == regime_var]
ratio_range = (
    var_ratio[var_ratio["scope"] != "pooled"]
    .groupby(["q", "horizon_days"])["length_ratio_to_ref"].agg(["min", "max"])
)
print("length ratio, pooled with per-release-year min–max:")
print(
    var_ratio[var_ratio["scope"] == "pooled"]
    .set_index(["q", "horizon_days"])[["area_ratio_to_ref", "length_ratio_to_ref"]]
    .join(ratio_range).round(3).to_string()
)

# %% [markdown]
# # Footprint length ratio vs. horizon
#
# One panel per share `q`: the bed-friction footprint length over the
# frictionless one, `length_km(regime_var) / length_km(regime_ref)`, against
# the horizon, pooled (line) with the per-release-year min–max (band). The
# ratio is the quantity the comparison asks for; two absolute lines hide it in
# the gap between them. `length_km` rather than area because it reads as a
# distance (the area ratio is its square, in the CSV).

# %%
length_ratio = (
    var_ratio.set_index(["scope", "q", "horizon_days"])["length_ratio_to_ref"].to_xarray()
)
length_ratio.attrs = {"long_name": "footprint length ratio"}
length_ratio["horizon_days"].attrs = {"long_name": "elapsed time T", "units": "d"}
per_year = length_ratio.sel(scope=[str(y) for y in years])

fig, axes = plt.subplots(1, len(quantiles), squeeze=False, layout="constrained")
# Print-ready width; one line needs no map aspect, so 0.45 of the page width
# per panel row reads as a normal wide chart.
fig.set_size_inches(FIGURE_WIDTH_IN, 0.45 * FIGURE_WIDTH_IN)
band_color = plt.rcParams["axes.prop_cycle"].by_key()["color"][0]
for ax, q in zip(axes.flat, quantiles):
    # The band takes the line's colour (the cycle's first) so line and band
    # read as one object; alpha keeps the line visible through it (as 035).
    ax.fill_between(
        length_ratio["horizon_days"], per_year.sel(q=q).min("scope"),
        per_year.sel(q=q).max("scope"), color=band_color, alpha=0.3,
    )
    length_ratio.sel(q=q, scope="pooled").plot.line(x="horizon_days", ax=ax)
    ax.set_title(
        f"{REGIME_LABELS.get(regime_var, regime_var)} / "
        f"{REGIME_LABELS.get(regime_ref, regime_ref)}, "
        f"area holding {100 * q:g} % of particle-time"
    )
    ax.set_xticks(time_horizons_days)
fig_path = figure_dir / f"BottomFriction_length_{season}_r{hex_radius}m_{weighting}.png"
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
# sit side by side (reference left) at the full page width, under one colour
# norm and one colourbar under the pair, so the same colour is the same value
# in both. `dilution`
# is the fixed log scale shared with 026 and 026a, below-floor hexes in the
# under colour; `percentage` is
# linear to the pooled 99th percentile. Extent is the Baltic crop of the
# parameters (the key's bounds if all four are 0). Coastline is the outline
# of the BSH HBMnoku wet cells, merged fine-first as in 024a (fine polygons
# plus the coarse ones outside the fine bbox, minus a 0.001° strip along the
# bbox where the grids meet in open water). Subbasin borders are
# the edges two HELCOM level-2 polygons share, without codes.

# %%
def to_hex_gdf(counts_df, key_df, released):
    """Per-target-hex relative density as a GeoDataFrame; ``released`` is
    the released particle-time (samples) over the same ages."""
    gdf = (
        counts_df.groupby("target_hex", as_index=False)["weight"].sum()
        .merge(
            key_df[["hex_id", "geometry", "area_m2", "water_area_m2"]],
            left_on="target_hex", right_on="hex_id",
        )
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )
    gdf["percentage"] = 100 * gdf["weight"] / gdf["weight"].sum()
    gdf["dilution"] = gdf["weight"] / released / (
        gdf["water_area_m2"].where(
            (gdf["water_area_m2"] > 0)
            & (gdf["water_area_m2"] >= dilution_min_water_fraction * gdf["area_m2"])
        ) / 1e6
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

# Subbasin borders as the edges two HELCOM polygons share (029's block): the
# polygon outlines trace every island of the coast they were digitised
# against, which would smear over the hexes. Clipped to the extent, as the
# hexes are.
subbasin_outline = (
    gpd.read_file(data_root / "helcom_subbasins_2022" / "HELCOM_subbasins_2022_level2.shp")
    .to_crs(4326).boundary
)
basin_borders = gpd.GeoSeries(
    [
        a.intersection(b)
        for a, b in itertools.combinations(subbasin_outline, 2)
        if a.intersects(b)
    ],
    crs="EPSG:4326",
).clip(box(domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max))
print(f"subbasin borders: {len(basin_borders)} shared edges between "
      f"{len(subbasin_outline)} HELCOM polygons")


def hex_map(gdf, column, ax, norm, cmap, title):
    """One relative-density panel on a plain lon/lat axis: hexes through
    ``norm`` and ``cmap`` (None: the default colormap), subbasin borders,
    coastline."""
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
    # Override: mid-grey thin borders (029's), so a subbasin border reads
    # neither as coastline (black) nor as a hex fill; under the coastline.
    basin_borders.plot(ax=ax, color="0.5", linewidth=basin_border_lw, zorder=1.5)
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(domain_lon_min, domain_lon_max)
    ax.set_ylim(domain_lat_min, domain_lat_max)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (domain_lat_min + domain_lat_max))))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)


# Released particle-time over [0, map_horizon_days), every release year pooled.
released_map = {
    regime: released_per_bin[regime].sum() * map_horizon_days / age_bin_days
    for regime in counts
}
gdfs = {
    regime: to_hex_gdf(
        c[c["age_bin"] * age_bin_days < map_horizon_days], key, released_map[regime],
    )
    for regime, c in counts.items()
}
for regime, g in gdfs.items():
    print(
        f"{regime}: remaining / released particle-time "
        f"{g['weight'].sum() / released_map[regime]:.4f}"
    )
both = pd.concat(gdfs.values())
for regime, g in gdfs.items():
    print(
        f"{regime}: {g.loc[g['dilution'] < dilution_floor, 'percentage'].sum():.3f} % of the "
        f"weight below the dilution floor (under colour), "
        f"{g.loc[g['dilution'].isna(), 'percentage'].sum():.3f} % on zero-water and sliver hexes (no dilution), "
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
    fig, axes = plt.subplots(1, 2, layout="constrained")
    # Height seed: tall enough that the maps are width-limited, so they span
    # the page width; the slack this leaves is measured and cut below.
    fig.set_size_inches(FIGURE_WIDTH_IN, FIGURE_WIDTH_IN)
    for ax, regime in zip(axes, regimes):
        hex_map(gdfs[regime], column, ax, norm, cmap, title=(
            f"{REGIME_LABELS.get(regime, regime)}, {map_horizon_days} d, {weighting_label}"
        ))
    # One colourbar under both maps, spanning the pair. pad 0.01: the default
    # is a fraction of the grid height and leaves a blank band; aspect 40
    # (default 20): at full page width the default bar is ~9 mm thick.
    fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), ax=axes, location="bottom",
        pad=0.01, aspect=40, label=COLUMN_LABEL[column], extend=extend,
        **colorbar_kwargs,
    )
    # The fixed-aspect maps fill their gridspec cells' width but not their
    # height: measure that vertical slack after one layout pass and cut it
    # from the figure height. The maps stay width-limited, so their size does
    # not change and the layout's fixed-inch padding needs no second pass.
    fig.draw_without_rendering()
    cell, drawn = axes[0].get_position(original=True), axes[0].get_position()
    assert np.isclose(drawn.width, cell.width, rtol=1e-3), "maps are height-limited; raise the seed"
    seed_height_in = fig.get_size_inches()[1]
    fig.set_size_inches(
        FIGURE_WIDTH_IN, seed_height_in - (cell.height - drawn.height) * seed_height_in
    )
    fig.draw_without_rendering()
    cell, drawn = axes[0].get_position(original=True), axes[0].get_position()
    print(f"{column}: {fig.get_size_inches()[1] * 25.4:.0f} mm tall, map "
          f"{drawn.width * FIGURE_WIDTH_IN * 25.4:.0f} mm wide, residual slack "
          f"{(cell.height - drawn.height) * fig.get_size_inches()[1] * 25.4:.1f} mm")
    fig_path = (
        figure_dir
        / f"BottomFriction_map_{season}_T{map_horizon_days}d_r{hex_radius}m_{weighting}_{column}.png"
    )
    fig.savefig(fig_path, dpi=FIGURE_DPI)
    print(f"wrote {fig_path}")
    plt.show()

# %% [markdown]
# # Stranded share
#
# One 024e partition per `(regime, year, month, member)`,
# `HexAgg_beaching_r<radius>m_<regime>_<year>_mMM_<member>.parquet`, columns
# per docs/beaching.md "Store schema". The release set is the counts one
# above narrowed to the `(year, release_doy)` both regimes' 024e partitions
# hold, so the stranded rows and the area rows describe the same releases
# (a release one 024e side lacks is dropped from both and reported). As 029,
#
# `stranded(T) = 100 · Σweight[beach_hex ≥ 0 & beach_age_bin < T / age_bin_days] / Σweight`,
#
# with Σweight the released real drifters (deposits + residual), pooled from
# the summed weights and per release year. Each partition's build attrs are
# taken per file before the concat (pandas drops attrs that differ between
# the frames) and must agree on the window, the binning and the rate-model
# point across both regimes: the comparison is one rate model on two drift
# regimes. A regime without 024e partitions at `member` is reported and the
# table skipped, so the area comparison still runs before 024d–e have.

# %%
# Keys every pooled partition of one member must share: the reduction window,
# the binning and the rate-model point. Provenance (builder, git_sha,
# slurm_job_id) and the partition's identity (regime, release_year,
# release_month) legitimately differ and are not compared.
POOLED_KEYS = (
    "max_age_days", "age_bin_days", "disp_bin_km", "w_c", "tau_strong_hours",
    "tau_calm_days", "delta", "trap_flat", "trap_wall", "band_m",
)


def pooled_attrs(attrs, store):
    """The build attrs every pooled partition shares, from one `df.attrs`
    per partition."""
    old = sum("max_age_days" not in a for a in attrs)
    if old:
        raise ValueError(
            f"{old} of {len(attrs)} {store} partition(s) carry no build attrs "
            f"(old schema, written before max_age_days was recorded) — rebuild "
            f"them with the current 024e reducer"
        )
    distinct = {tuple(a[k] for k in POOLED_KEYS) for a in attrs}
    if len(distinct) > 1:
        raise ValueError(
            f"pooled {store} partitions disagree on {POOLED_KEYS}: "
            f"{sorted(distinct)} — rebuild them at one setting"
        )
    return dict(zip(POOLED_KEYS, distinct.pop()))


def beaching_partitions(regime):
    """{release_year: [path, ...]} of the regime's in-season 024e partitions
    at `member` in the counts comparison's years."""
    part_re = re.compile(
        rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}})_m(\d{{2}})"
        rf"_{re.escape(member)}\.parquet$"
    )
    out = {}
    for f in sorted(store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet")):
        m = part_re.search(f.name)
        if m is None or int(m[1]) not in years or int(m[2]) not in season_months:
            continue
        out.setdefault(int(m[1]), []).append(f)
    return out


beaching_parts = {regime: beaching_partitions(regime) for regime in dict.fromkeys(regimes)}
for regime, by_year in beaching_parts.items():
    print(f"{regime}: 024e {member} partitions for {sorted(by_year)} "
          f"({sum(len(v) for v in by_year.values())} monthly)")
    if not by_year:
        print(f"MISSING: 024e {regime} {member}, season {season} — the stranded-share "
              f"table is NOT computed; run 024d/024e for {regime}")
beaching_ok = all(beaching_parts.values())

# %%
# Each partition is reduced at once to per-`(release_doy, stranding bin)`
# weight sums, the stranding bin being `beach_age_bin` where `beach_hex ≥ 0`
# and -1 for the never-beached residual: additive, and all the horizons need.
if beaching_ok:
    beaching = {}
    beaching_attrs = []
    for regime, parts_by_year in beaching_parts.items():
        reduced = []
        for year, files in parts_by_year.items():
            for f in files:
                df = pd.read_parquet(
                    f, columns=["release_doy", "beach_hex", "beach_age_bin", "weight"]
                )
                beaching_attrs.append(df.attrs)
                df = df[df["release_doy"].isin(common_doys[year])]
                reduced.append(
                    df.assign(stranded_bin=df["beach_age_bin"].where(df["beach_hex"] >= 0, -1))
                    .groupby(["release_doy", "stranded_bin"], as_index=False)["weight"].sum()
                    .assign(release_year=year)
                )
        beaching[regime] = pd.concat(reduced, ignore_index=True)

    window = pooled_attrs(beaching_attrs, f"024e {member} ({', '.join(beaching_parts)})")
    assert window["age_bin_days"] == age_bin_days, (
        f"024e age_bin_days {window['age_bin_days']} != {age_bin_days}"
    )
    assert max(time_horizons_days) <= window["max_age_days"], (
        f"horizon {max(time_horizons_days)} d is past the 024e reduction window "
        f"({window['max_age_days']} d) — rebuild 024e or lower the horizon"
    )

    # Narrow the counts release set to the (year, release_doy) both 024e
    # sides hold.
    release_sets = {
        regime: set(zip(b["release_year"], b["release_doy"])) for regime, b in beaching.items()
    }
    beaching_releases = set.intersection(*release_sets.values())
    if not beaching_releases:
        raise ValueError(f"no common {season} release in the 024e partitions of {regimes}")
    for regime, own in release_sets.items():
        dropped = sorted(own - beaching_releases)
        lacking = sorted((y, d) for y in years for d in common_doys[y] if (y, d) not in own)
        if dropped:
            print(f"{regime}: {len(dropped)} 024e releases absent from the other "
                  f"regime, dropped: {dropped}")
        if lacking:
            print(f"{regime}: {len(lacking)} common counts releases without 024e rows: {lacking}")
    beaching = {
        regime: b[[pair in beaching_releases for pair in zip(b["release_year"], b["release_doy"])]]
        for regime, b in beaching.items()
    }
    print(f"stranded share on {len(beaching_releases)} common {season} releases; "
          f"window [0, {window['max_age_days']}) d, band {window['band_m']} m")

# %%
if beaching_ok:
    rows = []
    for regime, b in beaching.items():
        for h in time_horizons_days:
            is_stranded = (b["stranded_bin"] >= 0) & (b["stranded_bin"] < h // age_bin_days)
            per_year = (
                b[["release_year", "weight"]]
                .assign(stranded=b["weight"].where(is_stranded, 0.0))
                .groupby("release_year")[["weight", "stranded"]].sum()
            )
            per_scope = pd.concat(
                [per_year.rename(index=str), per_year.sum().to_frame("pooled").T]
            )
            for scope, r in per_scope.iterrows():
                rows.append(
                    {"regime": regime, "scope": scope, "horizon_days": h,
                     "released_weight": r["weight"],
                     "stranded_pct": 100 * r["stranded"] / r["weight"]}
                )
    stranded = pd.DataFrame(rows)
    ref = stranded.loc[
        stranded["regime"] == regime_ref, ["scope", "horizon_days", "stranded_pct"]
    ]
    stranded = stranded.merge(ref, on=["scope", "horizon_days"], suffixes=("", "_ref"))
    stranded["stranded_pct_minus_ref"] = stranded["stranded_pct"] - stranded.pop(
        "stranded_pct_ref"
    )

    stranded_csv_path = export_dir / f"034_stranded_{season}_r{hex_radius}m_{member}.csv"
    stranded.to_csv(stranded_csv_path, index=False, float_format="%.6g")
    print(f"wrote {stranded_csv_path}")
    stranded_range = (
        stranded[stranded["scope"] != "pooled"]
        .groupby(["horizon_days", "regime"])["stranded_pct"].agg(["min", "max"])
    )
    print(f"stranded % of released, season {season}, member {member}; "
          "pooled with per-release-year min–max:")
    print(
        stranded[stranded["scope"] == "pooled"]
        .set_index(["horizon_days", "regime"])[["stranded_pct", "stranded_pct_minus_ref"]]
        .join(stranded_range).round(1).to_string()
    )

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
        f"sum(weight)={g['weight'].sum():,.0f}, "
        f"top-10 hexes hold {g['percentage'].nlargest(10).sum():.2f} %, "
        f"outside the map extent {g.loc[outside, 'percentage'].sum():.3f} %"
    )
