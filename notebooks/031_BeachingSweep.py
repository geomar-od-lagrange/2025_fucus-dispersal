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
# # Beaching parameter sweep: reporting a range, not a number
#
# Parquet-only consumer that pools the `w_half` members of the beaching store
# written by `024d_BuildBeaching`. In the Baltic the beaching scheme can
# dominate the answer, so the headline stranding number is only meaningful as
# a **range over the rate parameters** — this notebook produces that range and
# the pattern metric that actually discriminates between members.
#
# `w_half` is the onshore-Stokes half-saturation in `s(w) = w/(w + w_half)`,
# and with `trap` degenerate the rate is `τ = τ0/s(w_onshore)` — so `w_half`
# and `τ0` are the whole rate model. They are **not independent**: rewriting
# `τ = τ0 + τ0·w_half/w` shows `τ0` as an additive floor and the *product*
# `τ0·w_half` as the weak-wave coefficient. Where `w ≪ w_half` only the
# product matters, so the two trade off along a ridge and the total beached
# fraction alone cannot separate them.
#
# What *does* separate them is **spatial selectivity**: large `w_half` keeps
# `g` far from saturation, so stranding concentrates where onshore waves are
# strong; small `w_half` saturates `g → 1` everywhere and the map collapses
# to a pure near-shore-residence field. So this notebook reports, per member:
#
# 1. **total beached fraction** — the range (degenerate along the ridge);
# 2. **concentration (Gini) of stranded weight across coastal hexes** — the
#    discriminating statistic;
# 3. **where-stranded maps** side by side across members.

# %%
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
output_root = "../output"
regime = "surface_stokes"
hex_radius = 6000
age_bin_days = 10
# 0 = pool all monthly partitions across years; 1..12 = that month only.
release_month = 0

# w_half members to compare (m/s), comma-separated. Each must have been built
# by 024d; missing members are reported and skipped rather than raising.
w_half_csv = "0.0125,0.025,0.05,0.1,0.2"
# The member treated as the baseline in the narrative summary.
w_half_baseline = 0.05

# Map extent (degrees). The 024a key tiles the whole BSH domain including the
# North Sea, which is empty for Fucus and costs ~35% of panel height; cropping
# to the Baltic proper buys ~1.35x px/km for free. Set all four to 0 to fall
# back to the key's full bounds.
extent_lon_min = 9.0
extent_lon_max = 30.7
extent_lat_min = 53.0
extent_lat_max = 66.0

cmap = "viridis"
panel_height_in = 6
fig_dpi_scale = 3

# %% [markdown]
# # Parse parameters

# %%
output_root = Path(output_root)
w_half_values = [float(x) for x in w_half_csv.split(",") if x]
mpl.rcParams["figure.dpi"] = fig_dpi_scale * mpl.rcParamsDefault["figure.dpi"]
# Hex seam stroke. edgecolor="face" means this is not a visible outline -- it
# closes the ~1 px anti-aliasing seam between adjacent polygons so the grid
# reads as a continuous field. The seam is a fixed PIXEL artifact, so a fixed
# point width makes the resulting hex dilation DPI-invariant (17.5% of hex
# width at every dpi). Pinning it to ~1 px instead lets dilation fall as
# resolution rises: ~8.6% at fig_dpi_scale=3 on the Baltic crop.
hex_seam_lw = 1.1 * 72 / (100 * fig_dpi_scale)

figure_dir = output_root / "Figures" / "031"
figure_dir.mkdir(parents=True, exist_ok=True)

store_root = output_root / "HexAggregates"
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")

month_suffix = f"_m{release_month:02d}" if release_month else ""
month_re = rf"_m{release_month:02d}" if release_month else r"_m\d{2}"


# %% [markdown]
# # Pool each sweep member
#
# One member = every `(year, month)` partition at a given `w_half`. Members
# are additive over release_doy/month/year exactly as a single member is.

# %%
def load_member(w_half):
    wh_suffix = f"_wh{w_half:g}".replace(".", "p")
    part_re = re.compile(
        rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}}){month_re}"
        rf"{re.escape(wh_suffix)}\.parquet$"
    )
    files = [
        f for f in sorted(store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet"))
        if part_re.search(f.name)
    ]
    if not files:
        return None, 0
    return pd.concat(
        [pd.read_parquet(f).reset_index(drop=True) for f in files], ignore_index=True
    ), len(files)


members = {}
for wh in w_half_values:
    df, n = load_member(wh)
    if df is None:
        print(f"w_half={wh:g}: no partitions found — skipped")
        continue
    members[wh] = df
    print(f"w_half={wh:g}: {n} partitions, {len(df):,} rows")
if not members:
    raise FileNotFoundError(
        f"no sweep members found for regime {regime!r} at {store_root} — run 024d."
    )

# %% [markdown]
# # Per-member statistics
#
# `beached_fraction` is stranded weight over all released weight (the residual
# rows carry the never-beached remainder, so the denominator is the full
# release pool). `gini` measures how unevenly the stranded weight is spread
# over the coastal hexes that receive any: 0 = uniform, → 1 = concentrated.


# %%
def gini(values):
    """Gini concentration of a non-negative weight vector."""
    v = np.sort(np.asarray(values, dtype="float64"))
    v = v[v > 0]
    if v.size == 0:
        return np.nan
    n = v.size
    return float((2.0 * np.arange(1, n + 1) - n - 1).dot(v) / (n * v.sum()))


rows = []
for wh, df in members.items():
    beached = df[df["beach_hex"] >= 0]
    per_hex = beached.groupby("beach_hex")["weight"].sum()
    total = float(df["weight"].sum())
    rows.append({
        "w_half": wh,
        "beached_fraction": float(beached["weight"].sum()) / max(total, 1.0),
        "gini": gini(per_hex.to_numpy()),
        "beach_hexes": int(per_hex.size),
        "median_age_days": float(
            (beached.groupby("beach_age_bin")["weight"].sum().pipe(
                lambda s: np.interp(0.5, s.cumsum() / s.sum(), s.index.to_numpy())
            ) + 0.5) * age_bin_days
        ) if len(beached) else np.nan,
    })
stats = pd.DataFrame(rows).set_index("w_half").sort_index()
print(stats.to_string(float_format=lambda v: f"{v:,.4f}"))

# %% [markdown]
# # The range, and what breaks the degeneracy
#
# Left: total beached fraction against `w_half` — the headline number's
# sensitivity. Right: concentration of stranded weight — the statistic that
# distinguishes a wave-selective member from a residence-driven one even where
# totals coincide.

# %%
fig, axes = plt.subplots(1, 2, layout="constrained")
stats["beached_fraction"].plot(ax=axes[0], marker="o", logx=True)
axes[0].set_ylabel("beached fraction of released weight")
axes[0].set_xlabel("w_half (m/s)")
stats["gini"].plot(ax=axes[1], marker="o", logx=True)
axes[1].set_ylabel("Gini concentration of stranded weight")
axes[1].set_xlabel("w_half (m/s)")
fig_path = figure_dir / f"BeachingSweepStats_{regime}_r{hex_radius}m{month_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Where-stranded maps across members
#
# A shared `LogNorm` across every member, so the change is in the pattern and
# not in the colour scale.

# %%
if any((extent_lon_min, extent_lon_max, extent_lat_min, extent_lat_max)):
    lon_min, lat_min, lon_max, lat_max = (
        extent_lon_min, extent_lat_min, extent_lon_max, extent_lat_max
    )
else:
    lon_min, lat_min, lon_max, lat_max = key.total_bounds
extent = [lon_min, lon_max, lat_min, lat_max]
coast = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
).clip(box(lon_min, lat_min, lon_max, lat_max))
domain_aspect = (
    (lon_max - lon_min) * np.cos(np.radians(0.5 * (lat_min + lat_max)))
) / (lat_max - lat_min)


def hex_gdf(df):
    grp = df[df["beach_hex"] >= 0].groupby("beach_hex")["weight"].sum()
    return (
        grp.rename("value").reset_index()
        .merge(key[["hex_id", "geometry"]], left_on="beach_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


gdfs = {wh: hex_gdf(df) for wh, df in members.items()}
shared = pd.concat([g["value"] for g in gdfs.values() if not g.empty])
vmax = float(shared.max())
norm = LogNorm(vmin=max(float(shared[shared > 0].min()), vmax / 1e4), vmax=vmax)

ncols = len(gdfs)
fig, axes = plt.subplots(
    1, ncols, figsize=(panel_height_in * domain_aspect * ncols, panel_height_in),
    layout="constrained", squeeze=False,
)
for ax, (wh, g) in zip(axes[0], sorted(gdfs.items())):
    if not g.empty:
        g.plot(ax=ax, column="value", cmap=cmap, norm=norm, legend=True,
               edgecolor="face", linewidth=hex_seam_lw, zorder=1)
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (extent[2] + extent[3]))))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"w_half = {wh:g} m/s")
fig_path = figure_dir / f"BeachingSweepMaps_{regime}_r{hex_radius}m{month_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Summary

# %%
lo, hi = stats["beached_fraction"].min(), stats["beached_fraction"].max()
print(f"regime={regime}, hex_radius={hex_radius} m, "
      + (f"month={release_month}, " if release_month else "all months, ")
      + f"{len(stats)} sweep members")
print(f"  beached fraction range: {100 * lo:.1f}% .. {100 * hi:.1f}% "
      f"(spread {100 * (hi - lo):.1f} points over w_half "
      f"{stats.index.min():g}..{stats.index.max():g} m/s)")
if w_half_baseline in stats.index:
    print(f"  baseline w_half={w_half_baseline:g}: "
          f"{100 * stats.loc[w_half_baseline, 'beached_fraction']:.1f}%")
print(f"  Gini range: {stats['gini'].min():.3f} .. {stats['gini'].max():.3f} "
      "(higher = stranding concentrated on fewer, wave-exposed hexes)")
