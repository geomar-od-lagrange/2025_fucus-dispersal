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
# # Subbasin connectivity matrix (POC)
#
# Proof-of-concept consumer of the 024c connectivity store. Reads
# ``HexAgg_connectivity_r<radius>m_<regime>_<year>.parquet`` files, pools
# across all available years, and for two release-month scopes (all-year and
# Aug/Sep) prints the subbasin→subbasin residence matrix and draws linear and
# log-scale heatmaps.
#
# Connectivity here is **residence** (particle-timesteps in the target
# subbasin), not particle flux — see ``plans/subbasin_connectivity.md`` for
# the semantics. The ``origin_subbasin``/``target_subbasin`` ids come from the
# HELCOM centroid assignment baked into 024a; id -1 = unnamed/outside and is
# excluded from the matrix (reported as a dropped fraction). No Dask cluster;
# parquet-only.

# %%
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the hex-aggregate store.
output_root = "../output"
# Which regime's connectivity partitions to read; one regime per run.
regime = "surface"
# Hex radius of the store (built by 024a/024c). Must match files on disk.
hex_radius = 6000
# Age-bin granularity of the connectivity store (must match 024). Pooled away
# here (we sum all age bins), but kept so the parameters cell matches sibling
# notebooks (026, 027) and is readable as the store-contract value.
age_bin_days = 10
# Figure DPI as a multiple of the matplotlib default. >1 sharpens the
# inline/embedded raster panels — the one plotting default we override here
# (rationale in docs/visualisations.md).
fig_dpi_scale = 2

# %% [markdown]
# # Parse parameters

# %%
output_root = Path(output_root)
# Set from the stock default so re-running this cell is idempotent.
mpl.rcParams["figure.dpi"] = fig_dpi_scale * mpl.rcParamsDefault["figure.dpi"]

# PNG figures land under the outputs tree (outside the repo), one subdir per
# notebook stage. savefig inherits figure.dpi (set above), so saved panels are
# as sharp as the inline ones — no dpi= kwarg needed.
figure_dir = output_root / "Figures" / "028"
figure_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Release-month scopes
#
# Derived here (not in the parameters cell — the plan requires primitives only
# in the parameters cell). ``[]`` means all months; ``[8, 9]`` keeps
# August/September only.

# %%
scopes = {
    "all_year": [],
    "aug_sep": [8, 9],
}

# %% [markdown]
# # Read subbasin id→name from the key sidecar + pool connectivity across years
#
# Layout: ``output_root/HexAggregates/HexAgg_key_r<radius>m.json`` carries the
# ``subbasin_id_to_name`` map (JSON string keys → cast to int). The connectivity
# partitions are ``HexAgg_connectivity_r<radius>m_<regime>_<year>.parquet``.
# The release year is parsed from each filename so the
# ``release_doy → month`` conversion is leap-correct (same pattern as 026a).

# %%
store_root = output_root / "HexAggregates"
key_json_path = store_root / f"HexAgg_key_r{hex_radius}m.json"
subbasin_id_to_name = {
    int(k): v
    for k, v in json.loads(key_json_path.read_text())["subbasin_id_to_name"].items()
}

# Match exactly four year digits after the regime so `surface` does not also
# glob `surface_stokes_*` files (prefix collision).
_CONN_RE = re.compile(rf"HexAgg_connectivity_r{hex_radius}m_{regime}_(\d{{4}})\.parquet$")
conn_files = sorted(
    store_root.glob(
        f"HexAgg_connectivity_r{hex_radius}m_{regime}_[0-9][0-9][0-9][0-9].parquet"
    )
)
if not conn_files:
    raise FileNotFoundError(
        f"no connectivity partitions for regime {regime!r} at {store_root} — run 024c."
    )

parts = []
for f in conn_files:
    year = int(_CONN_RE.search(f.name).group(1))
    df = pd.read_parquet(f).reset_index(drop=True)
    df["release_month"] = pd.to_datetime(
        (year * 1000 + df["release_doy"].astype("int32")).astype(str), format="%Y%j"
    ).dt.month
    parts.append(df)
conn = pd.concat(parts, ignore_index=True)
print(f"pooled {len(conn_files)} year partition(s); {len(conn):,} rows total")

# %% [markdown]
# # Per-scope: matrix + heatmaps

# %%
for scope_name, months in scopes.items():
    # --- filter to scope months -------------------------------------------------
    if months:
        scope_df = conn[conn["release_month"].isin(months)].copy()
        print(f"\nscope={scope_name!r} (months {months}): {len(scope_df):,} rows")
    else:
        scope_df = conn.copy()
        print(f"\nscope={scope_name!r} (all months): {len(scope_df):,} rows")

    # --- pool over release_doy and age_bin → (origin_subbasin, target_subbasin) ---
    pooled = (
        scope_df.groupby(["origin_subbasin", "target_subbasin"])["n_obs"]
        .sum()
        .reset_index()
    )

    # --- report -1 (unnamed/outside) fraction dropped ---------------------------
    total_obs = int(pooled["n_obs"].sum())
    unnamed_mask = (pooled["origin_subbasin"] < 0) | (pooled["target_subbasin"] < 0)
    unnamed_obs = int(pooled.loc[unnamed_mask, "n_obs"].sum())
    frac_dropped = unnamed_obs / total_obs if total_obs > 0 else float("nan")
    print(
        f"  n_obs in -1 origin/target: {unnamed_obs:,} / {total_obs:,} "
        f"({frac_dropped:.1%} dropped from matrix)"
    )

    # restrict to named subbasins
    named = pooled[(pooled["origin_subbasin"] >= 0) & (pooled["target_subbasin"] >= 0)].copy()

    # --- build square matrix ordered by id ascending ----------------------------
    all_ids = sorted(
        set(named["origin_subbasin"].unique()) | set(named["target_subbasin"].unique())
    )
    if not all_ids:
        print(f"  no named-subbasin connectivity in scope {scope_name!r} — skipping matrix/plots")
        continue
    id_names = [subbasin_id_to_name[i] for i in all_ids]

    pivot = named.pivot(
        index="origin_subbasin", columns="target_subbasin", values="n_obs"
    ).reindex(index=all_ids, columns=all_ids).fillna(0)
    # Rename index and columns to subbasin names for display.
    matrix_df = pivot.copy()
    matrix_df.index = pd.Index(id_names, name="origin")
    matrix_df.columns = pd.Index(id_names, name="target")

    print(f"  matrix shape: {matrix_df.shape}")
    print(matrix_df.to_string())

    # --- heatmap: linear and log ------------------------------------------------
    matrix_values = matrix_df.to_numpy(dtype=float)

    for scale in ("lin", "log"):
        fig, ax = plt.subplots(layout="constrained")

        if scale == "log":
            # LogNorm is required for log-scale display (not a cosmetic override):
            # the diagonal self-connectivity dominates by orders of magnitude,
            # so a linear scale washes out off-diagonal structure. Zeros are
            # masked to NaN so empty cells render blank (imshow treats NaN as
            # transparent / uses the axes facecolor).
            masked = np.where(matrix_values > 0, matrix_values, np.nan)
            im = ax.imshow(masked, norm=LogNorm())
        else:
            im = ax.imshow(matrix_values)

        fig.colorbar(im, ax=ax)

        # Tick labels set at plot time only (AGENTS.md: don't mutate coords
        # or DataFrame index for presentation purposes).
        n = len(id_names)
        ax.set_xticks(range(n))
        ax.set_xticklabels(id_names, rotation=90, ha="right")
        ax.set_yticks(range(n))
        ax.set_yticklabels(id_names)
        ax.set_xlabel("target subbasin")
        ax.set_ylabel("origin subbasin")
        ax.set_title(f"connectivity — {scope_name} — {scale}")

        fig_path = (
            figure_dir
            / f"SubbasinConnectivityMatrix_{regime}_r{hex_radius}m_{scope_name}_{scale}.png"
        )
        fig.savefig(fig_path)
        print(f"  wrote {fig_path}")
        plt.show()

# %% [markdown]
# # Validation / summary

# %%
print(f"regime={regime}, hex_radius={hex_radius} m, age_bin_days={age_bin_days}")
print(f"years pooled: {len(conn_files)}")
for scope_name, months in scopes.items():
    if months:
        scope_df = conn[conn["release_month"].isin(months)]
    else:
        scope_df = conn

    pooled = (
        scope_df.groupby(["origin_subbasin", "target_subbasin"])["n_obs"]
        .sum()
        .reset_index()
    )
    total_obs = int(pooled["n_obs"].sum())
    unnamed_obs = int(
        pooled.loc[
            (pooled["origin_subbasin"] < 0) | (pooled["target_subbasin"] < 0), "n_obs"
        ].sum()
    )
    named = pooled[(pooled["origin_subbasin"] >= 0) & (pooled["target_subbasin"] >= 0)]
    all_ids = sorted(
        set(named["origin_subbasin"].unique()) | set(named["target_subbasin"].unique())
    )
    matrix_obs = int(named["n_obs"].sum())

    # Diagonal (self-retention): origin_subbasin == target_subbasin rows.
    diag_obs = int(
        named.loc[named["origin_subbasin"] == named["target_subbasin"], "n_obs"].sum()
    )
    diag_share = diag_obs / matrix_obs if matrix_obs > 0 else float("nan")
    frac_dropped = unnamed_obs / total_obs if total_obs > 0 else float("nan")

    print(
        f"\nscope={scope_name!r}: "
        f"matrix {len(all_ids)}×{len(all_ids)}, "
        f"total n_obs={total_obs:,}, "
        f"diagonal (retention) share={diag_share:.1%}, "
        f"dropped -1 fraction={frac_dropped:.1%}"
    )
