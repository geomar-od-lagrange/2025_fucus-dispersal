---
jupyter:
  jupytext:
    cell_metadata_filter: tags,-all
    formats: py:percent,md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Build hex connectivity partition

Re-aggregate one `(regime, release_year)` of the 024 counts store into a
HELCOM subbasin→subbasin **residence** connectivity table:
`(origin_subbasin, target_subbasin, release_doy, age_bin) → n_obs`.

This is a pure re-partition of the counts store — every `(trajectory, obs)`
particle-timestep is relabelled from per-hex to per-subbasin and summed.
No trajectory zarrs, no Dask cluster: reads parquet, writes parquet.

Hard prerequisites: the counts partition and the hex key + sidecar built
by 024a/024 for the matching `(regime, release_year, hex_radius)`.

```python
import json
from pathlib import Path

import pandas as pd
```

# Parameters

```python tags=["parameters"]
# Read root of the hex-aggregate store and write root for the connectivity partition.
output_root = "../output"

# One (regime, release_year) per run.
regime = "surface"
release_year = 2019

# Hex radius (must match an existing key + counts file).
hex_radius = 6000

# Age-bin granularity — must match the counts file being read.
# 024c groups the already-materialised age_bin column from counts; it does
# not re-bin. Kept here only for cross-notebook param-block consistency.
age_bin_days = 10
```

# Derived layout / resolve prerequisites

```python
output_root = Path(output_root)
store_root = output_root / "HexAggregates"
store_root.mkdir(parents=True, exist_ok=True)

key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
meta_path = key_path.with_suffix(".json")
counts_path = store_root / f"HexAgg_counts_r{hex_radius}m_{regime}_{release_year}.parquet"
connectivity_path = (
    store_root / f"HexAgg_connectivity_r{hex_radius}m_{regime}_{release_year}.parquet"
)

missing = [p for p in (key_path, meta_path, counts_path) if not p.exists()]
if missing:
    raise FileNotFoundError(
        "Required inputs missing — run 024a_BuildHexKey and 024_BuildHexAggregates first.\n"
        + "".join(f"  expected: {p}\n" for p in missing)
    )

print(f"counts   ← {counts_path}")
print(f"key      ← {key_path}")
print(f"sidecar  ← {meta_path}")
print(f"output   → {connectivity_path}")
```

# Read counts + key

```python
counts = pd.read_parquet(counts_path)
print(f"counts: {len(counts):,} rows")

key = pd.read_parquet(key_path, columns=["hex_id", "helcom_subbasin"])
subbasin_id_to_name = {
    int(k): v
    for k, v in json.loads(meta_path.read_text())["subbasin_id_to_name"].items()
}
print(f"key: {len(key):,} hexes, {len(subbasin_id_to_name)} named subbasins")
```

# Map release_hex / target_hex → subbasin

`release_hex` and `target_hex` are mapped to `helcom_subbasin` via the key.
Two unnamed states both collapse to the single sentinel -1:
- hexes absent from the key (land-seed / NaN, `hex_id == -1`) → lookup
  yields NaN → filled to -1.
- in-key hexes outside every named polygon (`helcom_subbasin == -1`,
  the `_outside` category) → already -1, unchanged.

Filter with `>= 0` at read time when only named subbasins are wanted.

```python
hex_to_subbasin = key.set_index("hex_id")["helcom_subbasin"]

counts["origin_subbasin"] = (
    counts["release_hex"].map(hex_to_subbasin).fillna(-1).astype(int)
)
counts["target_subbasin"] = (
    counts["target_hex"].map(hex_to_subbasin).fillna(-1).astype(int)
)
```

# Group → connectivity partition

```python
conn = (
    counts
    .groupby(["origin_subbasin", "target_subbasin", "release_doy", "age_bin"])["n_obs"]
    .sum()
    .reset_index()
)
print(f"connectivity: {len(conn):,} rows")
```

# Write parquet

```python
conn.to_parquet(connectivity_path)
print(f"wrote {connectivity_path} ({connectivity_path.stat().st_size / 1e6:.2f} MB)")
```

# Validation

```python
# Sum-preserved invariant: every (trajectory, obs) timestep is re-labelled,
# not dropped or duplicated. Each hex maps to exactly one subbasin (centroid
# rule), so the groupby is a partition with no overlap.
assert conn["n_obs"].sum() == counts["n_obs"].sum(), (
    f"n_obs sum mismatch: conn {conn['n_obs'].sum()} != counts {counts['n_obs'].sum()}"
)
print("OK: n_obs sum preserved (re-partition invariant holds)")
```

```python
total_nobs = conn["n_obs"].sum()

named_origins = sorted(conn.loc[conn["origin_subbasin"] >= 0, "origin_subbasin"].unique())
named_targets = sorted(conn.loc[conn["target_subbasin"] >= 0, "target_subbasin"].unique())

frac_unnamed_origin = (
    conn.loc[conn["origin_subbasin"] < 0, "n_obs"].sum() / total_nobs
)
frac_unnamed_target = (
    conn.loc[conn["target_subbasin"] < 0, "n_obs"].sum() / total_nobs
)

diag_mask = conn["origin_subbasin"] == conn["target_subbasin"]
frac_diagonal = conn.loc[diag_mask & (conn["origin_subbasin"] >= 0), "n_obs"].sum() / total_nobs

print(f"regime={regime}, release_year={release_year}, hex_radius={hex_radius} m")
print(f"  connectivity rows:       {len(conn):,}")
print(f"  named origin subbasins:  {len(named_origins)}")
print(f"  named target subbasins:  {len(named_targets)}")
print(f"  n_obs fraction in -1 origin:  {frac_unnamed_origin:.3f}")
print(f"  n_obs fraction in -1 target:  {frac_unnamed_target:.3f}")
print(f"  within-subbasin diagonal share (origin==target, named): {frac_diagonal:.3f}")
```
