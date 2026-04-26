---
jupyter:
  jupytext:
    cell_metadata_filter: tags,-all
    formats: md,ipynb
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

# Build hex-aggregate counts partition

Aggregate one `(regime, release_year)` worth of trajectory zarrs into a
flat `(release_hex, age_bin, target_hex, release_doy) → n_obs` table.

The key file from `024a_BuildHexKey.md` is a hard prerequisite — its
sidecar JSON carries the `HexProj` parameters used here so labels match
the key one-for-one.

Same `(release_time, regime)` zarrs are aggregated additively (multiple
seeds = independent reseeded reruns, summed via groupby). Missing
release-date zarrs are tolerated: the build globs whatever is on disk.

```python
import json
import os
import re
import time
from pathlib import Path

import dask.dataframe as dd
import pandas as pd
import xarray as xr
from dask.distributed import Client

from hextraj import HexProj
from hextraj.hex_analysis import hex_connectivity_dask
```

```python
# Pattern: Fucus_BSH_YYYYMMDDTHHMMSS_{regime}_dt{N}min_seed{S}.
# `surface_stokes` must precede `surface` so the alternation matches the
# longer form first.
_ZARR_STEM_RE = re.compile(
    r"^Fucus_BSH_(\d{8}T\d{6})_(surface_stokes|surface|bottom)_dt\d+min_seed\d+$"
)


def parse_zarr_stem(path):
    """Parse a trajectory zarr filename into ``(release_time, regime)``."""
    m = _ZARR_STEM_RE.match(Path(path).stem)
    if m is None:
        raise ValueError(
            f"zarr filename does not match expected pattern: {Path(path).name!r}"
        )
    return pd.Timestamp(m.group(1)), m.group(2)
```

# Parameters

```python tags=["parameters"]
# Read root of trajectory zarrs and write root for the counts partition.
output_root = "../output"

# One (regime, release_year) per run.
regime = "surface"
release_year = 2019

# Hex radius (must match an existing key file built by 024a).
hex_radius = 6000

# Age-bin granularity. The set of bins that show up emerges from the
# zarrs — no upper cap; downstream consumers query by age at read time.
age_bin_days = 10
# Zarr output cadence.
output_dt_mins = 60
```

# Derived layout / key + projection

```python
output_root = Path(output_root)
store_root = output_root / "HexAggregates"
store_root.mkdir(parents=True, exist_ok=True)

key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
meta_path = key_path.with_suffix(".json")
counts_path = (
    store_root / f"HexAgg_counts_r{hex_radius}m_{regime}_{release_year}.parquet"
)

if not key_path.exists() or not meta_path.exists():
    raise FileNotFoundError(
        f"Key file or sidecar missing — run 024a_BuildHexKey.md first.\n"
        f"  expected: {key_path}\n  expected: {meta_path}"
    )

meta = json.loads(meta_path.read_text())
hp = HexProj(**meta["hex_proj"])
print(f"HexProj: {meta['hex_proj']}")
print(f"counts → {counts_path}")
```

# Dask cluster

```python
scheduler_file = os.environ.get("SCHEDULER_FILE")
if scheduler_file:
    for _ in range(60):
        if os.path.exists(scheduler_file):
            break
        time.sleep(1)
    client = Client(scheduler_file=scheduler_file)
else:
    client = Client(ip="0.0.0.0")
client
```

# Trajectory zarrs → counts

```python
zarrs = sorted(
    (output_root / f"Trajectories/{regime}/{release_year}").glob("*.zarr")
)
parsed = [(p, *parse_zarr_stem(p)) for p in zarrs]
for p, ts, fn_regime in parsed:
    assert ts.year == release_year, (ts, release_year, p)
    assert fn_regime == regime, (fn_regime, regime, p)

if not parsed:
    raise FileNotFoundError(
        f"no zarrs at {output_root}/Trajectories/{regime}/{release_year}/"
    )
release_doys = sorted({int(ts.dayofyear) for _, ts, _ in parsed})
print(f"{len(parsed)} zarrs, release_doys "
      f"{release_doys[0]}..{release_doys[-1]} ({len(release_doys)} unique)")
```

```python
def zarr_to_lazy_frame(path, release_doy):
    """Lazy (release_hex, target_hex, age_bin, release_doy) frame for one zarr.
    Land-seeded particles (zero first-step displacement) NaN out lon/lat and
    surface as INVALID_HEX_ID (-1); kept as a sentinel in both hex columns."""
    ds = xr.open_zarr(path)
    seeded_on_land = (
        (ds.lon.diff("obs").isel(obs=0, drop=True) == 0)
        & (ds.lat.diff("obs").isel(obs=0, drop=True) == 0)
    )
    ds = ds[["lon", "lat"]].where(~seeded_on_land)

    ddf = hex_connectivity_dask(ds, hp, traj_dim="trajectory").rename(
        columns={"from_id": "release_hex", "to_id": "target_hex"}
    )
    ddf["age_bin"] = ddf["obs"] * output_dt_mins // (60 * 24) // age_bin_days
    ddf["release_doy"] = release_doy
    return ddf[["release_hex", "release_doy", "age_bin", "target_hex"]]


t0 = time.time()
counts = (
    dd.concat([zarr_to_lazy_frame(p, ts.dayofyear) for p, ts, _ in parsed])
    .groupby(["release_hex", "release_doy", "age_bin", "target_hex"])
    .size().rename("n_obs").reset_index()
    .compute()
)
print(f"computed {len(counts):,} rows in {time.time() - t0:.1f}s")
```

```python
counts.to_parquet(counts_path)
print(f"wrote {counts_path} ({counts_path.stat().st_size / 1e6:.2f} MB)")
```

# Validation

```python
key_ids = set(pd.read_parquet(key_path, columns=["hex_id"])["hex_id"].astype(int))
seen = (
    set(counts["release_hex"].astype(int))
    | set(counts["target_hex"].astype(int))
)
unseen = seen - key_ids - {-1}  # -1 is the INVALID_HEX_ID sentinel.
assert not unseen, f"hex_ids missing from key.parquet: {sorted(unseen)[:10]} ..."
print(f"PASS: every release_hex/target_hex is in {key_path.name} (or -1).")
```

```python
print(f"regime={regime}, release_year={release_year}, hex_radius={hex_radius} m")
print(f"  rows:           {len(counts):,}")
print(f"  sum(n_obs):     {int(counts['n_obs'].sum()):,}")
print(f"  release_doys:   {counts['release_doy'].nunique()} "
      f"({counts['release_doy'].min()}..{counts['release_doy'].max()})")
print(f"  unique target:  {counts['target_hex'].nunique():,}")
print(f"  age_bins used:  {counts['age_bin'].min()}..{counts['age_bin'].max()}")
```
