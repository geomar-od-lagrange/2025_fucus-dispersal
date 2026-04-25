---
jupyter:
  jupytext:
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

# Time statistics

Global ensemble diagnostics across all regimes. Land-seed filter count
and final displacement distribution.

```python
import os
import time

import dask
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from dask.distributed import Client
from pathlib import Path
```

# Parameters

```python tags=["parameters"]
# Read root of trajectory zarrs (output_root/Trajectories/<regime>/...).
output_root = "../output"
```

# Dask cluster

Connect to an external scheduler when ``SCHEDULER_FILE`` is set (written
by the multi-task SLURM job). Otherwise spin up a local cluster on the
current node.

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

# Load all regimes

Layout assumption: ``output_root/Trajectories/<regime>/<release_year>/*.zarr``.
``regimes`` is the list of immediate subdirectories.

```python
output_root = Path(output_root)
trajectory_root = output_root / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")
```

```python
regime_datasets = {}
land_seeded_masks = {}
for regime in regimes:
    zarr_files = sorted((trajectory_root / regime).glob("**/*.zarr"))
    ds_raw = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")
    # Trajectories whose first step has zero displacement were seeded on
    # land and never advected; mask them across all obs.
    land_seeded = (
        (ds_raw.lon.diff("obs").isel(obs=0, drop=True) == 0)
        & (ds_raw.lat.diff("obs").isel(obs=0, drop=True) == 0)
    )
    ds_masked = ds_raw.where(~land_seeded)
    regime_datasets[regime] = ds_masked
    land_seeded_masks[regime] = land_seeded
regime_datasets
```

# Compute per-regime diagnostics (one shared pass per regime)

For each regime, land-seed count and first/last coords walk the same
trajectory graph; `dask.compute(*)` evaluates them in one go.

Last valid obs uses `ffill("obs")` + `isel(obs=-1)` to tolerate any
NaN-tail inside a chunk (trajectories that die before their chunk ends).

```python
per_regime = {}
for regime, ds in regime_datasets.items():
    n_total = ds.sizes["trajectory"]
    land_count_lazy = land_seeded_masks[regime].sum()
    first_lon_lazy = ds.lon.isel(obs=0, drop=True)
    first_lat_lazy = ds.lat.isel(obs=0, drop=True)
    last_lon_lazy = ds.lon.ffill("obs").isel(obs=-1, drop=True)
    last_lat_lazy = ds.lat.ffill("obs").isel(obs=-1, drop=True)
    land_count, first_lon, first_lat, last_lon, last_lat = dask.compute(
        land_count_lazy, first_lon_lazy, first_lat_lazy, last_lon_lazy, last_lat_lazy,
    )
    land = int(land_count)
    per_regime[regime] = dict(
        n_total=n_total,
        n_valid=n_total - land,
        first_lon=first_lon,
        first_lat=first_lat,
        last_lon=last_lon,
        last_lat=last_lat,
    )
```

# Per-regime summary

```python
print("\n".join(
    f"{r}: {per_regime[r]['n_total']} trajectories, "
    f"{per_regime[r]['n_total'] - per_regime[r]['n_valid']} land-seeded "
    f"({per_regime[r]['n_valid']} valid)"
    for r in regimes
))
```

# Land-seed count per regime

```python
land_seed = pd.Series(
    {r: int(per_regime[r]["n_total"] - per_regime[r]["n_valid"]) for r in regimes},
    name="land_seeded",
)
land_seed.plot.bar()
```

# Final displacement distribution

Great-circle approximation (111 km per degree lat). NaN last-lon (land-
seeded) drops out of the histogram automatically. `histtype="step"` so
the three regimes overlay without bar occlusion.

```python
finals = {}
for regime in regimes:
    r = per_regime[regime]
    dlat = r["last_lat"] - r["first_lat"]
    dlon = (r["last_lon"] - r["first_lon"]) * np.cos(np.deg2rad(r["first_lat"]))
    finals[regime] = (111.0 * np.sqrt(dlat ** 2 + dlon ** 2)).values

df_final = pd.DataFrame({r: pd.Series(v) for r, v in finals.items()}).dropna(how="all")
ax = df_final.plot.hist(bins=50, histtype="step")
ax.set_xlabel("Final displacement (km)")
```
