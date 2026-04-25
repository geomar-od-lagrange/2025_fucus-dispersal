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
    display_name: min_data (pixi)
    language: python
    name: min_data
---

<!-- #region papermill={"duration": 0.006373, "end_time": "2026-04-25T12:15:30.059363+00:00", "exception": false, "start_time": "2026-04-25T12:15:30.052990+00:00", "status": "completed"} -->
# Time statistics

Global ensemble diagnostics across all regimes. Land-seed filter count
and final displacement distribution.
<!-- #endregion -->

```python papermill={"duration": 0.486637, "end_time": "2026-04-25T12:15:30.550142+00:00", "exception": false, "start_time": "2026-04-25T12:15:30.063505+00:00", "status": "completed"}
import dask
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

import xarray as xr
```

<!-- #region papermill={"duration": 0.000762, "end_time": "2026-04-25T12:15:30.551940+00:00", "exception": false, "start_time": "2026-04-25T12:15:30.551178+00:00", "status": "completed"} -->
# Parameters
<!-- #endregion -->

```python papermill={"duration": 0.003395, "end_time": "2026-04-25T12:15:30.556062+00:00", "exception": false, "start_time": "2026-04-25T12:15:30.552667+00:00", "status": "completed"} tags=["parameters"]
output_root = "../output"
```

<!-- #region papermill={"duration": 0.000769, "end_time": "2026-04-25T12:15:30.557681+00:00", "exception": false, "start_time": "2026-04-25T12:15:30.556912+00:00", "status": "completed"} -->
# Dask cluster

Connect to an external scheduler when ``SCHEDULER_FILE`` is set (written
by the multi-task SLURM job). Otherwise spin up a local cluster on the
current node.
<!-- #endregion -->

```python papermill={"duration": 0.571013, "end_time": "2026-04-25T12:15:31.129430+00:00", "exception": false, "start_time": "2026-04-25T12:15:30.558417+00:00", "status": "completed"}
import os
import time
from dask.distributed import Client

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

<!-- #region papermill={"duration": 0.000827, "end_time": "2026-04-25T12:15:31.131372+00:00", "exception": false, "start_time": "2026-04-25T12:15:31.130545+00:00", "status": "completed"} -->
# Load all regimes
<!-- #endregion -->

```python papermill={"duration": 0.004188, "end_time": "2026-04-25T12:15:31.136372+00:00", "exception": false, "start_time": "2026-04-25T12:15:31.132184+00:00", "status": "completed"}
output_root = Path(output_root)
trajectory_root = output_root / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")
```

```python papermill={"duration": 0.452354, "end_time": "2026-04-25T12:15:31.589940+00:00", "exception": false, "start_time": "2026-04-25T12:15:31.137586+00:00", "status": "completed"}
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

<!-- #region papermill={"duration": 0.000889, "end_time": "2026-04-25T12:15:31.591966+00:00", "exception": false, "start_time": "2026-04-25T12:15:31.591077+00:00", "status": "completed"} -->
# Compute per-regime diagnostics (one shared pass per regime)

For each regime, land-seed count and first/last coords walk the same
trajectory graph; `dask.compute(*)` evaluates them in one go.

Last valid obs uses `ffill("obs")` + `isel(obs=-1)` to tolerate any
NaN-tail inside a chunk (trajectories that die before their chunk ends).
<!-- #endregion -->

```python papermill={"duration": 9.284953, "end_time": "2026-04-25T12:15:40.877747+00:00", "exception": false, "start_time": "2026-04-25T12:15:31.592794+00:00", "status": "completed"}
per_regime = {}
for regime, ds in regime_datasets.items():
    n_total = ds.sizes["trajectory"]
    lazy = dict(
        land_count=land_seeded_masks[regime].sum(),
        first_lon=ds.lon.isel(obs=0, drop=True),
        first_lat=ds.lat.isel(obs=0, drop=True),
        last_lon=ds.lon.ffill("obs").isel(obs=-1, drop=True),
        last_lat=ds.lat.ffill("obs").isel(obs=-1, drop=True),
    )
    results = dict(zip(lazy.keys(), dask.compute(*lazy.values())))
    land = int(results["land_count"])
    results["n_total"] = n_total
    results["n_valid"] = n_total - land
    per_regime[regime] = results
    print(
        f"{regime}: {n_total} trajectories, {land} land-seeded "
        f"({results['n_valid']} valid)"
    )
```

<!-- #region papermill={"duration": 0.001021, "end_time": "2026-04-25T12:15:40.880131+00:00", "exception": false, "start_time": "2026-04-25T12:15:40.879110+00:00", "status": "completed"} -->
# Land-seed count per regime
<!-- #endregion -->

```python papermill={"duration": 0.058603, "end_time": "2026-04-25T12:15:40.939681+00:00", "exception": false, "start_time": "2026-04-25T12:15:40.881078+00:00", "status": "completed"}
land_seed = pd.Series(
    {r: int(per_regime[r]["n_total"] - per_regime[r]["n_valid"]) for r in regimes},
    name="land_seeded",
)
land_seed.plot.bar()
```

<!-- #region papermill={"duration": 0.001039, "end_time": "2026-04-25T12:15:40.941890+00:00", "exception": false, "start_time": "2026-04-25T12:15:40.940851+00:00", "status": "completed"} -->
# Final displacement distribution

Great-circle approximation (111 km per degree lat). NaN last-lon (land-
seeded) drops out of the histogram automatically. `histtype="step"` so
the three regimes overlay without bar occlusion.
<!-- #endregion -->

```python papermill={"duration": 0.075999, "end_time": "2026-04-25T12:15:41.018891+00:00", "exception": false, "start_time": "2026-04-25T12:15:40.942892+00:00", "status": "completed"}
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
