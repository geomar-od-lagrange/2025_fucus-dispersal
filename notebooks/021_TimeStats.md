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

# Time statistics

Global ensemble diagnostics across all regimes. Land-seed filter count,
trajectory lifetimes, alive fraction vs age, final displacement
distribution.

```python
import dask
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

from helpers import load_trajectories, mask_land_seeded
```

# Parameters

```python tags=["parameters"]
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"

panel_size = 4
```

# Dask cluster

```python
from dask.distributed import Client
client = Client(ip="0.0.0.0")
client
```

# Load all regimes

```python
base_path = Path(base_path)
trajectory_root = base_path / "output" / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")
```

```python
regime_datasets = {}
land_seeded_masks = {}
for regime in regimes:
    ds_raw, _ = load_trajectories(trajectory_root / regime)
    ds_masked, land_seeded = mask_land_seeded(ds_raw)
    regime_datasets[regime] = ds_masked
    land_seeded_masks[regime] = land_seeded
regime_datasets
```

# Compute per-regime diagnostics (one shared pass per regime)

For each regime, lifetime / alive-fraction / first-last coords all walk
the same trajectory graph; `dask.compute(*)` evaluates them in one go.

Last valid obs uses `ffill("obs")` + `isel(obs=-1)` to tolerate any
NaN-tail inside a chunk (trajectories that die before their chunk ends).

```python
per_regime = {}
for regime, ds in regime_datasets.items():
    n_total = ds.sizes["trajectory"]
    lazy = dict(
        land_count=land_seeded_masks[regime].sum(),
        lifetime=ds.lon.notnull().sum("obs"),
        alive_numerator=ds.lon.notnull().sum("trajectory"),
        first_lon=ds.lon.isel(obs=0, drop=True),
        first_lat=ds.lat.isel(obs=0, drop=True),
        last_lon=ds.lon.ffill("obs").isel(obs=-1, drop=True),
        last_lat=ds.lat.ffill("obs").isel(obs=-1, drop=True),
    )
    results = dict(zip(lazy.keys(), dask.compute(*lazy.values())))
    land = int(results["land_count"])
    results["n_total"] = n_total
    results["n_valid"] = n_total - land
    results["alive"] = results["alive_numerator"] / results["n_valid"]
    per_regime[regime] = results
    print(
        f"{regime}: {n_total} trajectories, {land} land-seeded "
        f"({results['n_valid']} valid)"
    )
```

# Land-seed count per regime

```python
land_seed = pd.Series(
    {r: int(per_regime[r]["n_total"] - per_regime[r]["n_valid"]) for r in regimes},
    name="land_seeded",
)
land_seed.plot.bar(figsize=(panel_size, panel_size))
```

# Lifetime distribution

Valid-obs count per trajectory; land-seeded (all-NaN) contribute 0 and
are dropped from the histogram.

```python
df_life = pd.DataFrame({
    regime: pd.Series(per_regime[regime]["lifetime"].values)
    for regime in regimes
})
df_life = df_life.where(df_life > 0)
df_life.plot.hist(bins=50, alpha=0.5, figsize=(panel_size * 2, panel_size))
```

# Alive fraction vs age

```python
da_alive = xr.concat(
    [per_regime[r]["alive"].expand_dims(regime=[r]) for r in regimes],
    dim="regime",
)
da_alive.to_pandas().T.plot.line(figsize=(panel_size * 2, panel_size))
```

# Final displacement distribution

Great-circle approximation (111 km per degree lat). NaN last-lon (land-
seeded) drops out of the histogram automatically.

```python
finals = {}
for regime in regimes:
    r = per_regime[regime]
    dlat = r["last_lat"] - r["first_lat"]
    dlon = (r["last_lon"] - r["first_lon"]) * np.cos(np.deg2rad(r["first_lat"]))
    finals[regime] = (111.0 * np.sqrt(dlat ** 2 + dlon ** 2)).values

df_final = pd.DataFrame({r: pd.Series(v) for r, v in finals.items()}).dropna(how="all")
df_final.plot.hist(bins=50, alpha=0.5, figsize=(panel_size * 2, panel_size))
```
