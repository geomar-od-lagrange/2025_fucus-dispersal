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
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
```

# Parameters

```python tags=["parameters"]
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"
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
land_seed_counts = {}
for regime in regimes:
    zarr_files = sorted((trajectory_root / regime).glob("**/*.zarr"))
    ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")
    dlon0 = ds.lon.diff("obs").isel(obs=0)
    dlat0 = ds.lat.diff("obs").isel(obs=0)
    land = (
        ((dlon0 == 0) & (dlat0 == 0)).drop_vars("obs", errors="ignore").compute()
    )
    land_seed_counts[regime] = int(land.sum())
    regime_datasets[regime] = ds.isel(trajectory=~land)
    print(f"{regime}: {ds.sizes['trajectory']} trajectories, {land_seed_counts[regime]} land-seeded")
```

# Dask cluster

```python
from dask.distributed import Client
client = Client()
client
```

# Land-seed count per regime

```python
land_seed = pd.Series(land_seed_counts, name="land_seeded")
land_seed.plot.bar()
```

# Lifetime distribution

Max valid obs index per trajectory.

```python
lifetimes = {
    regime: ds.lon.notnull().sum("obs").compute().values
    for regime, ds in regime_datasets.items()
}
df_life = pd.DataFrame({r: pd.Series(v) for r, v in lifetimes.items()})
df_life.plot.hist(bins=50, alpha=0.5)
```

# Alive fraction vs age

```python
alive_curves = {
    regime: ds.lon.notnull().mean("trajectory").compute()
    for regime, ds in regime_datasets.items()
}
da_alive = xr.concat(
    [v.expand_dims(regime=[k]) for k, v in alive_curves.items()], dim="regime"
)
da_alive.plot.line(x="obs", hue="regime")
```

# Final displacement distribution

Great-circle approximation (111 km per degree lat).

```python
finals = {}
for regime, ds in regime_datasets.items():
    last_valid = ds.lon.notnull().sum("obs").compute() - 1
    last_lon = ds.lon.isel(obs=last_valid).compute()
    last_lat = ds.lat.isel(obs=last_valid).compute()
    lon0 = ds.lon.isel(obs=0).compute()
    lat0 = ds.lat.isel(obs=0).compute()
    dlat = last_lat - lat0
    dlon = (last_lon - lon0) * np.cos(np.deg2rad(lat0))
    finals[regime] = 111.0 * np.sqrt(dlat ** 2 + dlon ** 2).values

df_final = pd.DataFrame({r: pd.Series(v) for r, v in finals.items()})
df_final.plot.hist(bins=50, alpha=0.5)
```
