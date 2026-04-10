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

# Quick trajectory plot

Raw trajectories from all completed experiment zarr stores. One color per
experiment type. Subsampled for speed.

```python
import warnings
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)
```

# Parameters

```python
trajectory_dir = Path("/gxfs_work/geomar/smomw122/2025_fucus-dispersal/output/Trajectories")
max_trajs_per_ds = 200
dpi = 300
```

# Load all zarr stores

```python
zarr_files = sorted(trajectory_dir.rglob("*.zarr"))
print(f"Found {len(zarr_files)} zarr stores")

datasets = {}
for zf in zarr_files:
    name = zf.stem
    # Extract experiment type from name: ..._surface_..., ..._bottom_..., ..._surface_stokes_...
    if "surface_stokes" in name:
        exp_type = "surface_stokes"
    elif "surface" in name:
        exp_type = "surface"
    elif "bottom" in name:
        exp_type = "bottom"
    else:
        exp_type = "unknown"

    try:
        ds = xr.open_zarr(zf, decode_times=False)
        # Drop invalid time values (NaT / fill values from oversized chunks)
        # then decode manually
        if "time" in ds:
            raw_time = ds.time.values
            valid_time = np.abs(raw_time) < 1e15  # filter obviously bogus values
            if not valid_time.all():
                print(f"  {name}: filtering {(~valid_time).sum()} invalid time values")
                ds = ds.where(valid_time)
        n_traj = ds.sizes["trajectory"]
        datasets[name] = {"ds": ds, "type": exp_type, "n_traj": n_traj, "path": zf}
        print(f"  {name}: {exp_type}, {n_traj} trajectories")
    except Exception as e:
        print(f"  {name}: SKIPPED ({e})")
```

# Plot

```python
colors = {
    "surface": "tab:blue",
    "bottom": "tab:orange",
    "surface_stokes": "tab:green",
    "unknown": "tab:gray",
}

plt.rcParams["figure.dpi"] = dpi
fig, ax = plt.subplots(figsize=(14, 10))

for name, info in datasets.items():
    ds = info["ds"]
    exp_type = info["type"]
    n = min(info["n_traj"], max_trajs_per_ds)
    color = colors[exp_type]

    # Subsample trajectories
    idx = np.linspace(0, info["n_traj"] - 1, n, dtype=int)
    lon = ds.lon.isel(trajectory=idx).values
    lat = ds.lat.isel(trajectory=idx).values

    for i in range(n):
        valid = ~np.isnan(lon[i])
        if valid.sum() > 1:
            ax.plot(lon[i][valid], lat[i][valid],
                    color=color, linewidth=0.3, alpha=0.3)

# Legend: one entry per type
for exp_type, color in colors.items():
    if any(d["type"] == exp_type for d in datasets.values()):
        ax.plot([], [], color=color, linewidth=1.5, label=exp_type)

ax.legend(loc="upper right", fontsize=10)
ax.set_xlabel("lon")
ax.set_ylabel("lat")
ax.set_aspect(1 / np.cos(np.radians(55)))
ax.set_title(f"All trajectories ({max_trajs_per_ds} per dataset, {len(datasets)} datasets)")
fig.tight_layout()
plt.show()
```
