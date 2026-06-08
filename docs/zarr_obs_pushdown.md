# Reading release-time metadata: push the obs slice into the zarr read

Trajectory zarrs are chunked `(trajectory=10000, obs=1000)` (see
[seeding.md](seeding.md) / `notebooks/010_FucusDispersal`). The
post-processing notebooks need per-trajectory **release-time** quantities —
release subbasin, release quarter, the on-land test, release position — all
derived from `obs=0` (the land test also needs `obs=1`).

The obvious idiom reads ~1000× more than it returns and exhausts worker
memory at the full-sweep scale (60M+ trajectories):

```python
ds = xr.concat([xr.open_zarr(z) for z in files], dim="trajectory")
subbasin = release_subbasin(ds.lon.isel(obs=0), ds.lat.isel(obs=0), ...)  # DON'T
```

## Why it OOMs

`ds.isel(obs=0)` on a dask-backed concat is **not** pushed into the read. Each
dask block *is* a full `(10000, 1000)` obs-chunk, so the read task produces a
40 MB array as its **result**, and the `obs=0` slice is a *separate downstream
task*. Those full-chunk results linger in worker memory waiting to be sliced;
under concurrency many co-reside and the allocator keeps the freed ones as
"unmanaged" memory. The worker climbs to the 80 % pause threshold (≈51 GiB of
64), pauses, and the computation deadlocks — whether the `obs=0` value is
`.compute()`-ed standalone or merely fused into a larger `dask.compute`.

## The fix — slice obs inside each per-file read

```python
edge = xr.concat(
    [
        xr.open_zarr(z)[["lon", "lat", "time"]]
        .isel(obs=slice(0, 2))   # 0 for release, 1 for the land test
        .chunk(obs=2)
        for z in files
    ],
    dim="trajectory",
)
lon0 = edge.lon.isel(obs=0, drop=True)
on_land = (edge.lon.diff("obs").isel(obs=0, drop=True) == 0) & (
    edge.lat.diff("obs").isel(obs=0, drop=True) == 0
)
```

Slicing `isel(obs=slice(0, k)).chunk(obs=k)` **before** the concat fuses the
slice into each zarr read: the task result is `(10000, k)` (~tens of KB), not
`(10000, 1000)`. Zarr still decompresses the full chunk internally, but that
buffer is transient inside the single read call and freed immediately — it
never becomes a lingering task result. Measured drop: bytes-crossing-into-dask
fall by ~the obs length (~1000×); peak memory stays flat.

The full-`obs` `ds` is still opened normally and kept lazy — it's needed for
the actual reductions (distance vs. time, density/age histograms) and for the
bounded sampled subset in the raw-trajectory plots. Only the release-time
metadata goes through `edge`.

## Where this applies

| Notebook | Uses the `edge` pushdown |
|----------|--------------------------|
| `020_RawTrajectories` | subbasin + quarter for stratified sampling |
| `022_DispersalDistance` | subbasin, quarter, on-land, release lon/lat (distance + DE filter) |
| `023_Heatmaps` | subbasin, quarter, on-land (kept lazy for `hist_by` masks) |

`021_TimeStats` and `024_BuildHexAggregates` still use a raw `isel(obs=0)`
edge read, but there the full `obs` extent is read *anyway* for the actual
reduction (time statistics, hex aggregation), so the `obs=0` slice fuses into a
streaming reduction rather than standing alone — lower risk. Convert them to
the `edge` pushdown if they ever start pausing at ~51 GiB.

## Cross-references

- [seeding.md](seeding.md) — trajectory layout and the `(trajectory, obs)`
  chunking these reads slice against.
- [distance_calculation.md](distance_calculation.md) — `022`'s great-circle
  distance now reads release position from the attached `release_lon`/
  `release_lat` coords.
