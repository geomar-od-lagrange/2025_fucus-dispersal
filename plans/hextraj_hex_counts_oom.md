# hextraj: `hex_counts` OOMs on dask-backed inputs

## Problem

`hex_counts(hex_ids, reduce_dims=..., hp=...)` calls `hex_ids.values.ravel()`
when reducing over all dimensions, which eagerly materialises the entire dask
array on a single worker.  For large trajectory datasets this is hundreds of
GiB and causes a `MemoryError` from dask distributed:

```
MemoryError: Task 'finalize-hlgfinalizecompute-...' has 757.25 GiB worth of
input dependencies, but worker ... has memory_limit set to 32.00 GiB.
```

The relevant code path in `hex_analysis.py`:

```python
# reduce_dims == all dims  →  full materialisation
hex_array = hex_ids.values.ravel()          # <-- pulls entire dask array
counts = pd.Series(hex_array).value_counts(sort=False)
```

## Workaround (applied in 024_HexHeatmaps)

Flatten to a `dask.dataframe` Series and do distributed `value_counts()`, then
pass only the small aggregated result to `hex_counts` for geometry:

```python
flat = hex_ids.data.ravel()                 # stays lazy (dask array)
vc = dd.from_dask_array(flat, columns="hex_id")["hex_id"].value_counts().compute()
vc = vc[vc.index >= 0]
gdf = hex_counts(pd.Series(vc.index, dtype=np.int64), hp=hp)
gdf["count"] = vc.reindex(gdf.index).values
```

## Suggested fix in hextraj

When `hex_ids` is dask-backed and `reduce_dims` covers all dimensions,
`hex_counts` should do distributed counting instead of materialising.  Minimal
change: detect `dask.is_dask_collection(hex_ids)` and use the
`dask.dataframe.value_counts` path shown above before calling
`_build_counts_geodataframe`.

The partial-reduction branch (groupby over `keep_dims`) has the same issue —
it iterates `hex_ids.groupby(keep_dims)` which also triggers compute.

## Action

- [ ] Open issue on hextraj repo with reproducer
- [ ] Optionally submit PR with the dask-aware path
