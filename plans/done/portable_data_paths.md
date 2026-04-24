# Portable data paths for viz notebooks

## Problem

Every viz notebook hard-codes a personal scratch path:

```python
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"
```

…and then reads input shapefiles from `base_path / "data" / ...`:

| File | Notebooks |
|------|-----------|
| `data/Fucus_location_shp/REDLIST_SIS_Macrophytes.shp` | 020, 022, 023 |
| `data/HELCOM_subbasins_2022_level2/HELCOM_subbasins_2022_level2.shp` | 020, 022, 023, 024 |

This breaks reproducibility for any other user on the cluster (or off it):
the path is tied to `smomw122`'s `gxfs_work`, which is not readable by
collaborators and not guaranteed to outlive the user account.

The shapefiles themselves are likely too large / not-our-IP to track in git
directly, so we need *some* mechanism to share them across user accounts on
the cluster.

## Options to evaluate

### A. Shared read-only mount

- Park a single canonical copy under a project-shared NESH path (e.g.
  `/gxfs_work/geomar/<project-share>/data/...`) with group-read perms.
- Notebook parameter becomes `data_root = os.environ.get("FUCUS_DATA_ROOT", "<shared default>")`.
- Pros: no per-user setup, no duplication, git stays clean.
- Cons: requires a project share to exist; permissions/quota coordination.

### B. Tracked manifest + materialise step

- Track a small YAML/TOML manifest in git: shapefile sources, URLs/DOIs,
  checksums, target subdirs.
- Add a pixi task (e.g. `pixi run fetch-data`) that reads the manifest and
  populates `$FUCUS_DATA_ROOT/data/...` on first use.
- Pros: fully reproducible from a fresh clone; works off-cluster too.
- Cons: needs upstream URLs that won't rot; some shapefiles may not have
  stable public sources (HELCOM yes, Fucus REDLIST shp ?).

### C. Hybrid

- Manifest in git for fetchable assets (HELCOM).
- Shared mount for the rest (Fucus shapefile if unfetchable).
- Single `data_root` env var resolves both.

## Decision

**Deferred.** Park here until we agree which option fits the cluster's
sharing story and which assets actually have stable upstream URLs.

## When picked up

Touch every viz notebook's parameters cell + the two shapefile reads in
020/022/023 and the HELCOM read in 024. Mechanical change once the
mechanism is chosen.
