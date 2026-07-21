# Binary trap term from Copernicus Coastal Zones 2018

Re-establish `trap(shore_type)` in the beaching rate as a **binary substrate
criterion** — sandy shores accept beaching, rocky shores do not — driven by
Copernicus Coastal Zones (CZ) 2018 land cover/land use.

## Why this unblocks the term

[docs/beaching.md](../docs/beaching.md) switched `trap` off for a stated
reason, not a structural one:

> the only typing available is BSH's `H0 ≤ 0` tidal-flat flag, which is not a
> retentiveness proxy for Baltic shores — the basin is effectively tide-free,
> so the flag fires essentially only in the German Bight.

CZ 2018 removes exactly that objection. Its level-4 nomenclature splits along
the axis we want, over the whole (EEA) Baltic coast:

| CZ class | code | binary |
|---|---|---|
| Sandy beaches | 6.2.1.1 | **sandy** |
| Shingle beaches | 6.2.1.2 | sandy (see open question) |
| Dunes | 6.2.2 | sandy |
| Sparse vegetation on sands | 6.1.1 | sandy |
| Intertidal flats | 7.2.3 | sandy |
| Salt marshes | 7.2.1 | sandy |
| Bare rocks and outcrops | 6.3.1.1 | **rocky** |
| Coastal cliffs | 6.3.1.2 | rocky |
| Sparse vegetation on rocks | 6.1.2 | rocky |
| Port areas and associated land | 1.2.3.x | rocky (artificial hard shore) |

The remaining ~60 classes (urban, cropland, forest, grassland, …) describe
hinterland, not shore face, and only matter where they happen to abut the
waterline — see "Attributing the shore face" below.

## Architecture: the hook already exists

`024d_BuildBeaching` keeps the term fully wired and degenerate:

```python
trap_flat = 1.0
trap_wall = 1.0
...
flat_at = rast["nearest_flat"][row, col].reshape(ntraj, nobs)
trap = np.where(flat_at, trap_flat, trap_wall).astype("float32")
```

`nearest_flat` is a boolean raster on the BSH grid, indexed by the nearest-land
index the EDT already computes. **The whole change is to swap what fills that
raster** — from `H0 ≤ 0` to a CZ-derived sandy mask — plus renaming the two
weights. No change to the rate model, the EDT, the Stokes sampler, or the
store schema beyond the label vocabulary.

Fortunate alignment: CZ ships in **EPSG:3035**, which is already the CRS the
`024d` raster is built in. No reprojection of the rasterisation target.

## Obtaining the data

No CLMS client library exists (`copernicusmarine` is Marine-only; the
`eea/clms-*-api-client-python` repos are per-product, for HR-SI/HR-WSI). We
drive the REST API directly, which needs `pyjwt`, `cryptography` and
`requests` added to `pixi.toml`.

Auth is an RFC 7523 JWT-bearer grant: sign a ≤1 h assertion with the service
key's RSA private key, POST it to `token_uri`, get a 1 h bearer token.

The flow below was **verified live** against a Baltic-bbox request before
this plan was written, so the endpoint names, axis order and response fields
are observed rather than read off the docs (which are wrong on axis order).

| step | endpoint | note |
|---|---|---|
| token | `POST https://land.copernicus.eu/@@oauth2-token` | `grant_type=urn:ietf:params:oauth:grant-type:jwt-bearer`, RS256, 1 h |
| request | `POST /api/@datarequest_post` | dataset UID + download-info ID + `BoundingBox` |
| poll | `GET /api/@datarequest_status_get?TaskID=…` | one task; `@datarequest_search` returns all |
| fetch | `GET` the `DownloadURL` | public FME results endpoint, no auth |

Terminal states are `Finished_ok` (carries `DownloadURL`, `FileSize`) or a
failure state in `Message`. The cut runs on an FME queue and takes O(10 min);
CLMS also emails the link, so a timed-out poll is recoverable rather than lost.

Coastal Zones 2018 (vector): UID `b818093444c343f0b67ed7e114a18951`,
download-info ID `ba7a3dba-916e-42fb-8fe0-11ed68a02bc6`.

`BoundingBox` axis order is **`[W, N, E, S]`** — the docs label it `[N,E,S,W]`
but their own France example decodes as west, north, east, south, and the
echoed request for our Baltic box confirms it.

`@get-download-file-urls` (a direct link to the 2 GB prepackaged Europe zip)
returns 404 on this deployment; the bbox cut is the working route.

**Service key location.** `.clms.json` at repo root, gitignored. It holds an
RSA private key — never commit, never echo.

## Attributing the shore face

The BSH coastline geojsons are 8 wet-area polygons, not segmented shore, and
`024d` deliberately does not use them (they drop ~36 % of release points).
The EDT already yields, per water cell, the nearest **land cell index**. So:

1. Rasterise the CZ sandy/rocky mask onto the BSH 3035 grid.
2. A land cell's class is the CZ polygon covering its centre.
3. Land cells whose covering polygon is a hinterland class (forest, urban,
   cropland) carry no shore information — resolve by nearest **coastal** CZ
   polygon within a cutoff, else fall through to the default below.

Step 3 is the real design question: CZ has a 0.5 ha MMU and a **10 m minimum
mapping width**, while BSH grid cells are ~1 km. A beach narrower than 10 m is
not mapped at all, and a mapped beach is routinely far below one grid cell. So
the rasterisation is a *presence* question — "is there any sandy shore polygon
along this cell's waterline?" — not a majority-area question. Plan to
rasterise with `all_touched`-style presence semantics against the coastal
class set, not centroid sampling.

## Coverage gap: Russia

CZ covers **EEA38 + UK — Russia is excluded**. Our domain runs to 30.7 °E, so
two stretches have no CZ polygons at all:

- **Kaliningrad oblast** (~19.9–22.9 °E) — and it splits both the Vistula and
  Curonian spits, which are major sandy features whose Polish/Lithuanian
  halves *are* covered. A silent default here would put a hard discontinuity
  mid-spit.
- **Eastern Gulf of Finland** (~28–30.7 °E, Leningrad oblast).

**Decision (interim): fill uncovered coast with `sandy`.** Chosen to unblock
the pipeline, to be replaced by the HELCOM BRISK layers below.

Defensibility differs sharply between the two gaps, and the interim fill
should not be read as equally supported across both:

- *Kaliningrad* — well supported. The Vistula and Curonian spits are
  aeolian sand barriers; "sandy" is very likely the correct class, and it
  also removes the mid-spit discontinuity against the covered Polish and
  Lithuanian halves.
- *Eastern Gulf of Finland* — weakly supported, plausibly wrong. This coast
  is glacially scoured, and its Finnish/Estonian neighbours grade toward
  rock and skerry. A blanket sandy fill here biases beaching **upward** in
  precisely the domain corner most distant from the release areas.

`024d` must therefore emit the covered/uncovered split as a diagnostic — and
`029`/`031` must be able to report beached weight landing on filled coast
separately — in the same spirit as the existing extrapolation diagnostics.
A silent default is not acceptable.

**Target replacement: HELCOM BRISK, Russia only.** CZ grounds the mask; BRISK
fills the Russian coast and nothing else. The layering is *not* the other way
round — see "Why BRISK is not the base layer" below.

BRISK splits each habitat across geometry types, so no single layer is
spatially complete. All four shore-type records are needed:

| record | title | lat range |
|---|---|---|
| `1a02f79a-2b45-4bac-96fc-27b72e71f447` | Sandy beaches polyline | 52.7 – 61.6 °N |
| `c00edbb4-f940-4b26-abc7-c13f796384bd` | Sandy beaches point | 59.2 – 65.9 °N |
| `fee38332-7d77-4818-9703-db921d84055e` | Rocky shores stone reefs polyline | 54.2 – 66.4 °N |
| `637a3bd2-e302-448f-9562-b7056712a69c` | Rocky shores stone reefs polygon | 53.8 – 66.4 °N |

The sandy polyline stopping at 61.6 °N reads as "Finland is missing"; Finnish
sandy shore is in the **point** layer instead. Both Russian gaps need the
union: Kaliningrad (~54.9–55.3 °N) falls in polyline range, the eastern Gulf
of Finland (~60 °N) straddles the polyline/point boundary.

Fetch via the MADS GP service, the same flow as
[`download_helcom_subbasins.sh`](../scripts/obtain/download_helcom_subbasins.sh).

### Why BRISK is not the base layer

Tempting, since it is Baltic-wide with no EEA38 hole. Rejected because its
completeness is only achieved by unioning three incompatible geometry
semantics. A **point** marks a beach *location* with no extent, so it cannot
answer "is there sandy shore along this cell's waterline?" the way a polygon
can — it can only assert presence somewhere nearby. Mixing point, polyline
and polygon into one shore mask means the mask's meaning varies with
latitude, which is precisely the "resolved morphology that isn't there"
failure [docs/beaching.md](../docs/beaching.md) criticises in the `H0` flag.

CZ, by contrast, is uniform polygon coverage with one nomenclature
throughout. It is the defensible base; BRISK is acceptable as a gap fill
because a gap fill only has to beat the alternative of a blanket constant.

Also worth carrying as a cross-check: **Luijendijk et al. 2018** global
sandy-beach classification (satellite-derived, 500 m alongshore) — an
independent third opinion where CZ and BRISK both cover.

## Open questions

- **Shingle beaches (6.2.1.2)** — sandy or rocky? Coarse clastic shores are
  reflective and mobile; for *Fucus* propagule retention they may behave more
  like rock than sand. Currently binned sandy; worth a sensitivity run.
- **Weight values.** The binary criterion as stated ("rocky shores don't
  allow beaching") implies `trap_rocky = 0`, i.e. an absorbing/reflecting
  dichotomy, not a rate ratio. That is a much stronger claim than the
  `trap_flat = 2.0` the old code contemplated, and it makes total beached
  weight sensitive to the sandy fraction. Recommend implementing as a ratio
  with `trap_rocky` a free parameter, swept in `031`, with 0 as one member —
  rather than hard-zeroing.
- **Does CZ resolve enough Baltic shore?** Needs a coverage count before
  committing: what fraction of BSH coastal land cells get a *coastal* CZ
  class at all, vs falling through to hinterland. If that fraction is low the
  whole approach is a presentation of resolved morphology that isn't there —
  the exact failure the doc calls out for the `H0` flag.

## Sequencing

1. `scripts/obtain/download_coastal_zones.sh` — auth + bbox request + poll +
   unpack, idempotent, matching the existing obtain-script conventions.
2. Coverage count (open question 3) — **gate**. If CZ doesn't resolve the
   Baltic shore, stop and report rather than proceed.
3. Derive a small Baltic-only sandy/rocky shore geojson for the data twin,
   same pattern as the derived release-points geojson.
4. Swap `nearest_flat` → `nearest_sandy` in `024d`; rename weights; add the
   coverage diagnostic.
5. Sweep `trap_rocky` in `031`; update `docs/beaching.md` provenance table —
   `trap` moves from "deliberately switched off" to active, and the defence
   in the doc must be rewritten, not just deleted.
6. Extend `ATTRIBUTION.md` (CZ is Copernicus/EEA) in the same pass.
