# Public-repo refactor: data bundle + NESH output layout + code cleanups

One plan covering the merge-to-main preparation on `wr/prod-prep`. Aimed
at clearing three entangled concerns in one go: (1) where the *input*
data lives (public twin repo vs. repo vs. NESH shared paths), (2) where
the *output* data lives on NESH, (3) notebook/script straightening that
was deferred while the viz work settled.

Session context:
- Licensing cleared for public redistribution of every input we use —
  see `plans/data_licensing_public_bundle.md`. The data twin can be a
  public repo, not a private one.
- No produced data on NESH needs preserving. Filenames, output layout,
  and notebook numbering are free to change.
- Workflow is monotonic up to publication: no branch maintenance, no
  data-version pinning needed (see §3).

Supersedes `plans/portable_data_paths.md` (already moved to
`plans/done/`).

## 1. Target layout

### Public main repo (this one) — code only

```
AGENTS.md                    # agent contract (symlinked to CLAUDE.md) — DONE
CLAUDE.md -> AGENTS.md       # symlink — DONE
.agents/skills/              # shared skills (jupytext) — DONE
.claude/skills -> ../.agents/skills   # compat symlink — DONE
.gitmodules                  # NEW: data/ submodule → data twin repo
scripts/
  fetch_data.sh              # NEW: one-shot setup wrapper (§3)
  obtain/
    helcom.sh                # NEW: download + convert HELCOM shapes
    fucus.sh                 # NEW: download MADS/SYKE macrophyte shapefile
    stokes_sample.sh         # NEW: thin wrapper around scripts/002 for demo window
    bsh_minimal.sh           # RENAMED from copy_minimal_bsh_data.sh; cluster-only
  001_prepare_sigma_files.py
  002_download_stokes.py
  003_prepare_2d_fields.py
  004_extract_coastline.py
  010_FucusDispersal_job.sh  # MERGED from the two twin scripts
  020..024_job.sh, submit_all_viz_jobs.sh  # updated base_path logic
data/                        # git submodule → data twin repo
ATTRIBUTION.md               # NEW: per-dataset blocks (see licensing plan)
LICENSE                      # code, unchanged
notebooks/                   # unchanged logic, renumbered per §5
plans/, pixi.toml, README.md, Notes.md
```

`min_data/` is **deleted**. Everything in it moves to the data twin
(inputs) or the NESH output tree (derived demo outputs).

### Public data twin repo (new, git-LFS)

Hosted at **<https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>**
(institutional GEOMAR GitLab, publicly clonable). Wired into the main
repo as a **git submodule** at `./data/`. Because the twin is public,
the recursive-clone hazard that would normally rule out submodules
doesn't apply — no auth wall to fail against.

```
HELCOM_subbasins_2022_level2/    # was gitignored here; lives only in twin
Fucus_location_shp/              # REDLIST_SIS_Macrophytes.{shp,dbf,shx,prj,...}
BSH_model_coastline/             # derived via scripts/004 but bundled for speed
  coastline.geojson
  coastline_always_wet.geojson
bsh_h0/                          # just the H0 files 025 needs for the key
  H0_file_fine.nc
  H0_file_coarse.nc
bsh_minimal/                     # demo subset for 003 / 010 / 025 smoke runs
  c_file_{fine,coarse}_2020/
  static_file_{fine,coarse}/
  t_file_{fine,coarse}_2020/
stokes_sample/
  baltic_stokes_20200101.nc
ATTRIBUTION.md                   # mirrors the main repo's, verbatim
README.md
.gitattributes                   # *.nc, *.shp, *.zarr → lfs
```

### NESH layout — outputs leave the repo

Currently on NESH:
```
<work>/2025_fucus-dispersal/         # repo clone + output/ in tree
```

Target:
```
<work>/2025_fucus-dispersal/           # just the repo clone (+ fetched data/)
<work>/2025_fucus_dispersal_outputs/   # all heavy outputs, outside the repo
    2d_fields/                         # from scripts/003
    sigma/                             # from scripts/001
    stokes/                            # from scripts/002 (full-year)
    Trajectories/<regime>/<year>/*.zarr   # from notebook 010
    HexAggregates/r{N}km_v1/           # from the aggregate-store notebook
    debug/                             # bottom-stationary audit etc.
```

Repo-name underscore mismatch (`fucus-dispersal` vs.
`fucus_dispersal_outputs`) is fine — they're distinct entities; the
outputs dir is the one that'll most often be read as a variable name.

## 2. Path convention across notebooks/scripts

Two explicit roots, no more `base_path`:

- `data_root` — inputs. Default `Path("./data")` (relative to repo).
  Used for reading HELCOM shapes, Fucus shapes, BSH coastline geojsons,
  H0 files.
- `output_root` — outputs. Default
  `os.environ.get("FUCUS_OUTPUT_ROOT", "./output")`. On NESH, set to
  `<work>/2025_fucus_dispersal_outputs/`.

Every notebook parameters cell takes both. Every job script exports
`FUCUS_OUTPUT_ROOT` once at the top and passes `data_root` through
papermill. Scripts/001–004 take them as CLI args (already close —
tidy the defaults).

## 3. Fetching inputs: submodule primary, `scripts/fetch_data.sh` wrapper

The submodule is the canonical path. A recursive clone plus LFS pull
populates `./data/` and pins the exact data SHA into the main repo's
history — no separate version file needed. Publication freezes both
via git tags + Zenodo DOIs.

Primary path:
```bash
git clone --recurse-submodules https://github.com/.../2025_fucus-dispersal.git
cd 2025_fucus-dispersal
git -C data lfs pull          # or: pixi run fetch-data
```

`scripts/fetch_data.sh` is a thin convenience wrapper that:
1. Runs `git submodule update --init data` (no-op if already there).
2. Runs `git -C data lfs pull`.
3. On submodule failure (e.g. GitLab temporarily unreachable), falls
   back to the `scripts/obtain/*.sh` chain to fetch from upstream
   public sources (HELCOM, MADS, Copernicus). Rebuild populates the
   same paths under `./data/` and the notebooks don't care which path
   filled them.
4. Skips tier-3 (BSH bulk) unless `WITH_BSH=1`; that's cluster-only.

The rebuild fallback also doubles as the **recipe of record** for the
data twin. Twin-side CI runs the same scripts on schedule (and on
main-repo recipe changes) to keep twin HEAD in sync with main HEAD —
so the submodule pointer on main is always advanceable to a fresh twin
commit when recipes change.

No `data/VERSION` file: the submodule SHA is the version, and "main
HEAD is the only supported vintage" is the workflow.

## 4. Code straightening (sequenced)

All code/structural changes deferred while the viz work completed;
pick them up in this order. No data on NESH needs preserving, so
renames are free.

1. **Reorder hex notebooks.** The current 025 (build aggregate store)
   produces exactly what the current 024 (hex heatmaps) aggregates
   on-the-fly from zarrs. Swap:
   - `024_BuildHexAggregates.md` (was 025): build key + counts parquet.
   - `025_HexHeatmaps.md` (was 024): read from parquet, render. Drops
     the dask-cluster machinery, the NaN-int sentinel workaround, and
     `counts_by_scope`; becomes a geopandas read + groupby-sum +
     `log_density_plot`.

2. **Fold 000 into 010 (or drop).** 010 reads
   `REDLIST_SIS_Macrophytes.geojson` which only 000 produces, and uses
   000's `CELLID = index`. The `cell_ID` particle attribute is never
   read downstream. Cheapest: 010 reads the `.shp` directly, drops
   `cell_ID`, and 000 moves to `notebooks/explore/` (the SH-coast
   splitting logic is still useful there).

3. **Merge the two 010 job scripts.** `010_FucusDispersal_bottom_job.sh`
   and `..._surface_job.sh` differ only in the trailing regime/vf loop.
   One script, regime matrix as argv.

4. **Unify `base_path` → `data_root` + `output_root`** (§2) across
   010/020/021/022/023/024/025. Mechanical; do after §4.1 renumbering.

5. **Pin output filename convention** and extract
   `parse_zarr_stem(path)` into `helpers.py`. The new 024 relies on
   parsing `Fucus_BSH_YYYYMMDD_{regime}_…zarr`; make it the one place
   that knows the format.

6. **Jupytext-ify 010** to match the rest of the tree. Lets the
   `# Parameters` cell carry `tags=["parameters"]` properly.

7. **Regime discovery everywhere.** `021` / `022` already do
   `sorted(p.name for p in trajectory_root.iterdir())`. Apply the same
   to the new 024 (aggregate build) so it doesn't hard-code the three
   names — this also unblocks future regime variants (e.g. bottom at
   vf=0.97).

Not doing now: experiment tracking (`plans/experiment_tracking.md`),
hextraj OOM upstream PR (`plans/hextraj_hex_counts_oom.md` — keep the
downstream workaround until upstream accepts a patch).

## 5. Attribution

Create `ATTRIBUTION.md` at repo root per the per-dataset blocks in
`plans/data_licensing_public_bundle.md` §"Attribution block drafts".
Same file copied verbatim into the data twin so the licence travels
with the blobs. `LICENSE` (code, MIT) unchanged; add one line to
`README.md` distinguishing code licence (MIT) from data terms
(per-dataset in ATTRIBUTION.md).

Two outstanding attribution inputs still needed:
- BSH minimal: enumerate the year/grid/file scope actually shipped in
  the twin once the subset is pinned.
- Fucus/REDLIST: enumerate species covered from the `.dbf` columns.

## 6. Migration sequence

Rough order; each step self-contained and mergeable.

**Already landed** (this branch):
- `AGENTS.md` + `CLAUDE.md` symlink at repo root.
- `.agents/skills/jupytext/SKILL.md` + `.claude/skills` symlink.
- `.gitignore` updated to the "ignore everything under `.{agents,claude}/`
  except `skills/`" pattern.
- `plans/portable_data_paths.md` moved to `plans/done/` (superseded).

**To do:**

1. Initialise the data twin at
   <https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>:
   add `.gitattributes` (LFS patterns), `README.md`, `ATTRIBUTION.md`.
   Push HELCOM + Fucus + BSH-coastline geojsons first (small, tier
   1 / derived).
2. Write `scripts/obtain/helcom.sh`, `scripts/obtain/fucus.sh`, plus
   `scripts/fetch_data.sh` (submodule + LFS + rebuild fallback, §3).
   This is the recipe that the twin's CI will also run.
3. Copy the BSH H0 and BSH minimal subsets into the twin. Add
   `scripts/obtain/bsh_minimal.sh` (renamed from
   `copy_minimal_bsh_data.sh`; cluster-only guarded).
4. On the main repo: `git submodule add https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data data`.
   Delete `min_data/`. Update 025 to read H0 from `data/bsh_h0/`.
   Update every notebook's `base_path` usage per §2.
5. On NESH: create `<work>/2025_fucus_dispersal_outputs/`, wipe the
   in-repo `output/` (nothing to preserve), update every job script
   to export `FUCUS_OUTPUT_ROOT`.
6. Do the code straightening (§4) — can be spread across follow-up PRs
   once the layout migration is in.
7. Once green on main, move this plan to `plans/done/` (with a
   one-liner pointing to the `docs/` doc that describes the resulting
   state — create `docs/` if it doesn't exist yet, per `AGENTS.md`
   conventions).

## 7. Open questions (don't block)

- Whether to drop the BSH minimal from the twin entirely and keep it
  cluster-only (saves LFS cost, loses laptop smoke-testing). Default
  is "include it"; revisit if LFS size becomes a problem.
- `FUCUS_OUTPUT_ROOT` default — env var only, or sibling-of-repo
  convenience (`../2025_fucus_dispersal_outputs/`)? Env var is more
  portable; sibling is zero-config. Proposing env var with the sibling
  path as a suggested value in `README.md`.
