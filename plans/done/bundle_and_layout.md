> **Done.** Outcomes folded into the per-topic docs
> ([seeding](../../docs/seeding.md),
> [hexbinning_and_connectivity](../../docs/hexbinning_and_connectivity.md),
> [2d_field_extraction](../../docs/2d_field_extraction.md));
> the standalone `docs/bundle_and_layout.md` was dropped per
> `plans/wrapup.md §6`. This file is retained as the historical plan record.

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
  obtain/                    # upstream-rebuild recipe, used by twin CI
    download_helcom_subbasins.sh
    download_fucus_shapefile.sh
    download_stokes_sample.sh
    copy_bsh_minimal.sh      # RENAMED from copy_minimal_bsh_data.sh
  fetch_data.sh              # invokes the obtain/ chain end-to-end
  001_prepare_sigma_files.py
  002_download_stokes.py
  003_prepare_2d_fields.py
  004_extract_coastline.py
  010_FucusDispersal_bottom_job.sh
  010_FucusDispersal_surface_job.sh
  020..024_job.sh, submit_all_viz_jobs.sh  # updated path args
data/                        # git submodule → data twin repo
docs/                        # NEW (per AGENTS.md): standalone state docs
ATTRIBUTION.md               # per-dataset blocks (see licensing plan)
LICENSE                      # code, unchanged
notebooks/                   # unchanged logic, renumbered per §5
plans/, pixi.toml, README.md, Notes.md
```

`min_data/` is **deleted**. Everything in it moves to the data twin
(inputs) or the NESH output tree (derived demo outputs).

`data/` keeps its generic name because it's a single submodule mount
point, not a bucket that mixes sources: provenance lives one level
down in the source-prefixed subdir names (§1 twin layout), so
`data/helcom_subbasins/…` already reads as "what data" without the
mount point needing to carry it.

`scripts/` (non-interactive batch: preprocessing + orchestration + obtain)
and `notebooks/` (scientific analysis with narrative) stay split as
they are. The 001–004 preprocessing scripts run unattended and don't
carry analytical narrative, so they belong with the job scripts, not
in `notebooks/`. Revisit if the count of non-notebook scripts grows
enough to warrant subdividing `scripts/`.

Primary setup path is the submodule — no convenience wrapper needed.
`scripts/fetch_data.sh` is the **upstream rebuild** entry point (invokes
`scripts/obtain/*`), not a submodule wrapper; it's what you run when the
twin is unreachable or you want to regenerate the bundled blobs from
source. README.md and AGENTS.md describe the submodule + `git lfs pull`
flow directly — that's self-explaining and doesn't need wrapping.

### Public data twin repo (new, git-LFS)

Hosted at **<https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>**
(institutional GEOMAR GitLab, publicly clonable). Wired into the main
repo as a **git submodule** at `./data/`. Because the twin is public,
the recursive-clone hazard that would normally rule out submodules
doesn't apply — no auth wall to fail against.

Directory names carry their source as a prefix (`helcom_…`, `fucus_…`,
`bsh_…`, `cmems_…`) so provenance is visible without cross-referencing
the attribution file.

```
helcom_subbasins/                # HELCOM subbasins 2022 level 2
fucus_redlist_shapefile/         # REDLIST_SIS_Macrophytes.{shp,dbf,shx,prj,...}
bsh_coastline/                   # derived via scripts/004 but bundled for speed
  coastline.geojson
  coastline_always_wet.geojson
bsh_minimal/                     # demo subset for 003 / 010 / 025 smoke runs;
                                 # mirrors the on-NESH BSH layout so notebooks
                                 # don't special-case the twin.
  static_file_fine/              # includes H0 (no separate bsh_h0/ dir)
  static_file_coarse/
  c_file_fine_2020/
  c_file_coarse_2020/
  t_file_fine_2020/
  t_file_coarse_2020/
cmems_stokes_sample/
  baltic_stokes_20200101.nc
derived/                         # 000-style one-shot preprocess outputs
  fucus_release_points.geojson   # built once by scripts/obtain, committed
  helcom_subbasins_simplified.geojson  # same idea, if 000 factors it out
ATTRIBUTION.md                   # mirrors the main repo's, verbatim
README.md
.gitignore                       # os cruft, editor files, build leftovers
.gitattributes                   # see below
```

`.gitattributes` LFS patterns must cover every binary blob the twin
carries: `*.nc`, `*.zarr/**`, `*.shp`, `*.dbf`, `*.shx`, `*.prj`,
`*.sbn`, `*.sbx`, `*.cpg`, `*.qix`, `*.tif`, `*.tiff`. The `.geojson`
files stay as plain text (diffable, small).

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
```

Repo-name underscore mismatch (`fucus-dispersal` vs.
`fucus_dispersal_outputs`) is fine — they're distinct entities; the
outputs dir is the one that'll most often be read as a variable name.

The large on-NESH BSH store (full multi-year `c_file_*`, `t_file_*`)
lives at its existing NESH path and is **not** under `./data/`. `./data/`
always contains the demo subset bundled by the twin; the big store is
addressed by a separate `bsh_root` parameter on the runs that need it.

## 2. Path convention across notebooks/scripts

Three explicit roots, passed as **parameters**, not environment
variables:

- `data_root` — inputs bundled in the twin. Default
  `Path("./data")` (relative to repo).
- `output_root` — heavy outputs. Default `Path("./output")` locally;
  on NESH set to `<work>/2025_fucus_dispersal_outputs/`.
- `bsh_root` — the large on-NESH BSH store. No default; only the
  stages that need it (003, 010) take it.

Every notebook parameters cell takes the roots it uses as plain
papermill parameters (primitives only, per AGENTS.md). Every job script
passes them explicitly on the papermill command line. No
`FUCUS_OUTPUT_ROOT` / `FUCUS_DATA_ROOT` env vars: the papermill
parameter is the single contract, visible in each notebook's
parameters cell.

`scripts/001..004` are CLI tools and take the same roots as `--data-root`
/ `--output-root` args (defaults matching the notebook defaults).

## 3. Fetching inputs

The submodule is the canonical path — self-explaining, no wrapper.

```bash
git clone --recurse-submodules https://github.com/.../2025_fucus-dispersal.git
cd 2025_fucus-dispersal
git -C data lfs pull
```

After a plain clone:
```bash
git submodule update --init data
git -C data lfs pull
```

README.md and AGENTS.md document this directly; both mention that
git-lfs must be installed locally first (e.g. `pixi global install
git-lfs`), then `git lfs install` once per user.

`scripts/fetch_data.sh` is only the **upstream rebuild** path: it runs
the `scripts/obtain/*` chain to reconstruct `./data/` from public
sources (HELCOM, MADS/SYKE, Copernicus), and is also what the twin's CI
runs on schedule to keep the twin's blobs in sync with the recipe. Use
it when the twin is unreachable or when the recipe changes and you need
to refresh locally.

Everything lives under `./data/`, always. There is no cluster-only data
*inside* `./data/` — the only cluster-only thing is the heavy on-NESH
BSH store, which is a separate path (§1 NESH, `bsh_root` in §2) and
never touched by this flow.

Version pinning: the main repo's submodule SHA records the exact twin
commit. Publication freezes both via git tags + Zenodo DOIs. No
`data/VERSION` file.

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

2. **Move 000 to a one-shot data-prep step.** 000 currently builds
   `REDLIST_SIS_Macrophytes.geojson` (plus SH-coast splitting) at run
   time. Instead, run 000 once, commit the resulting geojson to the
   twin under `data/derived/`, and have 010 read that pre-baked file
   directly. The `cell_ID` particle attribute is never read downstream
   — drop it. Same approach for any HELCOM subbasin simplification:
   bake once into the twin, not recomputed every run. The SH-coast
   splitting logic that doesn't fit the bake-once shape stays in
   `notebooks/explore/`.

3. **Keep the two 010 job scripts separate.** The surface and bottom
   runs will diverge (bottom Ekman, boundary layer, Fucus on rough
   substrate — bottom needs its own treatment, not a shared matrix
   sweep). No merge; the earlier proposal to collapse them is
   withdrawn.

4. **Unify `base_path` → `data_root` + `output_root` (+ `bsh_root`
   where relevant)** (§2) across 010/020/021/022/023/024/025.
   Mechanical; do after §4.1 renumbering.

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

Not doing now: experiment tracking (`plans/experiment_tracking.md`).
The hextraj OOM workaround is already fixed upstream — move
`plans/hextraj_hex_counts_oom.md` to `plans/done/` as part of this
migration.

## 5. Attribution

Create `ATTRIBUTION.md` at repo root per the per-dataset blocks in
`plans/data_licensing_public_bundle.md` §"Attribution block drafts".
Same file copied verbatim into the data twin so the licence travels
with the blobs. `LICENSE` (code, MIT) unchanged; add one line to
`README.md` distinguishing code licence (MIT) from data terms
(per-dataset in ATTRIBUTION.md).

**Why `ATTRIBUTION.md` is needed in the main repo even though the twin
ships the data.** The per-dataset licences (CC BY 4.0 for HELCOM and
MADS/SYKE REDLIST; the CMEMS and BSH terms) require attribution on
*both* redistribution and use-in-derived-products. The twin's copy
discharges redistribution. The main repo's copy discharges the
derived-products trigger: we commit freshly-executed `.ipynb` files
alongside their `.md` sources (per AGENTS.md), and those executed
notebooks contain figures rendered from every bundled dataset. A
repo-root `ATTRIBUTION.md` is the single place that satisfies the
licences for that derivative content. No per-plot annotations.

Scope covered: everything the pipeline *touches*, which is a superset
of what's redistributed in the twin — notably the full on-NESH BSH
store is used (003, 010) but not redistributed, and still needs
attribution here.

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
   add `.gitattributes` (LFS patterns — see §1 twin layout),
   `.gitignore`, `README.md`, `ATTRIBUTION.md`. Push
   `helcom_subbasins/` + `fucus_redlist_shapefile/` +
   `bsh_coastline/` geojsons first (small, derived).
2. Rename `scripts/obtain/*` to the 2-3-word forms
   (`download_helcom_subbasins.sh`, `download_fucus_shapefile.sh`,
   `download_stokes_sample.sh`, `copy_bsh_minimal.sh`). Rewrite
   `scripts/fetch_data.sh` to chain them (upstream rebuild only — no
   submodule / LFS wrapping). The twin's CI runs the same chain.
3. Run the 000-style one-shot preprocess locally, commit the resulting
   geojsons to the twin under `data/derived/`.
4. Copy the `bsh_minimal/` subset (including H0 under
   `static_file_{fine,coarse}/`) into the twin.
5. On the main repo: `git submodule add
   https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data data`.
   Delete `min_data/`. Update 010 and 025 to read pre-baked geojsons
   and H0 from the new layout. Unify every notebook's `base_path`
   usage to the three-root convention per §2.
6. On NESH: create `<work>/2025_fucus_dispersal_outputs/`, wipe the
   in-repo `output/` (nothing to preserve), update every job script to
   pass `--output-root` / `--data-root` / `--bsh-root` explicitly on
   the papermill command line.
7. Do the code straightening (§4) — can be spread across follow-up PRs
   once the layout migration is in.
8. Once green on main, move this plan to `plans/done/` (with a
   one-liner pointing to the `docs/` doc that describes the resulting
   state — create `docs/` if it doesn't exist yet, per `AGENTS.md`
   conventions).

