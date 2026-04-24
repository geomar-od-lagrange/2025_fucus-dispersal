# Agent guidelines for this project

(This file is `AGENTS.md` on disk and symlinked to `CLAUDE.md`.)

## Scope

Parcels-based Lagrangian dispersal study for *Fucus vesiculosus* in the
Baltic Sea, driven by BSH operational-model currents (HBMnoku) plus
optional Stokes drift from the CMEMS Baltic Wave Hindcast. The repo
covers the full pipeline: input preparation, Parcels runs, per-trajectory
visualisations, and a hex-aggregated dispersal store for interactive
source↔sink queries. Intended to support a published study.

## Principles

**This is pre-alpha research code.** No installed base, no users to migrate.
Internal API changes are free. Anything prefixed with `_` is internal.

**Greenfield mindset.** If the current shape is in the way of the right shape,
reshape it. Don't add workaround constraints when restructuring eliminates the
problem. Deletions, renames, and rewrites are the normal mode. **No produced
data on NESH needs preserving** — filenames, output layout, and notebook
numbering are free to change.

**"Already X" is not a justification.** When deciding a name, layout, or API
shape, don't defend an option on the grounds that "AGENTS.md already uses X",
"the notebooks already call it X", or "historical usage is X". That's inertia,
not an argument. State the merit-based case at the moment of the decision, or
concede there isn't one and change course. This applies equally to plan docs,
PR descriptions, and code review — scrub "already…" reasoning wherever it
appears.

**Be ruthless about dropping dead code.** Patch sparingly; rewrite when the
abstraction is wrong. Prefer clean parameter plumbing over clever hacks — no
monkey-patches, global state swaps, or closure tricks when passing a parameter
through the call chain is cleaner.

**Be diligent about follow-through.** When you touch a name, signature, or
path, grep for every reference and update them in the same pass. Don't leave
stale imports, dead references, or half-updated docs for a later cleanup step.

## Agent workflow

**Planning before code:** Write plans to `plans/*.md` before touching source.
Don't skip planning for complex changes; don't let implementation agents make
architectural decisions unguided.

**Model choice:** Use a lighter model for mechanical and verification tasks.
Reserve more capable models for architecture, design decisions, and judgment
calls.

**Always review after implementation:** A separate review agent should examine
the result. This catches both conceptual mistakes and quality issues.

**Experimental validation:** Use `tmp_*/` directories to prove ideas before
committing to architecture changes. Once validated, clean up or move to
permanent locations.

## Tooling

Use **pixi** for environment management. Run all commands with `pixi run
<command>` directly — this repo does not use `[tasks]` in `pixi.toml`.
Prefer conda packages; use pypi only when no suitable conda package exists.

**Shell.** Always work from the repo root. Don't `cd` into subdirectories
between commands — pass relative paths to tools (e.g.
`pixi run jupytext --sync --execute notebooks/foo.md`) or use a subshell
`(cd notebooks && ...)` for one-off needs. Don't pass `-C <path>` to
`git` for ordinary repo-root operations; it already operates on the cwd.

**`./data/` is a submodule — pwd matters.** `./data/` points at the data
twin repo (see §Data access), so git commands inside `data/` act on the
*submodule*, not the main repo. A plain `git status` in `data/` reports
only twin-side changes and ignores every modification in the main tree;
`git rm --cached -r data/` from inside `data/` resolves `data/data/` and
fails outright. Never leave cwd sitting inside `data/` between unrelated
commands — if you need to operate on the submodule (fetch, pull, the
rare twin push), use an explicit `(cd data && git <cmd>)` subshell so
cwd returns to the repo root automatically, or prefix a single command
with `git -C data <cmd>` (this is the one exception to the "no `-C`"
rule above).

## Pipeline stages

A single linear pipeline, named with 3-digit stage numbers and gaps of 10
for inserts. Region/year belong to the parameters cell, never in the
filename.

- `scripts/00x_*.py` / `notebooks/00x_*.md` — inputs and preprocessing
  (sigma cleanup, Stokes download, 2D-field extraction, BSH coastline
  extraction).
- `scripts/010_FucusDispersal_job.sh` + `notebooks/010_FucusDispersal.md`
  — Parcels runs, papermill-swept across (year, doy, regime,
  velocity_factor). Writes trajectory zarrs under `output_root/Trajectories/`.
- `notebooks/020..023_*.md` — per-trajectory visualisations (raw lines,
  time stats, distance vs. time, density + mean-age heatmaps).
- `notebooks/024_*.md` — build the hex-aggregated dispersal store
  (`key.parquet` + `counts/` partitions).
- `notebooks/025_*.md` — hex heatmaps driven by the aggregate store.

Naming is currently mid-migration (see `plans/bundle_and_layout.md` for
the target). Ordering may drift during the prod-prep refactor; when in
doubt, defer to the plan.

## Data access

Inputs live in the **data twin repo** at
<https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>, wired in
as a git submodule at `./data/`. Clone with
`git clone --recurse-submodules` or run `scripts/fetch_data.sh` after a
plain clone. The submodule carries HELCOM polygons, the Fucus shapefile,
BSH coastline geojsons, a BSH H0 slice, a minimal BSH demo subset, and
a Stokes sample — all under the attributions listed in
`ATTRIBUTION.md`.

`scripts/obtain/*.sh` is the canonical recipe for each input from
upstream public sources (HELCOM, MADS/SYKE, Copernicus, BSH). These are
the fallback when the submodule is unreachable, and the recipe the
twin's CI runs to rebuild its blobs. Keep them in sync with how the
twin was produced.

Cluster outputs live **outside** the repo tree on NESH, at
`<work>/2025_fucus_dispersal_outputs/`. Every job script sets `output_root`
as a shell variable and passes it to papermill via `-p output_root
${output_root}`; there is no `$FUCUS_OUTPUT_ROOT` environment variable.
The papermill parameter is the single contract — visible in each notebook's
parameters cell. Never commit output data, and don't leak `output/` into
the repo directory.

**BSH H0 semantics.** `H0` is *minus z of the sea floor* on a z-up axis
with `z = 0` at MSL — **not** water depth. Always-wet cells have `H0 > 0`
and `H0` coincides with depth; tidal-flat cells can have `H0 < 0`.
Anywhere you want depth, filter to `H0 > 0` first (or use the
always-wet mask). Mixing tidal flats into a mean will pull results
toward zero or flip their sign.

- Subset with `.sel()` / `.isel()`; keep loads lazy until needed.
- For CMEMS data, prefer `copernicusmarine.open_dataset` (or the
  existing `scripts/002_download_stokes.py` for bulk pulls) with
  minimal arguments.

## Conventions

### Code

**Be careful with generated or derived artifacts.** Some files (executed
notebooks under `notebooks_executed/`, derived geojsons, zarr stores)
are cached outputs, not hand-written source. Before editing them, check
how they were produced and whether changing the source is the right fix
instead.

**Shared helpers.** `notebooks/helpers.py` holds genuinely shared
utilities (`load_trajectories`, `mask_land_seeded`,
`attach_release_metadata`, `QUARTER_LABELS`). Don't factor *duplicated
setup* — parameter cells, path construction, Dask client boilerplate —
into `helpers.py`; notebooks should read standalone. Extract to
`helpers.py` only when the logic is non-trivial and reused verbatim.

### Notebooks

See the **jupytext skill** (`.agents/skills/jupytext/SKILL.md`) for the full
workflow: creating, syncing, executing, and fixing notebooks.

**Naming.** 3-digit stage number with 10-step gaps for inserts
(`000_…`, `010_…`, `020_…`). Region/date/regime belong to the
parameters cell, not the filename.

**Execution from the repo root.** Always invoke papermill with
`--cwd notebooks/` so `from helpers import …` resolves and relative
paths in the notebook (`data/foo.shp`) resolve under the repo.

- The `.md` is the source of truth — always edit it, never the `.ipynb`.
  Commit **both** the `.md` and the freshly-executed `.ipynb` so rendered
  figures show on GitHub. Execute with
  `pixi run jupytext --sync --execute <nb>.md`. Papermill only for parameter
  injection or when streaming progress matters. Do not use
  `jupyter nbconvert --execute`.
- When using **papermill** from the repo root, always pass `--cwd notebooks/`
  so relative paths in the notebook resolve against the notebook's
  directory rather than wherever the shell happens to be.
- Markdown cells for narrative; clean code cells for execution.
- Well-scoped cells — don't mix imports, parameters, and calculations.
- Every notebook must have one early parameters cell tagged `"parameters"`
  containing only primitive assignments (`int`, `float`, `str`, `bool`,
  `None`). All calculations, transformations, and derived values belong in
  subsequent cells. This keeps notebooks papermill-compatible for parameter
  sweeps.
- **Never write summary cells with prose that assumes results.** Summary cells
  must compute and print dynamically.
- Use xarray, pandas etc. _public_ API. Example: `ds.lon.isel(trajectory=0)`
  instead of `ds.lon.values[0, :]`.
- Prefer lazy `.where(mask)` on `(trajectory, obs)` arrays over eager
  subsetting when building per-scope views — fuses block-wise with the
  trajectory graph; eager fancy indexing doesn't.
- Dask client bootstrap: `Client(scheduler_file=os.environ["SCHEDULER_FILE"])`
  when the SLURM job scripts set it, else `Client(ip="0.0.0.0")`. Use a
  single RNG seed per notebook run (`np.random.default_rng(seed)`) and
  print it for reproducibility.
- After fixing bugs, rerun immediately without asking.

### Plotting

Default to vanilla matplotlib / xarray / cartopy. **Don't pass styling
kwargs that override defaults until a real need shows up.** Concretely:

- No `color=`, `facecolor=`, `edgecolor=` overrides. When categories
  must be distinguished, take colours from the default cycle
  (`plt.rcParams["axes.prop_cycle"].by_key()["color"]`) by index — never
  hardcode names like `"tab:green"` or `"lightgray"`.
- No `cmap=` overrides — use the default colormap.
- No `figsize=` or `dpi=` overrides — use the matplotlib defaults.
- No `cbar_kwargs={"shrink": ...}` / `{"aspect": ...}` or similar size
  knobs — let xarray/matplotlib place the colorbar.
- No `linewidth=`, `markersize=`, `linestyle=` overrides.
- Axis / colorbar labels: leave whatever xarray derives from
  `long_name` / `units` attributes. Only call `ax.set_title(...)` for
  context the data doesn't carry.
- Cartopy: `ax.add_feature(cfeature.LAND)` and
  `ax.add_feature(cfeature.COASTLINE)` with no styling kwargs. Use
  Natural Earth or OSM tiles when a basemap is needed.
- Prefer xarray's built-in plotting (`.plot`, `.plot.line`, `.plot.scatter`,
  faceting via `col=` / `row=` / `hue=`) over raw matplotlib.

When a default genuinely doesn't read, prefer reshaping the data
(faceting, hue) before reaching for explicit styling. If you do override
a default, leave a one-line comment explaining why.

Some existing notebooks (notably the hex heatmaps) still carry overrides
from before this convention. They'll be pruned during the viz cleanup
tracked in `plans/bundle_and_layout.md` §4.

### Documentation

`docs/*.md` contains standalone documentation for the current state of the
code. Each doc should make sense on its own without referencing previous
implementations, changelogs, or development history. Explain design choices by
comparing alternatives and their trade-offs, not by narrating what changed.
Git history is the changelog; docs describe what *is*, not what *was*.

`plans/*.md` describe intent before implementation. When a plan is
implemented: write a corresponding `docs/` file, move the plan to
`plans/done/`, and add a one-liner at the top pointing to the doc.
Agents get context by reading `docs/*.md` (what is) + open `plans/*.md`
(what's next).

Use markdown relative links when referencing other files in `plans/` and
`docs/`. Example: `[bundle_and_layout.md](../plans/bundle_and_layout.md)`
from a doc, `[seafloor-location-H0-semantics.md](seafloor-location-H0-semantics.md)`
within `plans/`.

## Attribution

All redistributed data carries per-dataset terms listed in
`ATTRIBUTION.md` at the repo root. When adding a new input, check its
licence first and extend `ATTRIBUTION.md` (and the data twin's copy) in
the same PR. Code licence (MIT) applies only to code, not to data;
`README.md` makes the distinction explicit.
