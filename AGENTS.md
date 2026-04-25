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

**Dev-tool installs.** For tools that don't touch the project env
(git-lfs is the canonical example), recommend `pixi global install
<tool>` from conda-forge. Don't suggest `pip`, `brew`, or `apt` — not
even as "or" fallbacks — in README, docs, or chat replies. Tools that
*do* need project dependencies (jupytext, papermill, jupyter, python)
stay in the project env; invoke them via `pixi run <tool>` so they see
the kernel and the full dep tree.

**Always invoke through `pixi run`.** Don't resolve a python executable
once (e.g. `pixi run which python`) and then call that binary directly
in subsequent commands. Conda packages like GDAL, proj, and cartopy
ship `etc/conda/activate.d/*.sh` scripts that set `PROJ_DATA` /
`GDAL_DATA` and similar; bypassing `pixi run` skips those hooks and
produces cryptic env-setup failures that look like missing data files.
Same applies to `jupyter`, `jupytext`, `papermill`, `python` itself —
prefix every invocation with `pixi run`.

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
- `scripts/010_FucusDispersal_{surface,bottom}_job.sh` +
  `notebooks/010_FucusDispersal.md` — Parcels runs, papermill-swept
  across (year, doy, regime). Writes trajectory zarrs under
  `output_root/Trajectories/<regime>/<year>/`.
- `notebooks/020..023_*.md` — per-trajectory visualisations (raw lines,
  time stats, distance vs. time, density + mean-age heatmaps).
- `notebooks/024_*.md` — build the hex-aggregated dispersal store
  (`key.parquet` + `counts/` partitions).
- `notebooks/025_*.md` — hex heatmaps driven by the aggregate store.

## Data access

Inputs live in the **data twin repo** at
<https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>, wired in
as a git submodule at `./data/`. Clone with
`git clone --recurse-submodules` or run `scripts/fetch_data.sh` after a
plain clone. Top-level layout (`<source>_<dataset>`):
`helcom_subbasins_2022/` (HELCOM polygons),
`helcom_fucus_redlist/` (Fucus shapefile + the derived release-points
geojson baked by 000), `bsh_hbmnoku_static/` (BSH static grids
including H0 in `static_file_{fine,coarse}/` plus the wet-cell
coastline geojsons), `bsh_hbmnoku_demo/` (one-day BSH HBMnoku c/h/t/z
demo subset), and `cmems_stokes_sample/` (one-day Stokes drift) — all
under the attributions listed in `ATTRIBUTION.md`.

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

**BSH H0 semantics.** See [docs/h0_semantics.md](docs/h0_semantics.md).

- Subset with `.sel()` / `.isel()`; keep loads lazy until needed.
- For CMEMS data, prefer `copernicusmarine.open_dataset` (or the
  existing `notebooks/002_download_stokes.py` for bulk pulls) with
  minimal arguments.

## Conventions

### Code

**Be careful with generated or derived artifacts.** Some files (executed
notebooks under `notebooks_executed/`, derived geojsons, zarr stores)
are cached outputs, not hand-written source. Before editing them, check
how they were produced and whether changing the source is the right fix
instead.

**Trust derivations; let downstream raise.** Don't add a "plausibility
check" that re-reads the source to verify a value precomputed for
performance — that defeats the precompute. The error path catches
mismatches.

**Inline at the consumer.** Don't pre-allocate arrays that are sliced
once each. Move the draw, lookup, or computation into the loop body
where it's used — pre-allocation implies a shared contract that
doesn't exist.

**Notebook-local utilities.** Notebooks own their utilities. Four
tests gate any extraction:

1. **Default: notebooks own their utilities.** A shared module isn't
   forbidden, but the bar is really high. Three copies of a 12-line
   function isn't enough — the candidate must (a) clearly pass the
   idiom and one-concern tests below, and (b) be duplicated in enough
   places that copies actively drift or bug fixes touch N callsites.
   Cross-notebook duplication is the default state; standalone
   reading is worth that price.
2. **Idiom test before `def`-ing.** Would a reader recognize the body
   as an idiom faster than the call? `sorted(path.glob("*.zarr")) +
   xr.concat(...)` reads as itself; the call `load_trajectories(path)`
   is a layer of indirection over an idiom. ~12 lines of STRtree +
   NaN handling earn a name; 2–3 lines of xarray idioms don't. When
   in doubt, inline.
3. **One concern per function.** Don't bundle unrelated assignments
   under one name. If you're reaching for `attach_metadata`,
   `prepare_dataset`, `setup_*`, you're probably bundling — split
   into one function per real concern, or inline both.
4. **Don't mutate coords for presentation.** If labels only matter at
   plot time, set them at plot time (`FacetGrid.set_title`,
   `ax.set_xticklabels`). Mutating a possibly-dask-backed coord via
   `.values` + comprehension is eager and brittle, and the dataset
   doesn't need to know about your tick labels.

### Notebooks

See the **jupytext skill** (`.agents/skills/jupytext/SKILL.md`) for the full
workflow: creating, syncing, executing, and fixing notebooks.

**Naming.** 3-digit stage number with 10-step gaps for inserts
(`000_…`, `010_…`, `020_…`). Region/date/regime belong to the
parameters cell, not the filename.

**Execution from the repo root.** Always invoke papermill with
`--cwd notebooks/` so relative paths in the notebook (`../data/foo.shp`,
`../output/...`) resolve against the notebook's directory.

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
a default, leave a one-line comment explaining why and document the
rationale in [docs/visualisations.md](docs/visualisations.md).

### Documentation

`docs/*.md` contains standalone documentation for the current state of the
code. Each doc should make sense on its own without referencing previous
implementations, changelogs, or development history. Explain design choices by
comparing alternatives and their trade-offs, not by narrating what changed.
Git history is the changelog; docs describe what *is*, not what *was*.

**Keep docs concise.** Lead with what the code does (file path + a
small table or a few-line code block), pair with a short *why* recap.
Default 50–150 lines per doc; hit 200+ only if the topic is genuinely
large. Cut alternatives narratives to "Considered X, rejected because
Y" one-liners. Front-matter docs that bloat past this drift fast and
bury the implementation pointer that matters.

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

**Keep `ATTRIBUTION.md` minimal.** One short block per dataset: source,
catalog/DOI URL, and the literal attribution string the licence
requires. ≤5 lines per dataset. Skip introductions, licence-name/URL/
text, derivative-works reasoning, species/column lists, contact emails
(unless contact *is* the authoritative attribution channel — BSH's
`opmod@bsh.de` is the one current exception). Detailed licence
reasoning belongs in `plans/`, not the front-matter file.
