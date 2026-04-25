# 2025_fucus-dispersal

Parcels-based Lagrangian dispersal study for *Fucus vesiculosus* in the
Baltic Sea, driven by BSH operational-model currents (HBMnoku) and
optional Stokes drift from the CMEMS Baltic Wave Hindcast. Pipeline
covers input preparation, Parcels runs, per-trajectory visualisations,
and a hex-aggregated dispersal store for source↔sink queries.

## Repo layout

- `notebooks/` — pipeline as paired jupytext `.md` / `.ipynb` files,
  numbered `000` (preprocessing) through `025` (visualisations).
- `scripts/` — shell job scripts (`*_job.sh`), Python preprocessing
  scripts (`001`–`004`), and `obtain/` — the canonical recipes for
  fetching each input from upstream public sources.
- `data/` — git submodule pointing at the **data twin repo** at
  <https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>. The
  twin carries the curated input bundle (HELCOM polygons, Fucus
  REDLIST shapefile, BSH HBMnoku statics + demo subset, CMEMS Stokes
  sample). Throughout the codebase the term "twin" refers to this
  submodule.
- `docs/` — current-state documentation.
- `plans/` — open work plans (implemented plans move to
  `plans/done/`).

## Setup

Clone with the data twin attached:

```sh
git clone --recurse-submodules https://github.com/geomar-od-lagrange/2025_fucus-dispersal.git
cd 2025_fucus-dispersal
git -C data lfs pull
```

After a plain clone:

```sh
git submodule update --init data
git -C data lfs pull
```

`git-lfs` must be installed locally first (e.g. `pixi global install
git-lfs`), then `git lfs install` once per user. If the twin is
unreachable or you want to rebuild its blobs from upstream sources, run
`scripts/fetch_data.sh`.

Environment management uses **pixi**:

```sh
pixi install
```

Run all commands through `pixi run <command>` — see `AGENTS.md`.

## Licensing

Code is licensed under MIT (see `LICENSE`). The redistributed input
data carries per-dataset terms — see `ATTRIBUTION.md` for source,
catalog/DOI, and attribution requirements.
