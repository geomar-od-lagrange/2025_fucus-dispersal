# 2025_fucus-dispersal

Parcels-based Lagrangian dispersal study for *Fucus vesiculosus* in the
Baltic Sea, driven by BSH operational-model currents (HBMnoku) plus
optional Stokes drift from CMEMS wave hindcasts. Pipeline covers input
preparation, Parcels runs, per-trajectory visualisations, and a
hex-aggregated dispersal store for source↔sink queries.

`./data/` is a git submodule pointing at the **data twin** repo at
<https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>; "twin"
throughout the codebase means this submodule. Agent and contributor
conventions live in `AGENTS.md`; methodology in `docs/`; open work in
`plans/`.

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
