# Run extension 2020–2025: restart handoff

Release years **2020–2025** are being simulated for the two headline regimes,
`bottom` and `surface_stokes` — 73 releases per full year and 29 for 2025
(doys 1…141; the last 220 d window that the 2D fields cover closes
2025-12-27), so **394 releases per regime, 788 in total**. The Stokes forcing
was extended past the CMEMS hindcast with the NRT product (`002`), and eight
`surface_stokes` 2D-field stems missing on **2021-03-04** (fine + coarse × the
four 6-hourly slots) were rebuilt by `003` job 23836382 before the 2021 run
could start. Everything downstream — the `024`/`024b`/`024c` aggregates, the
`024d` sidecars, the `024e/f/g` reducers, the `025…031` consumers and the
`032` headline — is already queued as one SLURM dependency chain
(`scripts/submit_downstream_chain.sh`, tag `fxea1a333f`, chain log
`…/logs/chain/chain_20260915T134040.log`, git `d02897e`).

## The job graph, as read from `scontrol` at 2026-09-15T16:06 CEST

| stage | job(s) | what | state | dependency |
|---|---|---|---|---|
| 010 | 23836364 / 23836392 / 23836393 / 23836394 / 23836401 / 23836403 | `bottom` 2020 / 21 / 22 / 23 / 24 / 25 | RUNNING since 12:43–13:33, wall limit 1-12:00:00 | none |
| 010 | 23836395 / 23836459 / 23836397 / 23836398 / 23836402 / 23836404 | `surface_stokes` 2020 / 21 / 22 / 23 / 24 / 25 | RUNNING since 13:03–14:00, wall limit 1-12:00:00 | none |
| A | 23836462…23836497 (36 jobs) | `024` counts, `024b` distance, `024c` connectivity — one triple per (regime, year) | PENDING | each triple `afterany` its own 010 job; the `024c` of each triple `afterok` its `024` |
| B | 23836498 | `024d` beaching sidecars, `surface_stokes` only | PENDING | `afterany` all six `surface_stokes` 010 jobs |
| C | 23836499 | `024e`/`024f`/`024g` reducers, every member | PENDING | `afterok:23836498` |
| D | 23836533 | `025…031` consumers, 80 units / 12 tasks | PENDING | `afterok` all 36 stage-A jobs **and** 23836499 |
| E | 23836534 | `032` headline regeneration | PENDING | `afterok:23836533` |

`--ntasks` on B and C (292 each) is a **concurrency cap, not a work count**:
both jobs rebuild their work list from what is on disk at run time, so B will
process 686 zarrs (292 existing + 394 new) and C 2 825 cells, `xargs -P` many
rounds deep.

## Has a stage landed?

Counts against `output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs`:

| stage | check | expected |
|---|---|---|
| 010 | `ls -d $O/Trajectories/<regime>/<year>/*.zarr \| wc -l` | 73 per regime-year, **29 for 2025** |
| A | `ls $O/HexAggregates/HexAgg_{counts,distance,connectivity}_r6000m_*` | 24 each (12 now + 12 new; one parquet per regime-year, `surface` included) |
| B | `ls -d $O/BeachingForcing/surface_stokes/*/*.zarr \| wc -l` | **686** |
| C | `HexAgg_beaching_*` 1 582 (14 members), `HexAgg_survocc_*` 791 (7), `HexAgg_survconn_*` 452 (4) | 113 monthly partitions per member: 48 + 60 + 5 (2025 has only Jan–May) |
| D | `$O/Figures/{025,026,026a,026b,027,028,029,030,031}`, `$O/Exports/{025,026,026a,026b,027,028,029}` | all present, PNG mtimes after D ran |
| E | `$O/Headline/HEADLINE.md` + the eleven numbered PNGs | mtimes after E ran |

Logs: the 010 runs write `slurm-<jobid>.out` into the **repo root**; every
chained stage writes under `$O/logs/` (`aggregates/`, `024d/`, `024efg/`,
`025-031/`, `headline/`, `chain/`). Each multi-task stage ends with a tally
line — `024d done: <nok> OK, <nfail> FAIL of <njobs>`, and the same shape for
`024efg done:` and `025-031 done:`. A non-zero FAIL count names the failed
units on the preceding `FAIL …` lines.

## If a stage failed

The chain is **resumable and idempotent**. Rerun

```
scripts/submit_downstream_chain.sh 23836364 23836392 23836393 23836394 \
    23836401 23836403 23836395 23836459 23836397 23836398 23836402 23836404
```

and it adopts every stage still queued under its own job name instead of
duplicating it; `ADOPT_A` / `ADOPT_B` / `ADOPT_C` take ids for stages queued
under names it cannot reconstruct (which is how the current graph was
assembled — stages A–C were adopted, D and E freshly submitted). A stage that
has already *finished* is no longer in the queue and would be resubmitted, so
resume while the chain is still pending. Every stage script honours
`LIST_ONLY=1` to print its work list without touching the queue, and
`scripts/submit_010.sh`, `submit_aggregates.sh`, `submit_024d.sh`,
`submit_024efg.sh` re-submit a single unit.

Failure modes seen on 2026-09-15:

- **Slurm "Duplicate jobid" controller fault.** Killed 23836396
  (`surface_stokes` 2021) six minutes in, mid-setup; `sacct` reports only
  `CANCELLED 0:15`, no message in the job log. The culprit nodes are in
  `scripts/node-blacklist.txt` (clk352, clk407, clk423, clk468, clk584,
  srp205, srp243, srp266); 23836459 is the resubmission. **The blacklist is
  not applied to the twelve running 010 jobs** — several are allocated on
  listed nodes — so a repeat is possible.
- **Socket timeouts where the job was created anyway.** `sbatch` returns
  "Socket timed out on send/recv operation" after the controller has already
  made the job. `robust_sbatch` (`scripts/slurm_submit.sh`) re-queries the
  queue by job name and adopts it rather than submitting a second copy.
- **srun CPU-bind launch failure.** Every per-release step must be launched
  `srun -n1 -N1 --exact --no-kill …`; without `-N1` the step inherits the
  job's node count and dies on the CPU binding.

## Data caveats

The 2025 forcing is not hole-free and the cause is upstream: **25 raw BSH
c-files in the source archive are corrupt** (HDF read errors) and need a
re-delivery from BSH. The resulting hole table, the hole-free release
boundary and which release years interpolate across holes are in
[../docs/job_scripts.md](../docs/job_scripts.md) §"The Parcels runs (010)" —
read it there rather than duplicating it. Short version: hole-free to
2025-01-12T12:00, so **every 2025 release and 41 of the 2024 releases**
interpolate across scattered 6 h gaps.

Second caveat for the consumers: 2025 stops at doy 141, so a pooled season
draws on six release years for DJF/MAM but only 2020–2024 (plus 2016–2019)
for JJA/SON. Say which years a seasonal panel actually pools.

## Quota

At writing: `$WORK` **26.441 T used** of 29 T soft / 30 T hard (10.69 M
files). Measured on disk, a full regime-year costs 205 G (`bottom`) and 110 G
(`surface_stokes`), and a year of `024d` sidecars 8.1 G, so the sweep adds
≈ **1.7 T** (1.1 T bottom + 0.6 T surface_stokes + 44 G sidecars) and lands
around 28.2 T — inside the soft limit with under a terabyte of slack.

**Stop line:** do not start another full-regime sweep, a second `024d` window
or a re-run of an existing regime-year without freeing space first. The
obvious candidate is `Trajectories/surface/2016…2019`: the `surface` (no
Stokes) regime is not extended and no consumer reads it
(`025-031_ConsumersAll_job.sh` defaults to `surface_stokes bottom`).

## Open items

- **`surface` stays 2016–2019.** Its `024`/`024b`/`024c` stores keep the old
  year range while the two headline regimes go to 2025. No consumer reads
  `surface`, so this is a store-inventory asymmetry, not a figure one.
- **Beaching is `surface_stokes`-only**, so `029`/`030`/`031` and the
  survival-weighted `028` cover one regime whatever the year range.
- [analysis_wrapup.md](analysis_wrapup.md) still lists §4 (biomass maps,
  blocked on the Estonian production figures) and §7 (bottom-velocity sweep,
  deferred — it needs a kernel change plus fresh runs).
- ~~Refresh the trap sensitivity.~~ Done: the trap-axis table in
  [../docs/beaching.md](../docs/beaching.md) matches the 2016–2025 pooled
  Headline table B (0.769 / 0.725 / 0.668 / 0.416).
