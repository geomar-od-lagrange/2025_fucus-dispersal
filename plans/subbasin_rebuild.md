# Rebuild after the shoreward subbasin assignment

`024a` now assigns `helcom_subbasin` by largest overlap with a shoreward
snap instead of centroid-within (rule and rationale:
[../docs/hexbinning_and_connectivity.md](../docs/hexbinning_and_connectivity.md)).
The change is committed; **nothing downstream has been rebuilt**.

**Folded into [horizons_and_weighting.md](horizons_and_weighting.md)
step 5**, which runs this rebuild on SLURM together with the 180 d
reducer rerun. The node-local recipe below is kept for SLURM outages.

## Everything below is STALE on NESH right now

| artefact | why |
|---|---|
| `HexAggregates/HexAgg_key_r6000m.parquet` + `.json` | written by the old centroid rule |
| `HexAgg_connectivity_r6000m_*` (24 partitions) | `origin/target_subbasin` from the old key |
| `HexAgg_survconn_r6000m_*` (452 partitions, 4 members) | ditto |
| `Figures/026a`, `026b`, `028`, `029`, `033` and their `Exports/` | consume the subbasin column |

`025`, `026`, `027`, `030`, `031` and the `024`/`024b`/`024d`/`024e`/`024f`
stores do **not** key on subbasin and stay valid. `032` reads
`Exports/028`, so it follows 028.

Until the rebuild runs, every subbasin-level number in the study still
describes the ~50 % of stranded weight and ~20–40 % of `n_obs` that the
old rule happened to name.

## Order and commands

SLURM was down when this was written, so the recipe is node-local:
`xargs -P` + `pixi run papermill`, one private `JUPYTER_RUNTIME_DIR` per
unit, three attempts (the ZMQ port race the job scripts document). With
SLURM back, prefer `scripts/submit_downstream_chain.sh`, which wires the
same ordering as a dependency chain.

Run everything from the repo root.

```sh
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

# 1. key (~5-10 min, single process). Second argument is the snap cap in
#    hex radii; 1.0 is production.
bash scripts/024a_BuildHexKey_job.sh 6000 1.0

# 2. 024c connectivity, 24 (regime, year) partitions, ~1-5 min each
for r in surface_stokes bottom; do for y in $(seq 2016 2025); do echo "$r $y"; done; done > /tmp/024c.txt
for y in $(seq 2016 2019); do echo "surface $y"; done >> /tmp/024c.txt
xargs -a /tmp/024c.txt -P 6 -n 1 bash -c '
  read -r regime year <<< "$1"
  export JUPYTER_RUNTIME_DIR="${TMPDIR:-/tmp}/jrt-024c-$$"; mkdir -p "$JUPYTER_RUNTIME_DIR"
  for a in 1 2 3; do bash scripts/024c_BuildHexConnectivity_job.sh "$regime" "$year" 6000 && exit 0
    sleep $((a*30)); done; echo "FAIL 024c $regime $year"' _

# 3. 024g survival connectivity, 113 (year, month) cells x 4 members on
#    disk. Measured ~270 s per cell under the SLURM fan-out, so budget
#    ~3 h at -P 12 for all four members; the production member alone
#    (step_wc0p1_ts3_tcinf_tf1_tw0) is ~45 min and unblocks 028.
#    The job script reads its rate model from the environment and
#    derives the member tag from it -- do not spell the tag by hand.
#    production member (TRAP_FLAT=1, TRAP_WALL=0):
SLURM_NTASKS=12 SLURM_CPUS_PER_TASK=1 \
  W_C=0.1 TAU_STRONG_HOURS=3.0 TAU_CALM_DAYS=0 DELTA=0 TRAP_FLAT=1 TRAP_WALL=0 \
  bash scripts/024g_BuildSurvivalConnectivity_job.sh surface_stokes 6000
#    the other three: TRAP_FLAT=1 TRAP_WALL=1 (tag step_wc0p1_ts3_tcinf),
#    TRAP_WALL=0.25, TRAP_WALL=0.5.
#    NOTE: that script srun's each cell. Off SLURM, drop the `srun ... --exact`
#    prefix from its xargs body (or copy the body into a node-local driver).

# 4. subbasin-keyed consumers. 026b is the long pole (~18 min per unit,
#    ~50 GiB RSS -> -P 2 on a 251 GB node); the rest are minutes.
#    026a/026b: {surface_stokes, bottom} x 5 seasons
#    028: surface_stokes unweighted + production member, bottom unweighted, x 5 seasons
#    029: surface_stokes x production member x 5 seasons
#    033: once
STAGES="026a 026b 028 029" LIST_ONLY=1 bash scripts/025-031_ConsumersAll_job.sh > /tmp/units.txt
#    then drive /tmp/units.txt with the same xargs pattern, calling the
#    per-stage scripts/0NN_*_job.sh directly (they are plain bash +
#    papermill; only the 025-031 driver needs SLURM).
bash scripts/033_ShorelineTypeMap_job.sh

# 5. headline report, which reads Exports/028
bash scripts/headline_report_job.sh
```

## Verify as you go

- `024a`: the validation cell prints the overlap/snap split, the residual
  `-1` class and the snap-cap sensitivity. The residual `-1` hexes that
  hold water must lie on the North Sea / Skagerrak side — if a Baltic
  coastal hex is still `-1`, the snap cap is too tight.
- `024c`: its own gate is `conn.n_obs.sum() == counts.n_obs.sum()`; the
  printed `n_obs fraction in -1 origin/target` is the number that must
  fall sharply against the old run.
- `028`: the "X % of n_obs dropped as -1" line. Old values, for the
  before/after record: 79.3 % (`surface_stokes`, ALL, unweighted),
  62.8 % (same, `step_wc0p1_ts3_tcinf_tf1_tw0`, `w_obs`), 60.3 %
  (`bottom`, SON, unweighted). Part of what remains is land-seeded
  particles (`release_hex == -1` in counts), which the subbasin rule
  cannot fix and 024g's population excludes by construction — so the
  weighted 028 should land much lower than the unweighted one.

## Then re-read the conclusions

Two products can move materially and must be re-read, not assumed:

- `029`'s per-subbasin source-vs-sink table — half the stranded weight
  was unattributed, so the transport direction it reports was drawn on a
  biased half of the deposit.
- `028`'s emission-fraction matrices and the self-retention numbers.
  Self-retention is a ratio within named basins; adding coastal hexes
  adds both diagonal and off-diagonal mass, so the direction of the
  change is not predictable from the dropped fraction alone.
