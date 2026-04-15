# Job Script — Action Plan

## Fixes for this PR

### Fix broken notebook path

`notebooks/FucusDispersal.ipynb` → `notebooks/010_FucusDispersal.ipynb`.

### Remove dead `mkdir` for output/Trajectories

The notebook ensures the output dir exists via absolute HPC paths. The job script's `mkdir -p output/Trajectories/${release_year}/` is unused.

### Remove redundant `mkdir -p notebooks_executed/`

The next line creates `notebooks_executed/TrajectoryCalc/${release_year}` with `-p`, which creates all parents.

### Make job name generic

Replace `19h1d0s087Fucus_dispersal` with `Fucus_dispersal`.

### Remove `srun &` + `wait` scaffolding

Single `srun` backgrounded then immediately `wait`ed — dead scaffolding. Remove the `&` and `wait`.

### Boolean params are correct — no change needed

Tested: papermill `-p` uses `ast.literal_eval`. Python-style `True`/`False` works. Lowercase `true`/`false` would break (passed as truthy strings).

### 120G per CPU is correct — no change needed

Actual usage for high-res data + trajectories. Scale out across multi-tenant nodes.

## Restructure for parallel execution (next step)

Deferred until after notebook performance improvements. The pattern:

```bash
#!/bin/bash
#SBATCH --job-name=Fucus_dispersal
#SBATCH --ntasks=N_PARALLEL
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=48:00:00
#SBATCH --partition=base

module load gcc12-env/12.3.0
module load singularity/3.11.5

N_PARALLEL="${SLURM_NTASKS}"
NOTEBOOK="notebooks/010_FucusDispersal.ipynb"
CONTAINER="parcels-container_2024.10.07-7af7fd0.sif"

# Fixed parameters
MAX_AGE_DAYS=220
CALC_DT_MINS=15
OUTPUT_DT_MINS=60
PARTICLES_PER_CELL=10
REPEATED_RELEASE=True
REPEATED_RELEASE_DT_DAYS=7

# Generate one srun command per parameter combination
generate_commands() {
    for release_year in 2018 2019 2020; do
        for half in 1 2; do
            if [ "$half" -eq 1 ]; then
                first_month=1; first_day=1; last_month=6; last_day=25
            else
                first_month=7; first_day=1; last_month=12; last_day=30
            fi
            for release_depth in 0; do
                for speed in 1 0.97 0.87; do
                    tag="y${release_year}_h${half}_d${release_depth}_s${speed}"
                    outdir="notebooks_executed/TrajectoryCalc/${release_year}"
                    mkdir -p "${outdir}"

                    echo "srun --ntasks=1 --exact \
                        singularity run \
                            -B /sfs -B /gxfs_work -B \$PWD:/work \
                            --pwd /work \
                            ${CONTAINER} \
                            bash -c \". /opt/conda/etc/profile.d/conda.sh \
                                && conda activate base \
                                && papermill --cwd notebooks/ \
                                    ${NOTEBOOK} \
                                    ${outdir}/Fucus_${tag}.ipynb \
                                    -p release_year ${release_year} \
                                    -p first_release_month ${first_month} \
                                    -p first_release_day ${first_day} \
                                    -p last_release_month ${last_month} \
                                    -p last_release_day ${last_day} \
                                    -p max_age_days ${MAX_AGE_DAYS} \
                                    -p calc_dt_mins ${CALC_DT_MINS} \
                                    -p output_dt_mins ${OUTPUT_DT_MINS} \
                                    -p relative_release_depth ${release_depth} \
                                    -p relative_particle_speed ${speed} \
                                    -p particles_per_cell ${PARTICLES_PER_CELL} \
                                    -p repeated_release ${REPEATED_RELEASE} \
                                    -p repeated_release_dt_days ${REPEATED_RELEASE_DT_DAYS} \
                                    -p is_papermill True \
                                    -k python\""
                done
            done
        done
    done
}

generate_commands | xargs -P "${N_PARALLEL}" -I{} bash -c '{}'

jobinfo
```

### SLURM directives for parallel mode

| Directive | Old | New | Reason |
|---|---|---|---|
| `--ntasks` | 1 | N_PARALLEL | One task slot per concurrent `srun --exact` |
| `--cpus-per-task` | 1 | 1 (keep) | Each task is single-CPU |
| `--mem-per-cpu` | 120G | 120G (keep) | Actual usage |
| `--time` | 48:00:00 | Adjust down | Wall time = ceil(scenarios / N) * single_task_time |

Use `srun --exact` (not `--exclusive`) so each step claims exactly 1 task slot. Needs testing on the cluster.

### Scenario parameters

Leave for later — notebook performance work comes first, parameterization follows naturally. Bottom roughness scenarios will need a new notebook parameter.
