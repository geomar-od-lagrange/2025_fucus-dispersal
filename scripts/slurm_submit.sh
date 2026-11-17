#!/bin/bash
# Sourceable helper: submit to a flaky slurmctld without ever duplicating a job.
#
#   source scripts/slurm_submit.sh
#   read -r id how <<< "$(robust_sbatch --job-name=chain-024d ... script.sh)"
#
# robust_sbatch prints one line, "<job id> <adopted|submitted>", on stdout and
# everything else on stderr -- two fields rather than a global, because callers
# capture it in $(...) and a variable set in that subshell never comes back.
#
# Why this exists: when the controller is loaded, `sbatch` can fail with
# "Socket timed out on send/recv operation" *after* the job has been created.
# A naive retry then queues the same work twice. So every submission carries a
# unique --job-name, and the name -- not the exit status -- is the source of
# truth: before each retry, and once before the first attempt, the queue is
# searched for that name. A name that is already there is adopted (its id read
# back) instead of resubmitted, which also makes a whole chain of these
# idempotent: rerun it after a partial failure and it picks up what is queued.
#
# The lookup itself can fail the same way, so it retries too -- and if it never
# succeeds, robust_sbatch gives up rather than guessing. A failed lookup means
# "unknown", and resubmitting on unknown is exactly the duplicate this avoids.
#
# Every call to sbatch and squeue is wrapped in `timeout`, because an
# unresponsive controller hangs the client rather than returning an error.
#
# Environment:
#   SBATCH_RETRY_DELAYS   backoff between attempts (default "15 30 60 120 120")
#                         -- one fewer entry than the number of attempts
#   SBATCH_CALL_TIMEOUT   seconds before an sbatch call is killed (default 120)
#   SQUEUE_CALL_TIMEOUT   seconds before a squeue call is killed (default 30)
#   SQUEUE_RETRIES        lookup attempts before "unknown" (default 6)

_robust_sbatch_seq=0

# Job id of the queued job named $1, empty if there is none.
# Returns non-zero only when squeue could not be reached at all -- an answer of
# "no such job" is a success with empty output.
slurm_id_for_name() {
    local name=$1 out attempt
    local retries=${SQUEUE_RETRIES:-6}
    for (( attempt = 1; attempt <= retries; attempt++ )); do
        if out=$(timeout "${SQUEUE_CALL_TIMEOUT:-30}" \
                 squeue --me -h -o '%i %j' -n "${name}" 2>/dev/null); then
            awk -v n="${name}" '$2 == n {print $1}' <<< "${out}" | sort | head -1
            return 0
        fi
        echo "  squeue lookup of ${name} failed (attempt ${attempt}/${retries})" >&2
        sleep 5
    done
    return 1
}

# robust_sbatch [sbatch args...] <script> [script args...]
# Prints "<job id> <adopted|submitted>" on stdout, everything else on stderr.
robust_sbatch() {
    local args=("$@") name="" a i
    for (( i = 0; i < ${#args[@]}; i++ )); do
        a=${args[i]}
        case "${a}" in
            --job-name=*) name=${a#--job-name=} ;;
            -J|--job-name) name=${args[i + 1]:-} ;;
        esac
    done
    if [ -z "${name}" ]; then
        # No caller-supplied name: synthesise one that is at least unique to
        # this invocation, so a timed-out-but-created job is still adoptable
        # within this call (across reruns it is not -- name your jobs).
        _robust_sbatch_seq=$(( _robust_sbatch_seq + 1 ))
        name="sbatch-$$-${_robust_sbatch_seq}"
        args=(--job-name="${name}" "${args[@]}")
    fi
    local has_parsable=""
    for a in "${args[@]}"; do
        [ "${a}" = "--parsable" ] && has_parsable=1
    done
    [ -n "${has_parsable}" ] || args=(--parsable "${args[@]}")

    local existing
    if ! existing=$(slurm_id_for_name "${name}"); then
        echo "cannot reach squeue: refusing to submit ${name} blind" >&2
        return 1
    fi
    if [ -n "${existing}" ]; then
        echo "adopted queued job ${existing} for ${name}" >&2
        echo "${existing} adopted"
        return 0
    fi

    local delays=(${SBATCH_RETRY_DELAYS:-15 30 60 120 120})
    local attempts=$(( ${#delays[@]} + 1 )) attempt out id rc
    for (( attempt = 1; attempt <= attempts; attempt++ )); do
        rc=0
        out=$(timeout "${SBATCH_CALL_TIMEOUT:-120}" sbatch "${args[@]}" 2>&1) || rc=$?
        # --parsable prints "<id>" or "<id>;<cluster>"; a plain sbatch prints
        # "Submitted batch job <id>". Take the first bare number either way.
        id=$(grep -oE '[0-9]+' <<< "${out}" | head -1)
        if [ "${rc}" -eq 0 ] && [ -n "${id}" ]; then
            echo "${id} submitted"
            return 0
        fi
        echo "  sbatch ${name} failed (attempt ${attempt}/${attempts}, rc=${rc}): ${out}" >&2
        # The attempt may still have created the job. Ask the queue, and adopt
        # what is there rather than submitting a second copy.
        if ! existing=$(slurm_id_for_name "${name}"); then
            echo "  sbatch failed and squeue is unreachable: ${name} is in an" \
                 "unknown state, not retrying" >&2
            return 1
        fi
        if [ -n "${existing}" ]; then
            echo "  the timed-out attempt landed after all: adopting ${existing}" >&2
            echo "${existing} adopted"
            return 0
        fi
        if [ "${attempt}" -lt "${attempts}" ]; then
            echo "  retrying ${name} in ${delays[attempt - 1]}s" >&2
            sleep "${delays[attempt - 1]}"
        fi
    done
    echo "gave up submitting ${name} after ${attempts} attempts" >&2
    return 1
}
