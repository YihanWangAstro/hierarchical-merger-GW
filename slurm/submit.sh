#!/bin/bash
# Submits the four stages to SLURM, the two long ones (Code A, Code C) as job
# arrays whose tasks run on a node each, and every stage waiting for the one
# before it.
#
#   slurm/submit.sh OUT PARAMS [options of the program]
#
#   slurm/submit.sh out/run data/params0_1000.csv --nreal0 200 --nrealA 5 --nrealB 8 --nrealC 20
#
# OUT is the directory of the tables, PARAMS the input table of Code 0.  What
# the jobs ask for is set in the environment:
#
#   PARTS=16          tasks of each array
#   CPUS=32           cores of every task
#   TIME=12:00:00     time limit of every job
#   MEM=16G           memory of the jobs that hold a whole table (Code 0, Code B, the merges)
#   BIN=bin/scatter_spacehub
#   CAPTURE=tidal     the stars are captured by tides: the tc stage in place of Code 0
#   SBATCH_PARTITION, SBATCH_ACCOUNT, SBATCH_QOS, ...   read by sbatch itself
#
# Run it again with the same arguments after a failure (a task that ran out of
# time, a node that went down): stages whose lists are there are left alone,
# and so are parts that were written to their end.
set -euo pipefail

if [ $# -lt 2 ]; then
    sed -n '2,23p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
OUT=$1; PARAMS=$2; shift 2
PARTS=${PARTS:-16}; CPUS=${CPUS:-32}; TIME=${TIME:-12:00:00}; MEM=${MEM:-16G}
BIN=${BIN:-bin/scatter_spacehub}
CAPTURE=${CAPTURE:-binary}
HERE=$(cd "$(dirname "$0")" && pwd)

[ -x "$BIN" ] || { echo "error: $BIN is not there: run make, or set BIN" >&2; exit 1; }
[ -r "$PARAMS" ] || { echo "error: cannot read $PARAMS" >&2; exit 1; }
mkdir -p "$OUT/log"
NAME=$("$BIN" --help | sed -n 's/.*lists_\(.*\)\.m.*/\1/p' | head -1)      # the integrator: names the files

after=""     # the job the next one waits for
# submit NAME-OF-JOB [sbatch options ...] -- command ...   (sets `after` to the new job)
submit() {
    local job=$1; shift
    local options=()
    while [ "$1" != "--" ]; do options+=("$1"); shift; done
    shift
    local id
    if ! id=$(sbatch --parsable ${after:+--dependency=afterok:$after} --job-name="bhs-$job" --time="$TIME" \
                     --output="$OUT/log/%x-%j.out" "${options[@]}" "$HERE/run.sbatch" "$@"); then
        refused "$job"
    fi
    after=${id%%;*}
    echo "submitted $job as job $after"
}
# sbatch said why it refused; say what the cluster would accept, as far as it tells.
refused() {
    {
        echo
        echo "error: sbatch refused the job $1 (its message is above); the jobs before it stay queued"
        echo "  asked for: partition ${SBATCH_PARTITION:-(the default)}, account ${SBATCH_ACCOUNT:-(your default)}," \
             "$CPUS cores, $TIME, $MEM for the jobs that hold a table"
        if command -v sacctmgr > /dev/null; then
            echo "  the accounts you may use, with the partitions and QOS they are tied to (empty: any):"
            echo "    account|partition|qos"
            sacctmgr -n -P show associations user="$USER" format=account,partition,qos 2> /dev/null | sort -u | sed 's/^/    /'
        fi
        if command -v sinfo > /dev/null; then
            echo "  the partitions, their time limit and the cores of their nodes:"
            sinfo -h -o "    %P %l %c" 2> /dev/null | sort -u
        fi
        echo "  set what fits and submit again, for example:"
        echo "    SBATCH_ACCOUNT=... SBATCH_PARTITION=... TIME=... slurm/submit_runs.sh"
    } >&2
    exit 1
}
done_already() { [ -s "$OUT/$1_$NAME.m" ] && echo "$1: its lists are there, left alone"; }

# Code 0 and Code B take minutes: one job each.
whole() {
    done_already "$1" || submit "$1" --cpus-per-task="$CPUS" --mem="$MEM" -- "$BIN" "$@" --dir "$OUT"
}
# Code A and Code C: an array of parts, then a small job that puts them together.
in_parts() {
    done_already "$1" && return
    submit "$1" --cpus-per-task="$CPUS" --array="1-$PARTS" --export="ALL,BHS_PARTS=$PARTS" \
           --output="$OUT/log/%x-%A_%a.out" -- "$BIN" "$@" --dir "$OUT"
    submit "$1-merge" --cpus-per-task=1 --mem="$MEM" -- "$BIN" "$@" --dir "$OUT" --merge "$PARTS"
}

if [ "$CAPTURE" = tidal ]; then
    # tc integrates nothing: one core, so that it need not wait for a whole node.
    done_already tc || submit tc --cpus-per-task=1 --mem=1G -- "$BIN" tc --in "$PARAMS" "$@" --dir "$OUT"
    in_parts codeA --in "$OUT/tc_$NAME.csv" "$@"
else
    whole    code0 --in "$PARAMS" "$@"
    in_parts codeA "$@"
fi
whole    codeB "$@"
in_parts codeC "$@"
submit   lists --cpus-per-task=1 --mem="$MEM" -- "$BIN" lists --dir "$OUT"
echo "the lists will be in $OUT/lists_$NAME.m; logs in $OUT/log; squeue -u \$USER shows the jobs"
