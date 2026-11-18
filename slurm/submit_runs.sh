#!/bin/bash
# Submits the two runs to SLURM:
#
#   binary   the black hole takes its star from a stellar binary (Code 0)  -> $OUT/binary
#   tidal    the black hole captures its star by tides (tc)                -> $OUT/tidal
#
#   slurm/submit_runs.sh             both
#   slurm/submit_runs.sh binary      one of them
#   slurm/submit_runs.sh tidal
#
# Build first: slurm/build.sh. What the runs are and what they ask of the
# cluster is set below; a value in the environment takes its place, as in
#
#   SBATCH_PARTITION=cpu CPUS=48 slurm/submit_runs.sh
#
# Run it again after a failure: what is finished is left alone.
set -euo pipefail
cd "$(dirname "$0")/.."

# ---- the runs --------------------------------------------------------------
PARAMS=${PARAMS:-data/params0_1000.csv}   # mass pairs and radii
OUT=${OUT:-out}                           # where the two directories go
NUM0=${NUM0:-100}                         # Code 0 scatterings per mass pair and radius
NUMA=${NUMA:-100}                         # Code A histories per black hole/star binary
NUMB=${NUMB:-1}                           # Code B scatterings per binary that reaches it
NUMC=${NUMC:-100}                         # Code C histories per black-hole binary

# ---- the cluster -----------------------------------------------------------
# export SBATCH_PARTITION=...  SBATCH_ACCOUNT=...  SBATCH_QOS=...   if yours asks for them
PARTS_BINARY=${PARTS_BINARY:-16}          # tasks of the job arrays of the binary run (a node each)
PARTS_TIDAL=${PARTS_TIDAL:-2}             # the tidal run is a hundred times smaller
export TIME=${TIME:-12:00:00}             # time limit of every job
export MEM=${MEM:-16G}                    # memory of the jobs that hold a whole table
# CPUS: cores of every task. Unset, it is the number every node of the partition has.
if [ -z "${CPUS:-}" ]; then
    CPUS=$(sinfo -h -o "%c" ${SBATCH_PARTITION:+-p "$SBATCH_PARTITION"} 2> /dev/null | sort -n | head -1 || true)
    [[ "$CPUS" =~ ^[0-9]+$ ]] || { echo "error: cannot ask sinfo for the cores of a node: set CPUS" >&2; exit 1; }
    echo "cores per task: $CPUS (what the nodes of ${SBATCH_PARTITION:-the cluster} have: set CPUS for another number)"
fi
export CPUS
echo "partition: ${SBATCH_PARTITION:-the default}; account: ${SBATCH_ACCOUNT:-your default}; time limit $TIME"

[ -x bin/scatter_spacehub ] || { echo "error: bin/scatter_spacehub is not there: run slurm/build.sh first" >&2; exit 1; }
which=${1:-both}
case "$which" in both|binary|tidal) ;; *) echo "error: binary, tidal or nothing for both, not '$which'" >&2; exit 1 ;; esac
numbers=(--nreal0 "$NUM0" --nrealA "$NUMA" --nrealB "$NUMB" --nrealC "$NUMC")

if [ "$which" != tidal ]; then
    echo "=== binary: $OUT/binary"
    PARTS=$PARTS_BINARY slurm/submit.sh "$OUT/binary" "$PARAMS" "${numbers[@]}"
fi
if [ "$which" != binary ]; then
    echo "=== tidal: $OUT/tidal"
    CAPTURE=tidal PARTS=$PARTS_TIDAL slurm/submit.sh "$OUT/tidal" "$PARAMS" "${numbers[@]}"
fi
echo
echo "squeue -u \$USER shows the jobs; each run ends with lists_spacehub.m in its directory"
