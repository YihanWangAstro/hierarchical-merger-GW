#!/bin/bash
# Submits the two runs to SLURM:
#
#   binary   the black hole takes its star from a stellar binary (Code 0)  -> $OUT/binary
#   tidal    the black hole captures its star by tides (tc)                -> $OUT/tidal
#
# The binary run is at the six radii of its input table, 0.11 to 2 pc; the
# tidal run takes the mass pairs of that table to all ten radii, 0.01 to 2 pc.
#
#   slurm/submit_runs.sh             both
#   slurm/submit_runs.sh binary      one of them
#   slurm/submit_runs.sh tidal
#
# Build first: slurm/build.sh. What the runs are and what they ask of the
# cluster is set below; a value in the environment takes its place, as in
#
#   CPUS=8 PARTS_BINARY=32 slurm/submit_runs.sh
#
# Run it again after a failure: what is finished is left alone.
set -euo pipefail
cd "$(dirname "$0")/.."

# ---- the runs --------------------------------------------------------------
PARAMS=${PARAMS:-data/params0_1000.csv}   # mass pairs, and the radii of the binary run
OUT=${OUT:-out}                           # where the two directories go
NUM0=${NUM0:-100}                         # Code 0 scatterings per mass pair and radius
NUMA=${NUMA:-100}                         # Code A histories per black hole/star binary
NUMB=${NUMB:-1}                           # Code B scatterings per binary that reaches it
NUMC=${NUMC:-100}                         # Code C histories per black-hole binary
# The tidal run has no Code 0 and so no NUM0: tides leave one black hole/star
# binary per mass pair and radius, always on the same orbit, where Code 0
# leaves some thirty. The size of its sample is the number of Code A
# histories of that one binary, NUM of the notebook "TC Code A". With 100 the
# run ends with some 200 black-hole binaries, too few: hence 100 times that.
NUMA_TIDAL=${NUMA_TIDAL:-10000}           # Code A histories per binary in the tidal run

# ---- the cluster -----------------------------------------------------------
export SBATCH_PARTITION=${SBATCH_PARTITION:-WorkC}   # the partition the jobs go to
# The account: SBATCH_ACCOUNT if it is set; otherwise the one SLURM ties you
# to on that partition, or to on any partition; failing that, your default.
# (SBATCH_QOS likewise, if your cluster asks for one.)
if [ -z "${SBATCH_ACCOUNT:-}" ] && command -v sacctmgr > /dev/null; then
    account=$(sacctmgr -n -P show associations user="$USER" format=account,partition 2> /dev/null |
              awk -F'|' -v p="$SBATCH_PARTITION" '$2 == p { print $1; found = 1; exit }
                                                  $2 == "" && any == "" { any = $1 }
                                                  END { if (!found && any != "") print any }' || true)
    [ -n "$account" ] && export SBATCH_ACCOUNT=$account
fi
CPUS=${CPUS:-64}                          # cores of every task; `node` for all the cores a node has
PARTS_BINARY=${PARTS_BINARY:-16}          # tasks of the job arrays of the binary run
PARTS_TIDAL=${PARTS_TIDAL:-16}            # and of the tidal run, which is about as large
# Time limit of every job. With 64 cores a task of the binary run is expected
# to take a quarter of an hour (estimated on other hardware, not measured
# here), and with 8 cores eight times that: six hours leave room for both. A
# job that runs into its limit is killed and its part has to be run again, so
# rather too long than too short; but the scheduler fits a short job into a
# gap more easily than a long one.
export TIME=${TIME:-06:00:00}
export MEM=${MEM:-16G}                    # memory of the jobs that hold a whole table
if [ "$CPUS" = node ]; then               # as many as every node of the partition has
    CPUS=$(sinfo -h -o "%c" -p "$SBATCH_PARTITION" 2> /dev/null | sort -n | head -1 || true)
    [[ "$CPUS" =~ ^[0-9]+$ ]] || { echo "error: cannot ask sinfo for the cores of a node: set CPUS to a number" >&2; exit 1; }
fi
[[ "$CPUS" =~ ^[0-9]+$ ]] || { echo "error: CPUS is a number of cores, or node: not '$CPUS'" >&2; exit 1; }
export CPUS
echo "cores per task: $CPUS"
echo "partition: $SBATCH_PARTITION; account: ${SBATCH_ACCOUNT:-your default}; time limit $TIME"

[ -x bin/scatter_spacehub ] || { echo "error: bin/scatter_spacehub is not there: run slurm/build.sh first" >&2; exit 1; }
# Built before the sources last changed (a git pull)? Then it is not the
# program these scripts go with: the tidal run would be at other radii.
if command -v make > /dev/null && ! make -q bin/scatter_spacehub 2> /dev/null; then
    echo "error: bin/scatter_spacehub is older than its sources: run slurm/build.sh again" >&2
    exit 1
fi
which=${1:-both}
case "$which" in both|binary|tidal) ;; *) echo "error: binary, tidal or nothing for both, not '$which'" >&2; exit 1 ;; esac

if [ "$which" != tidal ]; then
    echo "=== binary: $OUT/binary"
    PARTS=$PARTS_BINARY slurm/submit.sh "$OUT/binary" "$PARAMS" \
        --nreal0 "$NUM0" --nrealA "$NUMA" --nrealB "$NUMB" --nrealC "$NUMC"
fi
if [ "$which" != binary ]; then
    echo "=== tidal: $OUT/tidal"
    CAPTURE=tidal PARTS=$PARTS_TIDAL slurm/submit.sh "$OUT/tidal" "$PARAMS" \
        --nrealA "$NUMA_TIDAL" --nrealB "$NUMB" --nrealC "$NUMC"
fi
echo
echo "squeue -u \$USER shows the jobs; each run ends with lists_spacehub.m in its directory"
