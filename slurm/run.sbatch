#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
# Runs the command it is given on one node, with the cores of --cpus-per-task
# (the program takes its number of threads from the job).  For example, the
# whole pipeline in one job:
#
#   sbatch --cpus-per-task=64 --time=04:00:00 slurm/run.sbatch \
#       bin/scatter_spacehub all --in data/params0_1000.csv --dir out/run \
#       --nreal0 200 --nrealA 5 --nrealB 8 --nrealC 20
#
# As a task of a job array it runs one part of a stage: slurm/submit.sh
# starts such arrays and puts their parts together.
set -euo pipefail
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    exec "$@" --part "${SLURM_ARRAY_TASK_ID}/${BHS_PARTS:?BHS_PARTS, the number of parts, is not set}"
fi
exec "$@"
