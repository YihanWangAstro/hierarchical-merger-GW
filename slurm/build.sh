#!/bin/bash
# Builds the programs on a cluster.
#
#   slurm/build.sh                 with the g++ that is there
#   slurm/build.sh gcc/12.2.0      after `module load gcc/12.2.0`; `module avail gcc` lists what there is
#
# SpaceHub is taken from $SPACEHUB (its src directory), or from ../SpaceHub,
# which is cloned from GitHub if it is not there. A module that was loaded is
# written to slurm/modules.sh, and the jobs load it again before they run.
set -eo pipefail
cd "$(dirname "$0")/.."

# `module` is a shell function that a script does not always inherit.
find_module_command() {
    type module > /dev/null 2>&1 && return 0
    local init
    for init in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /usr/share/Modules/init/bash \
                /usr/share/lmod/lmod/init/bash /opt/ohpc/admin/lmod/lmod/init/bash; do
        [ -r "$init" ] && . "$init" && type module > /dev/null 2>&1 && return 0
    done
    return 1
}

: > slurm/modules.sh
for name in "$@"; do
    find_module_command || { echo "error: there is no \`module\` command here to load $name with" >&2; exit 1; }
    module load "$name" || { echo "error: cannot load the module $name: \`module avail\` lists what there is" >&2; exit 1; }
    echo "module load $name" >> slurm/modules.sh
done

CXX=${CXX:-g++}
command -v "$CXX" > /dev/null || { echo "error: no compiler $CXX: load a module, as in slurm/build.sh gcc/12.2.0" >&2; exit 1; }
echo "compiler: $("$CXX" --version | head -1)"

if [ -z "${SPACEHUB:-}" ]; then
    if [ ! -d ../SpaceHub/src ]; then
        echo "SpaceHub is not in ../SpaceHub: cloning it"
        git clone https://github.com/YihanWangAstro/SpaceHub.git ../SpaceHub
    fi
    SPACEHUB=$(cd ../SpaceHub/src && pwd)
fi
[ -r "$SPACEHUB/spaceHub.hpp" ] || { echo "error: no spaceHub.hpp in $SPACEHUB" >&2; exit 1; }
# Before the commit 1a0ca3fd one encounter in ten million never ends.
if git -C "$SPACEHUB" rev-parse > /dev/null 2>&1 && ! git -C "$SPACEHUB" merge-base --is-ancestor 1a0ca3fd HEAD 2> /dev/null; then
    echo "warning: SpaceHub in $SPACEHUB does not hold the commit 1a0ca3fd: update it (git pull)" >&2
fi
echo "SpaceHub: $SPACEHUB"

if ! make CXX="$CXX" SPACEHUB="$SPACEHUB"; then
    echo >&2
    echo "The build failed. If the compiler is old (C++17 is needed; GCC 9 or later is safest)," >&2
    echo "load a newer one: \`module avail gcc\`, then slurm/build.sh <that module>." >&2
    exit 1
fi

# A dozen scatterings: do the programs run here?
check=$(mktemp -d)
if bin/scatter_spacehub code0 --nreal 1 --nsys 1 --threads 1 --dir "$check" > "$check/log" 2>&1; then
    echo "built bin/scatter_spacehub, bin/scatter_bs and bin/scatter_leapfrog; a test of twelve scatterings ran"
    rm -rf "$check"
else
    echo "error: the programs were built but do not run here:" >&2
    tail -5 "$check/log" >&2
    exit 1
fi
[ -s slurm/modules.sh ] && echo "the jobs will load: $(sed 's/module load //' slurm/modules.sh | tr '\n' ' ')"
echo "next: slurm/submit_runs.sh"
