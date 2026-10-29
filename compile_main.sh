#!/bin/sh
set -eu

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 Nchains Nevents_pre Nevents output_directory combination variables_file" >&2
    exit 2
fi

project_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build_dir="${GAMMADDBAR_BUILD_DIR:-$project_dir/build}"
mpi_compiler="${MPICXX:-mpicxx}"

# CMake handles quoted ROOT/BAT flags and links the selected MPI installation.
# Reuse a configured build directory to preserve ROOT_CONFIG_DIR/BAT_CONFIG_DIR.
cmake -S "$project_dir" -B "$build_dir" -DMPI_CXX_COMPILER="$mpi_compiler"
cmake --build "$build_dir" --target GammaDDbar

# set -e prevents an unsuccessful build from launching a stale executable.
cp "$build_dir/GammaDDbar" "$project_dir/main.x"
time "$project_dir/main.x" "$1" "$2" "$3" "$4" "$5" "$project_dir/Variables/$6"
