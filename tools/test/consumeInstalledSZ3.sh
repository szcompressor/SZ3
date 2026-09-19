#!/bin/bash
# Every shape a downstream project takes when it consumes an installed SZ3, asserted against an
# install tree.
#
#   tools/test/consumeInstalledSZ3.sh <install-prefix> [<extra cmake prefix path>]
#
# The second argument is appended to CMAKE_PREFIX_PATH, for an HDF5 that is not on the default
# search path. CMAKE_PREFIX_PATH from the environment is passed through as well.
#
# Assert which shape built and what the export said, never just that cmake exited 0. No set -e:
# every check below is explicit, and the total at the bottom catches a section that stopped early.
set -u
PREFIX=$(cd "$1" && pwd)
EXTRA_PREFIX=${2:-}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK" || exit 1

# Under MSYS2 and Cygwin, cmake is a native Windows program while this shell deals in /d/... paths
# it cannot open. A prefix handed over unconverted is simply not searched, and cmake then finds SZ3
# wherever else it can.
native() { if command -v cygpath > /dev/null 2>&1; then cygpath -m "$1"; else echo "$1"; fi; }
LIBDIR=$PREFIX/lib
[ -d "$LIBDIR" ] || LIBDIR=$PREFIX/lib64
CMPFX=$(native "$PREFIX")
[ -n "$EXTRA_PREFIX" ] && CMPFX="$CMPFX;$(native "$(cd "$EXTRA_PREFIX" && pwd)")"
[ -n "${CMAKE_PREFIX_PATH:-}" ] && CMPFX="$CMPFX;$CMAKE_PREFIX_PATH"
# Windows resolves a DLL through PATH and has no RPATH for an install tree to be recorded in, so
# <prefix>/bin -- where the install rule puts the runtime artifact -- is the only way the consumer
# built below reaches hdf5sz3. POSIX form: bash splits PATH on ':', and the MSYS2 runtime converts
# the whole variable when it spawns a native program.
export PATH="$PREFIX/bin:$PATH"

pass=0; fail=0
ok()  { echo "PASS  $1"; pass=$((pass+1)); }
bad() { echo "FAIL  $1"; shift; for l in "$@"; do echo "        $l"; done; fail=$((fail+1)); }
# want <name> <expected substring> <file>
want() { if grep -qF "$2" "$3"; then ok "$1"; else bad "$1" "expected to find: $2" "got:" "$(head -10 "$3")"; fi; }

echo "=== prefix $PREFIX ==="
echo "    cmake prefix path $CMPFX"

# ---------------------------------------------------------------- 1. C++, unversioned targets
# target_link_libraries(app PRIVATE SZ3 hdf5sz3): the spelling with no SZ3:: namespace, which CMake
# turns into a bare -lhdf5sz3 unless SZ3Config.cmake defines the alias.
mkdir -p cxx
cat > cxx/CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(consumer CXX)
find_package(SZ3 REQUIRED)
add_executable(app main.cpp)
target_link_libraries(app PRIVATE SZ3 hdf5sz3)
EOF
cat > cxx/main.cpp <<'EOF'
#include <H5Z_SZ3.hpp>
#include <SZ3/api/sz.hpp>
#include <cstdio>
int main() {
    printf("%s\n", H5Z_SZ3_initialize() < 0 ? "INIT FAILED" : "INIT OK");
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    SZ3::Config conf(100);
    conf.absErrorBound = 1e-3;
    set_SZ3_conf_to_H5(dcpl, conf);
    printf("NFILTERS %d\n", (int)H5Pget_nfilters(dcpl));
    return 0;
}
EOF
if cmake -S cxx -B cxx/b -DCMAKE_PREFIX_PATH="$CMPFX" > cxx/cfg.log 2>&1; then
    ok "cxx-consumer-configures"
    if cmake --build cxx/b --parallel 4 > cxx/build.log 2>&1; then
        ok "cxx-consumer-builds"
        if ./cxx/b/app > cxx/run.log 2>&1; then
            ok "cxx-consumer-runs"
        else
            bad "cxx-consumer-runs" "exit $?" "$(tail -5 cxx/run.log)"
        fi
    else
        bad "cxx-consumer-builds" "$(tail -15 cxx/build.log)"
        bad "cxx-consumer-runs" "nothing was built to run"
        : > cxx/run.log
    fi
else
    bad "cxx-consumer-configures" "$(tail -15 cxx/cfg.log)"
    bad "cxx-consumer-builds" "nothing was configured to build"
    bad "cxx-consumer-runs" "nothing was built to run"
    : > cxx/run.log
fi
want "cxx-consumer-initializes-the-filter"   "INIT OK"    cxx/run.log
want "cxx-consumer-sets-one-filter-on-a-dcpl" "NFILTERS 1" cxx/run.log

# ---------------------------------------------------------------- 2. C only, namespaced target
# A C project never enables CXX, so anything SZ3Config.cmake resolves per-language has to cope.
# REQUIRED, and no if (SZ3_FOUND) around the executable: written the other way this builds nothing
# whenever SZ3 is not found, which is exactly when a C-only consumer is broken.
mkdir -p conly
cat > conly/CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(conly C)
find_package(SZ3 REQUIRED)
add_executable(c1 main.c)
target_link_libraries(c1 PRIVATE SZ3::hdf5sz3)
EOF
printf 'int main(void){return 0;}\n' > conly/main.c
if cmake -S conly -B conly/b -DCMAKE_PREFIX_PATH="$CMPFX" > conly/cfg.log 2>&1; then
    ok "c-only-consumer-configures"
    if cmake --build conly/b --parallel 4 > conly/build.log 2>&1; then
        ok "c-only-consumer-builds"
    else
        bad "c-only-consumer-builds" "$(tail -15 conly/build.log)"
    fi
else
    bad "c-only-consumer-configures" "$(tail -15 conly/cfg.log)"
    bad "c-only-consumer-builds" "nothing was configured to build"
fi
if [ -x conly/b/c1 ] || [ -x conly/b/c1.exe ]; then
    ok "c-only-consumer-produces-an-executable"
else
    bad "c-only-consumer-produces-an-executable" "no c1 under conly/b"
fi

# ---------------------------------------------------------------- 3. the exported include paths
# A build directory or a third party's include tree named here reaches every consumer, and points
# at something their machine need not have.
grep -ho 'INTERFACE_INCLUDE_DIRECTORIES "[^"]*"' "$LIBDIR"/cmake/SZ3/*.cmake 2>/dev/null \
  | sed 's/.*"\(.*\)"/\1/' | tr ';' '\n' | sed '/^$/d' | sort -u > inc.log
# Counted first: with no include directories to look at, the check below has nothing to reject and
# would pass on an export that declares none.
if [ -s inc.log ]; then
    ok "export-declares-include-directories"
    sed 's/^/        /' inc.log
else
    bad "export-declares-include-directories" "no INTERFACE_INCLUDE_DIRECTORIES under $LIBDIR/cmake/SZ3"
fi
# grep -v, not a !-inverted grep: this has to report the offending line, and set -e would exempt it.
if grep -v '^\${_IMPORT_PREFIX}/' inc.log > outside.log; then
    bad "exported-include-dirs-stay-in-prefix" "outside the install prefix:" "$(cat outside.log)"
else
    ok "exported-include-dirs-stay-in-prefix"
fi

# ---------------------------------------------------------------- 4. no HDF5 on the machine
# SZ3Config.cmake reports not-found and returns, rather than ending the consumer's configure or
# handing them a link line they cannot use. GROMACS's GMX_USE_SZ3=AUTO depends on this shape.
mkdir -p probe
cat > probe/CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(probe CXX)
find_package(SZ3 QUIET)
if (SZ3_FOUND)
    message(FATAL_ERROR "SZ3 reported itself found with no HDF5 to be had")
endif ()
EOF
if cmake -S probe -B probe/b -DCMAKE_PREFIX_PATH="$CMPFX" \
        -DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON > probe/cfg.log 2>&1; then
    ok "sz3-declines-when-hdf5-cannot-be-found"
else
    bad "sz3-declines-when-hdf5-cannot-be-found" "$(tail -15 probe/cfg.log)"
fi

echo
echo "  $pass passed, $fail failed"

# Raise this with the check it comes with. A section that stops early otherwise shows only as a
# smaller number at the bottom that nobody compares.
EXPECTED=11
ran=$((pass + fail))
if [ "$ran" -ne "$EXPECTED" ]; then
    echo "  the suite accounted for $ran checks, not $EXPECTED"
    exit 1
fi
exit $((fail > 0))
