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
# SZ3Config.cmake decides whether SZ3::hdf5sz3 exists by looking for exactly this file, and
# BUILD_H5Z_FILTER is off by default -- so a prefix without it is sound, not broken.
HAVE_FILTER=0
[ -f "$LIBDIR/cmake/SZ3/HDF5SZ3.cmake" ] && HAVE_FILTER=1
# ... but only when nothing else in the prefix says the filter was built. tools/H5Z-SZ3 installs its
# headers and its binary through install() rules separate from install(EXPORT), so either one
# arriving without HDF5SZ3.cmake means the export was dropped rather than never built.
# Matched on *hdf5sz3*, not lib*: MSVC names them hdf5sz3.dll and hdf5sz3.lib, with no lib prefix.
if [ "$HAVE_FILTER" = 0 ]; then
    traces=$(ls -d "$PREFIX"/include/hdf5_sz3 2>/dev/null
             ls "$PREFIX"/lib*/*hdf5sz3* "$PREFIX"/lib*/plugin/*hdf5sz3* "$PREFIX"/bin/*hdf5sz3* 2>/dev/null)
    if [ -n "$traces" ]; then
        echo "FAIL  the install carries the filter but the export does not declare it"
        echo "        missing: $LIBDIR/cmake/SZ3/HDF5SZ3.cmake"
        echo "        present:"
        echo "$traces" | sed 's/^/          /'
        echo "        SZ3::hdf5sz3 does not exist for any consumer, and every filter check below"
        echo "        would have been skipped as though this SZ3 had been built without the filter."
        exit 1
    fi
fi
CMPFX=$(native "$PREFIX")
[ -n "$EXTRA_PREFIX" ] && CMPFX="$CMPFX;$(native "$(cd "$EXTRA_PREFIX" && pwd)")"
[ -n "${CMAKE_PREFIX_PATH:-}" ] && CMPFX="$CMPFX;$CMAKE_PREFIX_PATH"
# Windows resolves a DLL through PATH and has no RPATH for an install tree to be recorded in, so
# <prefix>/bin -- where the install rule puts the runtime artifact -- is the only way the consumer
# built below reaches hdf5sz3. POSIX form: bash splits PATH on ':', and the MSYS2 runtime converts
# the whole variable when it spawns a native program.
export PATH="$PREFIX/bin:$PATH"

pass=0; fail=0; skip=0
ok()  { echo "PASS  $1"; pass=$((pass+1)); }
bad() { echo "FAIL  $1"; shift; for l in "$@"; do echo "        $l"; done; fail=$((fail+1)); }
skipped() { echo "SKIP  $1"; skip=$((skip+1)); }
# reason <text>, then the names it covers
skip_all() { reason=$1; shift; for n in "$@"; do skipped "$n ($reason)"; done; }
# want <name> <expected substring> <file>
want() { if grep -qF "$2" "$3"; then ok "$1"; else bad "$1" "expected to find: $2" "got:" "$(head -10 "$3")"; fi; }

echo "=== prefix $PREFIX ==="
echo "    cmake prefix path $CMPFX"

# ---------------------------------------------------------------- 1. C++, unversioned targets
# target_link_libraries(app PRIVATE SZ3 hdf5sz3): the spelling with no SZ3:: namespace, which CMake
# turns into a bare -lhdf5sz3 unless SZ3Config.cmake defines the alias.
if [ "$HAVE_FILTER" = 1 ]; then
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
    printf("%s\n", H5Zregister(H5PLget_plugin_info()) < 0 ? "INIT FAILED" : "INIT OK");
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    SZ3::Config conf(100);
    conf.absErrorBound = 1e-3;
    set_SZ3_conf_to_H5(dcpl, conf);
    printf("NFILTERS %d\n", (int)H5Pget_nfilters(dcpl));
    return 0;
}
EOF
# A multi-config generator defaults an unqualified build to Debug and writes it to b/<config>, so
# the checks below look for a binary that is not there. Naming the config and its output directory
# puts it at b/ either way: with no CMAKE_BUILD_TYPE here, a single-config generator ignores both.
if cmake -S cxx -B cxx/b -DCMAKE_PREFIX_PATH="$CMPFX" \
        -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE="$(native "$WORK")/cxx/b" > cxx/cfg.log 2>&1; then
    ok "cxx-consumer-configures"
    if cmake --build cxx/b --config Release --parallel 4 > cxx/build.log 2>&1; then
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
else
skip_all "this SZ3 was built without BUILD_H5Z_FILTER, so the export carries no hdf5sz3" \
    cxx-consumer-configures cxx-consumer-builds cxx-consumer-runs \
    cxx-consumer-initializes-the-filter cxx-consumer-sets-one-filter-on-a-dcpl
fi

# ---------------------------------------------------------------- 2. C only, namespaced target
# A C project never enables CXX, so anything SZ3Config.cmake resolves per-language has to cope.
# REQUIRED, and no if (SZ3_FOUND) around the executable: written the other way this builds nothing
# whenever SZ3 is not found, which is exactly when a C-only consumer is broken.
if [ "$HAVE_FILTER" = 1 ]; then
mkdir -p conly
cat > conly/CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(conly C)
find_package(SZ3 REQUIRED)
add_executable(c1 main.c)
target_link_libraries(c1 PRIVATE SZ3::hdf5sz3)
EOF
printf 'int main(void){return 0;}\n' > conly/main.c
if cmake -S conly -B conly/b -DCMAKE_PREFIX_PATH="$CMPFX" \
        -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE="$(native "$WORK")/conly/b" > conly/cfg.log 2>&1; then
    ok "c-only-consumer-configures"
    if cmake --build conly/b --config Release --parallel 4 > conly/build.log 2>&1; then
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
else
skip_all "this SZ3 was built without BUILD_H5Z_FILTER, so the export carries no hdf5sz3" \
    c-only-consumer-configures c-only-consumer-builds c-only-consumer-produces-an-executable
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
# grep -v, not a !-inverted grep: this has to report the offending line.
if grep -v '^\${_IMPORT_PREFIX}/' inc.log > outside.log; then
    bad "exported-include-dirs-stay-in-prefix" "outside the install prefix:" "$(cat outside.log)"
else
    ok "exported-include-dirs-stay-in-prefix"
fi

# ---------------------------------------------------------------- 4. no HDF5 on the machine
# SZ3Config.cmake reports not-found and returns, rather than ending the consumer's configure or
# handing them a link line they cannot use. GROMACS's GMX_USE_SZ3=AUTO depends on this shape.
# Only the filter makes SZ3 depend on HDF5 at all, so without it there is no declining to observe.
if [ "$HAVE_FILTER" = 1 ]; then
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
else
skip_all "this SZ3 was built without BUILD_H5Z_FILTER, so it depends on no HDF5" \
    sz3-declines-when-hdf5-cannot-be-found
fi

# ---------------------------------------------------------------- 5. no Zstd on the machine
# The other way SZ3Config.cmake declines. It resolves a system Zstd with find_library, which no
# CMAKE_DISABLE_FIND_PACKAGE_ reaches, so point the library search at an empty root instead.
# A bundled-Zstd build carries SZ3::zstd in the export and never runs that branch.
if grep -q "TARGET SZ3::zstd" "$LIBDIR/cmake/SZ3/SZ3Targets.cmake" 2>/dev/null ||
   [ -n "$(find "$LIBDIR" -name '*sz3_zstd*' 2>/dev/null)" ]; then
    skip_all "this SZ3 bundles its own Zstd, so it looks for none here" \
        sz3-declines-when-zstd-cannot-be-found
else
    mkdir -p zprobe/empty
    cp probe/CMakeLists.txt zprobe/CMakeLists.txt 2>/dev/null || cat > zprobe/CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(probe CXX)
find_package(SZ3 QUIET)
if (SZ3_FOUND)
    message(FATAL_ERROR "SZ3 reported itself found with no Zstd to be had")
endif ()
EOF
    if cmake -S zprobe -B zprobe/b -DCMAKE_PREFIX_PATH="$CMPFX" \
            -DCMAKE_FIND_ROOT_PATH="$(native "$WORK")/zprobe/empty" \
            -DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY > zprobe/cfg.log 2>&1; then
        ok "sz3-declines-when-zstd-cannot-be-found"
    else
        bad "sz3-declines-when-zstd-cannot-be-found" "$(tail -15 zprobe/cfg.log)"
    fi
fi

# ---------------------------------------------------------------- 6. no SZ3 internals exported
# Another SZ3 copy in the same process would take over any SZ3 function the filter exports.
if [ -f "$LIBDIR/libhdf5sz3.so" ] && command -v nm > /dev/null 2>&1; then
    nm -DC --defined-only "$LIBDIR/libhdf5sz3.so" | grep -E "^[0-9a-f]+ [TWVu] SZ3::" > exports.log
    if [ -s exports.log ]; then
        bad "hdf5sz3-exports-no-sz3-functions" "$(wc -l < exports.log) exported, for example:" "$(head -3 exports.log)"
    else
        ok "hdf5sz3-exports-no-sz3-functions"
    fi
else
    skipped "hdf5sz3-exports-no-sz3-functions (no shared ELF libhdf5sz3 here)"
fi

echo
echo "  $pass passed, $fail failed, $skip skipped"

# Raise this with the check it comes with. A section that stops early otherwise shows only as a
# smaller number at the bottom that nobody compares.
EXPECTED=13
ran=$((pass + fail + skip))
if [ "$ran" -ne "$EXPECTED" ]; then
    echo "  the suite accounted for $ran checks, not $EXPECTED"
    exit 1
fi
exit $((fail > 0))
