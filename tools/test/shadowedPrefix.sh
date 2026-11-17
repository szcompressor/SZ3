#!/bin/bash
# Another SZ3's headers on the include path of a build that is meant to use this one.
#
#   tools/test/shadowedPrefix.sh <hdf5 prefix> <older SZ3 source> <workdir> [scenario...]
#
# GROMACS #5821: HDF5 and sz3@3.2.0 shared one Spack view, GROMACS built its bundled SZ3 as a
# FetchContent(SYSTEM) subproject, and hdf5sz3 was compiled against the 3.2.0 headers because the
# view's include directory came first. Every scenario here puts an older SZ3 where a compiler can
# reach it ahead of the intended one, builds, and then asks the compiler's own dependency files
# which SZ3 headers each object was compiled from.
#
# The older SZ3 is installed twice: alone into <workdir>/old, which the cpath-* scenarios put on
# CPATH, and INTO <hdf5 prefix>, which the shared-* scenarios find HDF5 in. The second is removed
# again on exit, but a prefix other than a throwaway one must be opted into with
# SHADOW_MUTATE_PREFIX=1.
#
# A scenario passes when every object it checks was compiled from this tree's SZ3 headers. It
# checks hdf5sz3 always, and the parent's own target when SZ3 is a subproject that is not SYSTEM.
# Under SYSTEM, SZ3's headers reach the parent's target as -isystem, behind the HDF5 directory that
# target lists itself and behind CPATH: an order only the parent can change (docs/advanced-guide.md,
# "As a subproject"), so it is not asserted here. For the same reason there is no scenario for a
# SYSTEM subproject with the older SZ3 on CPATH, nor for a consumer of an installed SZ3.
#
# No set -e: every check is explicit, and the total at the bottom catches a section that stopped
# early. Unix Makefiles on purpose: Ninja consumes the dependency files this script reads.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=$(cd "$HERE/../.." && pwd)
HDF5_PREFIX=$(cd "$1" && pwd)
OLD_SRC=$(cd "$2" && pwd)
WORK=$3
shift 3
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
ALL="cpath-toplevel cpath-fetchcontent
shared-toplevel shared-toplevel-nosystem shared-add_subdirectory shared-add_subdirectory-system
shared-fetchcontent shared-fetchcontent-system"
SCENARIOS=${*:-$ALL}

pass=0; fail=0; skip=0
ok()  { echo "PASS  $1"; pass=$((pass+1)); }
bad() { echo "FAIL  $1"; shift; for l in "$@"; do echo "        $l"; done; fail=$((fail+1)); }
skipped() { echo "SKIP  $1"; skip=$((skip+1)); }

CUR_VER=$(sed -n 's/^project(SZ3 VERSION \([0-9.]*\)).*/\1/p' "$SRC/CMakeLists.txt")
OLD_VER=$(sed -n 's/^project(SZ3 VERSION \([0-9.]*\)).*/\1/p' "$OLD_SRC/CMakeLists.txt")
if [ -z "$CUR_VER" ] || [ -z "$OLD_VER" ] || [ "$CUR_VER" = "$OLD_VER" ]; then
    echo "need two different SZ3 versions: this tree is '$CUR_VER', $OLD_SRC is '$OLD_VER'"; exit 1
fi
case $HDF5_PREFIX in
"$WORK"/*) ;;
*) if [ "${SHADOW_MUTATE_PREFIX:-0}" != 1 ]; then
       echo "this installs SZ3 $OLD_VER into $HDF5_PREFIX; set SHADOW_MUTATE_PREFIX=1 to allow it"; exit 1
   fi ;;
esac
if [ -e "$HDF5_PREFIX/include/SZ3" ]; then
    echo "$HDF5_PREFIX already holds an SZ3, so it is not known which one the scenarios would meet"; exit 1
fi
CMAKE_VER=$(cmake --version | sed -n '1s/.* //p')
cmake_at_least() { [ "$(printf '%s\n%s\n' "$1" "$CMAKE_VER" | sort -V | head -1)" = "$1" ]; }

echo "=== shadowed prefix: SZ3 $CUR_VER ($SRC) against an older $OLD_VER, CMake $CMAKE_VER ==="

# ------------------------------------------------------------------ the older SZ3, built once
OLD=$WORK/old
rm -rf "$WORK/b-old" "$OLD"
if cmake -S "$OLD_SRC" -B "$WORK/b-old" -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
         -DBUILD_TESTING=OFF "-DCMAKE_INSTALL_PREFIX=$OLD" > "$WORK/old.log" 2>&1 \
   && cmake --build "$WORK/b-old" --parallel 4 >> "$WORK/old.log" 2>&1 \
   && cmake --install "$WORK/b-old" >> "$WORK/old.log" 2>&1 \
   && grep -q "define SZ3_VER  *\"$OLD_VER\"" "$OLD/include/SZ3/version.hpp"; then
    echo "    older SZ3 $OLD_VER installed alone in $OLD"
else
    echo "could not build and install the older SZ3, so no scenario has anything to be shadowed by:"
    tail -15 "$WORK/old.log" | sed 's/^/        /'
    exit 1
fi
SEEDED=0
seed() {
    [ "$SEEDED" = 1 ] && return 0
    if cmake --install "$WORK/b-old" --prefix "$HDF5_PREFIX" > "$WORK/seed.log" 2>&1 \
       && grep -q "define SZ3_VER  *\"$OLD_VER\"" "$HDF5_PREFIX/include/SZ3/version.hpp"; then
        # cmake --install --prefix rewrites the manifest, so it now lists what went into the prefix.
        cp "$WORK/b-old/install_manifest.txt" "$WORK/seed.manifest"
        SEEDED=1
        echo "    older SZ3 $OLD_VER installed into the HDF5 prefix $HDF5_PREFIX"
    else
        echo "could not install the older SZ3 into $HDF5_PREFIX:"; tail -5 "$WORK/seed.log"; exit 1
    fi
}
unseed() {
    [ "$SEEDED" = 1 ] || return 0
    tr -d '\r' < "$WORK/seed.manifest" | while read -r f; do rm -f "$f"; done
    rm -rf "$HDF5_PREFIX/include/SZ3" "$HDF5_PREFIX/include/SZ3c" "$HDF5_PREFIX"/lib*/cmake/SZ3
    rmdir "$HDF5_PREFIX"/lib*/cmake 2> /dev/null
    SEEDED=0
}
trap unseed EXIT

# ------------------------------------------------------------------ the parent project
# GROMACS's shape: it finds HDF5 itself, takes SZ3 as a subproject, and its own target lists the
# HDF5 include directory before it links the filter (src/gromacs/fileio/h5md/CMakeLists.txt).
mkdir -p "$WORK/parent"
cat > "$WORK/parent/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.19)
project(parent C CXX)
set(HOW "" CACHE STRING "how SZ3 is obtained")
set(SZ3_SRC "" CACHE PATH "SZ3 source tree")
find_package(HDF5 COMPONENTS C REQUIRED)
function(take_sz3)
    set(BUILD_H5Z_FILTER ON)
    set(BUILD_SZ3_BINARY OFF)
    set(BUILD_TESTING OFF)
    set(SZ3_INSTALL OFF)
    include(FetchContent)
    if(HOW STREQUAL "add_subdirectory")
        add_subdirectory(${SZ3_SRC} sz3 EXCLUDE_FROM_ALL)
    elseif(HOW STREQUAL "add_subdirectory-system")
        add_subdirectory(${SZ3_SRC} sz3 EXCLUDE_FROM_ALL SYSTEM)
    elseif(HOW STREQUAL "fetchcontent")
        FetchContent_Declare(SZ3 SOURCE_DIR ${SZ3_SRC})
        FetchContent_MakeAvailable(SZ3)
    elseif(HOW STREQUAL "fetchcontent-system")
        # What cmake/gmxManageSZ3.cmake does.
        FetchContent_Declare(SZ3 SOURCE_DIR ${SZ3_SRC} SYSTEM)
        FetchContent_MakeAvailable(SZ3)
    else()
        message(FATAL_ERROR "unknown HOW=${HOW}")
    endif()
endfunction()
take_sz3()
add_executable(app main.cpp)
target_include_directories(app SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(app PUBLIC ${HDF5_LIBRARIES})
target_link_libraries(app PRIVATE hdf5sz3)
EOF
cat > "$WORK/parent/main.cpp" <<'EOF'
#include <H5Z_SZ3.hpp>
#include <cstdio>
#include <cstring>
int main() {
    // What this file was compiled from, and what the filter library was.
    const H5Z_class2_t *cls = static_cast<const H5Z_class2_t *>(H5PLget_plugin_info());
    std::printf("app %s\nfilter %s\n", SZ3_VER, cls->name);
    SZ3::Config conf;
    std::printf("config %zu\n", sizeof(conf));
    return 0;
}
EOF

# Every SZ3 header the objects of <target> under <build> were compiled from, one per line.
headers_of() {
    find "$1" -path "*/$2.dir/*" \( -name '*.d' -o -name 'depend.make' \) -print0 2> /dev/null \
        | xargs -0 cat 2> /dev/null | tr ' \\' '\n\n' | grep '/SZ3/.*\.h' | sed 's/:$//' | sort -u
}
# check_headers <name> <build> <target> <intended root> <shadowing root>
check_headers() {
    headers_of "$2" "$3" > "$WORK/headers.txt"
    if [ ! -s "$WORK/headers.txt" ]; then
        bad "$1" "no dependency file under $2 names an SZ3 header, so nothing was checked"; return
    fi
    if grep -F "$5/" "$WORK/headers.txt" > "$WORK/shadowed.txt"; then
        bad "$1" "$(wc -l < "$WORK/shadowed.txt" | tr -d ' ') SZ3 headers came from the older SZ3, such as:" \
            "$(head -3 "$WORK/shadowed.txt")"
    elif grep -qF "$4/" "$WORK/headers.txt"; then
        ok "$1"
    else
        bad "$1" "the SZ3 headers came from neither tree:" "$(head -3 "$WORK/headers.txt")"
    fi
}
# The -I and -isystem directories of one object, in order, for the log.
include_line() {
    f=$(find "$1" -path "*/$2.dir/flags.make" | head -1)
    [ -n "$f" ] && sed -n 's/^CXX_INCLUDES = //p' "$f"
}

# run <scenario> <how> <shadow: cpath|shared> [cmake args...]
run() {
    name=$1; how=$2; shadow=$3; shift 3
    B=$WORK/b-$name
    rm -rf "$B"
    if [ "$shadow" = cpath ]; then
        unseed
        export CPATH=$OLD/include
        SHADOW=$OLD/include
    else
        seed
        unset CPATH
        SHADOW=$HDF5_PREFIX/include
    fi
    echo "--- $name"
    if [ "$how" = toplevel ]; then
        cmake -S "$SRC" -B "$B" -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DBUILD_H5Z_FILTER=ON \
              -DBUILD_SZ3_BINARY=OFF "-DCMAKE_PREFIX_PATH=$HDF5_PREFIX" "$@" > "$B.log" 2>&1
    else
        cmake -S "$WORK/parent" -B "$B" -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release "-DHOW=$how" \
              "-DSZ3_SRC=$SRC" "-DCMAKE_PREFIX_PATH=$HDF5_PREFIX" "$@" > "$B.log" 2>&1
    fi
    if [ $? != 0 ]; then
        bad "$name: configures" "$(tail -8 "$B.log")"; unset CPATH; return
    fi
    echo "        hdf5sz3: $(include_line "$B" hdf5sz3)"
    if cmake --build "$B" --target hdf5sz3 --parallel 4 > "$B.filter.log" 2>&1; then
        ok "$name: hdf5sz3 builds"
    else
        bad "$name: hdf5sz3 builds" "$(grep -m3 -E ' error[: ]' "$B.filter.log")"
    fi
    check_headers "$name: hdf5sz3 is compiled from this tree's SZ3 headers" "$B" hdf5sz3 "$SRC/include" "$SHADOW"
    case $how in
    toplevel | *-system) ;;  # no parent, or a parent whose include order is its own: see the top
    *)
        echo "        app:     $(include_line "$B" app)"
        if cmake --build "$B" --target app --parallel 4 > "$B.app.log" 2>&1; then
            check_headers "$name: the parent's own target is compiled from this tree's SZ3 headers" \
                          "$B" app "$SRC/include" "$SHADOW"
            "$B/app" > "$B.run.log" 2>&1
            if grep -qx "app $CUR_VER" "$B.run.log" && grep -q "^filter H5Z-SZ3-$CUR_VER " "$B.run.log"; then
                ok "$name: the parent's binary and the filter it links are both SZ3 $CUR_VER"
            else
                bad "$name: the parent's binary and the filter it links are both SZ3 $CUR_VER" "$(head -3 "$B.run.log")"
            fi
        else
            bad "$name: the parent's own target builds" "$(grep -m3 -E ' error[: ]' "$B.app.log")"
        fi ;;
    esac
    unset CPATH
}

for s in $SCENARIOS; do
    case $s in
    # ---- an older SZ3 on CPATH, as an environment module or `spack load` leaves it. CPATH is
    # searched after every -I and before every -isystem.
    cpath-toplevel)            run $s toplevel cpath ;;
    cpath-fetchcontent)        run $s fetchcontent cpath ;;
    # ---- an older SZ3 in the prefix HDF5 is found in: a Spack view, a conda environment, Homebrew.
    shared-toplevel)           run $s toplevel shared ;;
    # CMAKE_NO_SYSTEM_FROM_IMPORTED turns HDF5::HDF5's directory from -isystem into -I.
    shared-toplevel-nosystem)  run $s toplevel shared -DCMAKE_NO_SYSTEM_FROM_IMPORTED=ON ;;
    shared-add_subdirectory)   run $s add_subdirectory shared ;;
    shared-add_subdirectory-system)
        if cmake_at_least 3.25; then run $s add_subdirectory-system shared
        else skipped "$s (add_subdirectory SYSTEM needs CMake 3.25, this is $CMAKE_VER)"; fi ;;
    shared-fetchcontent)       run $s fetchcontent shared ;;
    # GROMACS #5821 itself.
    shared-fetchcontent-system)
        if cmake_at_least 3.25; then run $s fetchcontent-system shared
        else skipped "$s (FetchContent SYSTEM needs CMake 3.25, this is $CMAKE_VER)"; fi ;;
    *) bad "unknown scenario $s" "known: $(echo $ALL)" ;;
    esac
done

echo "== $pass passed, $fail failed, $skip skipped =="
[ "$fail" = 0 ]
