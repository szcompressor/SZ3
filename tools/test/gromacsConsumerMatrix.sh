#!/bin/bash
# Every value of GROMACS's GMX_USE_SZ3, asserted against an installed SZ3.
#
#   tools/test/gromacsConsumerMatrix.sh <gromacs source> <workdir> <sz3 prefix|-> <hdf5 prefix> \
#                                       <found|absent|declined> [extra cmake args...]
#
# The fifth argument says what the installed SZ3 does here: get found, not be installed at all, or
# decline because a dependency of its own cannot be resolved. Declining and being absent look the
# same to GROMACS, which is the point -- an SZ3 that ends the configure instead of declining
# breaks GMX_USE_SZ3=AUTO, whose contract is the fallback. They must not look the same to this
# script, so it asks the installed SZ3 itself before believing anything GROMACS says about it.
#
# Assert which library was wired in, never just that the configure exited 0.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
GMXSRC=$(cd "$1" && pwd)
WORK=$2
SZ3_PREFIX=$3
HDF5_PREFIX=$(cd "$4" && pwd)
EXPECT=$5
shift 5
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)

# Which library GROMACS should end up with. An SZ3 that declines and an SZ3 that was never
# installed both leave it its own copy, which is why the two are named apart here.
case $EXPECT in
found)    WIRED=external ;;
declined) WIRED=internal ;;
absent)   WIRED=internal ;;
*)        echo "unknown expectation $EXPECT: want found, absent or declined"; exit 1 ;;
esac
if [ "$EXPECT" = absent ] && [ "$SZ3_PREFIX" != "-" ]; then
    echo "expecting no SZ3 to be installed, but $SZ3_PREFIX was given as one to find"
    exit 1
fi
if [ "$EXPECT" != absent ] && [ "$SZ3_PREFIX" = "-" ]; then
    echo "expecting an installed SZ3 to be $EXPECT, with no prefix for it to be installed in"
    exit 1
fi

PREFIX_PATH=$HDF5_PREFIX
SZ3_VERSION=
if [ "$SZ3_PREFIX" != "-" ]; then
    SZ3_PREFIX=$(cd "$SZ3_PREFIX" && pwd)
    PREFIX_PATH="$HDF5_PREFIX;$SZ3_PREFIX"
    SZ3_VERSION=$(sed -n 's/^set(PACKAGE_VERSION "\(.*\)").*/\1/p' \
                  "$SZ3_PREFIX"/lib*/cmake/SZ3/SZ3ConfigVersion.cmake 2>/dev/null | head -1)
    # Required, not skipped: empty, the version assertion matches nothing and the check that
    # INTERNAL ignores this prefix drops out silently, which reads as a pass.
    if [ -z "$SZ3_VERSION" ]; then
        echo "no version in $SZ3_PREFIX/lib*/cmake/SZ3/SZ3ConfigVersion.cmake"
        exit 1
    fi
fi

# The cheapest GROMACS that still generates the h5md link line: no SIMD, no GPU, bundled FFT.
# BUILD_TESTING is on only because h5md-test is the target whose link line names the library.
COMMON=(-DCMAKE_BUILD_TYPE=Release
        -DGMX_SIMD=None
        -DGMX_FFT_LIBRARY=fftpack
        -DGMX_GPU=OFF
        -DGMX_OPENMP=OFF
        -DGMX_BUILD_OWN_FFTW=OFF
        -DBUILD_TESTING=ON
        -DGMX_USE_HDF5=ON
        "-DCMAKE_PREFIX_PATH=$PREFIX_PATH")

LINKTXT=src/gromacs/fileio/h5md/tests/CMakeFiles/h5md-test.dir/link.txt

pass=0; fail=0
ok()  { echo "PASS  $1"; pass=$((pass+1)); }
bad() { echo "FAIL  $1"; shift; for l in "$@"; do echo "        $l"; done; fail=$((fail+1)); }
# want <name> <expected substring> <file>
want() { if grep -qF "$2" "$3"; then ok "$1"; else bad "$1" "expected to find: $2" "got:" "$(grep -i sz3 "$3" | head -5)"; fi; }
notwant() { if grep -qF "$2" "$3"; then bad "$1" "did not expect: $2" "$(grep -nF "$2" "$3" | head -3)"; else ok "$1"; fi; }
# Which SZ3 a link token is: the installed one, or one GROMACS built in its own tree. A relative
# path is a build-tree library by construction, so only the install prefix has to be named here.
from_prefix() { while read -r t; do case $t in "$SZ3_PREFIX"/*) echo "$t";; esac; done; }
elsewhere() { while read -r t; do case $t in "$SZ3_PREFIX"/*) ;; *) echo "$t";; esac; done; }
oneline() { if [ -s "$1" ]; then tr '\n' ' ' < "$1" | sed 's/ *$//'; else printf '(none)'; fi; }

# Every link-line check below rests on one question -- does this token link an SZ3 library --
# so make the answer prove itself both ways before anything is asserted with it.
"$HERE"/sz3LinkTokens.sh --self-test || exit 1

EXTRA=("$@")
echo "=== GMX_USE_SZ3 matrix: $GMXSRC ==="
echo "    hdf5    $HDF5_PREFIX"
echo "    sz3     ${SZ3_PREFIX} ${SZ3_VERSION:+(version $SZ3_VERSION)}"
echo "    expect  $EXPECT, so GROMACS should wire in the $WIRED SZ3"
echo "    extra   ${EXTRA[*]:-none}"

# GROMACS asks with find_package(SZ3 ... QUIET), which prints neither the answer nor the reason, so
# every expectation below rests on a premise no GROMACS log can confirm. Put the same question to
# the installed SZ3 directly and stop here if the answer is not the one this run was set up to get:
# a scenario that failed to arrange itself reads, ten lines later, like the defect it looks for.
mkdir -p "$WORK/probe"
cat > "$WORK/probe/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.18)
# C and CXX, as GROMACS enables both: SZ3Config.cmake resolves OpenMP for whichever it is given.
project(sz3probe C CXX)
find_package(SZ3 QUIET)
if(SZ3_FOUND)
    set(answer found)
elseif(SZ3_NOT_FOUND_MESSAGE)
    set(answer declined)
elseif(SZ3_CONSIDERED_CONFIGS)
    # An SZ3Config.cmake that was neither accepted nor declined -- too old a version, say.
    set(answer rejected)
else()
    set(answer absent)
endif()
string(STRIP "${answer} ${SZ3_NOT_FOUND_MESSAGE}" sz3_answer)
message(STATUS "answer: ${sz3_answer}")
if(NOT answer STREQUAL SZ3_EXPECTED_ANSWER)
    message(FATAL_ERROR "the installed SZ3 answers ${answer}, not ${SZ3_EXPECTED_ANSWER}")
endif()
EOF
rm -rf "$WORK/probe/b"
cmake -S "$WORK/probe" -B "$WORK/probe/b" --no-warn-unused-cli \
      "-DCMAKE_PREFIX_PATH=$PREFIX_PATH" "-DSZ3_EXPECTED_ANSWER=$EXPECT" "${EXTRA[@]}" \
      > "$WORK/probe.log" 2>&1
probe_rc=$?
echo "    asked   $(sed -n 's/^-- answer: //p' "$WORK/probe.log" | head -1)"
if [ "$probe_rc" != 0 ]; then
    echo "this run's premise does not hold, so nothing below would mean what it says:"
    tail -6 "$WORK/probe.log" | sed 's/^/        /'
    exit 1
fi

for mode in EXTERNAL AUTO INTERNAL OFF; do
    rm -rf "$WORK/b-$mode"
    start=$(date +%s)
    cmake -S "$GMXSRC" -B "$WORK/b-$mode" "${COMMON[@]}" "-DGMX_USE_SZ3=$mode" "${EXTRA[@]}" \
          > "$WORK/$mode.log" 2>&1
    rc=$?
    echo "--- GMX_USE_SZ3=$mode: cmake exit $rc, $(( $(date +%s) - start ))s"
    LINE=$WORK/b-$mode/$LINKTXT
    : > "$WORK/$mode.link"
    [ -f "$LINE" ] && "$HERE"/sz3LinkTokens.sh "$LINE" > "$WORK/$mode.link"
    from_prefix < "$WORK/$mode.link" > "$WORK/$mode.installed"
    elsewhere   < "$WORK/$mode.link" > "$WORK/$mode.own"

    case "$mode:$WIRED" in
    EXTERNAL:external)
        if [ "$rc" = 0 ]; then ok "EXTERNAL configures"; else bad "EXTERNAL configures" "$(tail -15 "$WORK/$mode.log")"; fi
        want "EXTERNAL reports the installed SZ3 $SZ3_VERSION" \
             "Found external SZ3 library (found version $SZ3_VERSION)" "$WORK/$mode.log"
        if [ -s "$WORK/$mode.installed" ] && [ ! -s "$WORK/$mode.own" ]; then
            ok "EXTERNAL links the installed libhdf5sz3 by path"
            sed 's/^/        /' "$WORK/$mode.installed"
        else
            bad "EXTERNAL links the installed libhdf5sz3 by path" \
                "the unversioned hdf5sz3 target has to come out of the export, not as a bare -l" \
                "SZ3 libraries linked: $(oneline "$WORK/$mode.link")"
        fi
        ;;
    EXTERNAL:internal)
        # GROMACS's own error, not a crash inside SZ3Config.cmake: the caller asked for external.
        if [ "$rc" = 0 ]; then
            bad "EXTERNAL fails when no SZ3 can be used" "cmake exited 0"
        else
            ok "EXTERNAL fails when no SZ3 can be used"
        fi
        want "EXTERNAL fails with GROMACS's own message" \
             "Building with an external SZ3 library was selected" "$WORK/$mode.log"
        ;;
    AUTO:external)
        if [ "$rc" = 0 ]; then ok "AUTO configures"; else bad "AUTO configures" "$(tail -15 "$WORK/$mode.log")"; fi
        want "AUTO takes the installed SZ3" \
             "Found external SZ3 library (found version $SZ3_VERSION)" "$WORK/$mode.log"
        notwant "AUTO does not also build the internal copy" "Using internal SZ3 library" "$WORK/$mode.log"
        ;;
    AUTO:internal)
        if [ "$rc" = 0 ]; then
            ok "AUTO still configures when SZ3 cannot be used"
        else
            bad "AUTO still configures when SZ3 cannot be used" \
                "AUTO's contract is the internal fallback; SZ3 must decline, not end the configure" \
                "$(tail -15 "$WORK/$mode.log")"
        fi
        want "AUTO falls back to the internal SZ3" "Using internal SZ3 library" "$WORK/$mode.log"
        ;;
    INTERNAL:*)
        if [ "$rc" = 0 ]; then ok "INTERNAL configures"; else bad "INTERNAL configures" "$(tail -15 "$WORK/$mode.log")"; fi
        want "INTERNAL uses GROMACS's own copy" "Using internal SZ3 library" "$WORK/$mode.log"
        # Absence on its own would also be satisfied by a build that links no SZ3 whatsoever.
        if [ -s "$WORK/$mode.own" ]; then
            ok "INTERNAL links the SZ3 in GROMACS's own tree"
            sed 's/^/        /' "$WORK/$mode.own"
        else
            bad "INTERNAL links the SZ3 in GROMACS's own tree" \
                "SZ3 libraries linked: $(oneline "$WORK/$mode.link")"
        fi
        # Only when there is an installed SZ3 to ignore: a check that cannot fail is not a check.
        # Any mention of the prefix, not just a library from it: a -L or an rpath into it would
        # still let the link, or the loaded binary, reach the copy INTERNAL is meant to skip.
        if [ -n "$SZ3_VERSION" ]; then
            if [ -f "$LINE" ] && grep -qF "$SZ3_PREFIX" "$LINE"; then
                bad "INTERNAL ignores the installed SZ3" \
                    "$(tr ' ' '\n' < "$LINE" | grep -F "$SZ3_PREFIX" | head -5)"
            else
                ok "INTERNAL ignores the installed SZ3"
            fi
        fi
        ;;
    OFF:*)
        if [ "$rc" = 0 ]; then ok "OFF configures"; else bad "OFF configures" "$(tail -15 "$WORK/$mode.log")"; fi
        if [ -s "$WORK/$mode.link" ]; then
            bad "OFF links no SZ3 at all" "SZ3 libraries linked: $(oneline "$WORK/$mode.link")"
        else
            ok "OFF links no SZ3 at all"
        fi
        ;;
    *)
        # Unreachable while $WIRED is validated above; here so a mode added to the loop is not
        # silently unasserted.
        bad "unhandled GMX_USE_SZ3 value" "$mode:$WIRED"
        ;;
    esac
done

echo "== $pass passed, $fail failed =="
[ "$fail" = 0 ]
