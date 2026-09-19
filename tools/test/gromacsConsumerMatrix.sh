#!/bin/bash
# Every value of GROMACS's GMX_USE_SZ3, asserted against an installed SZ3.
#
#   tools/test/gromacsConsumerMatrix.sh <gromacs source> <workdir> <sz3 prefix|-> <hdf5 prefix> \
#                                       <external|internal> [extra cmake args...]
#
# The fifth argument says what an installed SZ3 is expected to do here: be found (external), or
# decline so GROMACS falls back to its own copy (internal). Declining and being absent look the
# same to GROMACS, which is the point -- an SZ3 that ends the configure instead of declining
# breaks GMX_USE_SZ3=AUTO, whose contract is the fallback.
#
# Assert which library was wired in, never just that the configure exited 0.
set -u
GMXSRC=$(cd "$1" && pwd)
WORK=$2
SZ3_PREFIX=$3
HDF5_PREFIX=$(cd "$4" && pwd)
EXPECT=$5
shift 5
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)

if [ "$SZ3_PREFIX" = "-" ] && [ "$EXPECT" = "external" ]; then
    echo "expecting SZ3 to be found, with no SZ3 prefix to find it in"
    exit 1
fi

PREFIX_PATH=$HDF5_PREFIX
SZ3_VERSION=
if [ "$SZ3_PREFIX" != "-" ]; then
    SZ3_PREFIX=$(cd "$SZ3_PREFIX" && pwd)
    PREFIX_PATH="$HDF5_PREFIX;$SZ3_PREFIX"
    SZ3_VERSION=$(sed -n 's/^set(PACKAGE_VERSION "\(.*\)").*/\1/p' \
                  "$SZ3_PREFIX"/lib*/cmake/SZ3/SZ3ConfigVersion.cmake 2>/dev/null | head -1)
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

EXTRA=("$@")
echo "=== GMX_USE_SZ3 matrix: $GMXSRC ==="
echo "    hdf5    $HDF5_PREFIX"
echo "    sz3     ${SZ3_PREFIX} ${SZ3_VERSION:+(version $SZ3_VERSION)}"
echo "    expect  $EXPECT"
echo "    extra   ${EXTRA[*]:-none}"

for mode in EXTERNAL AUTO INTERNAL OFF; do
    rm -rf "$WORK/b-$mode"
    start=$(date +%s)
    cmake -S "$GMXSRC" -B "$WORK/b-$mode" "${COMMON[@]}" "-DGMX_USE_SZ3=$mode" "${EXTRA[@]}" \
          > "$WORK/$mode.log" 2>&1
    rc=$?
    echo "--- GMX_USE_SZ3=$mode: cmake exit $rc, $(( $(date +%s) - start ))s"
    LINE=$WORK/b-$mode/$LINKTXT
    : > "$WORK/$mode.link"
    if [ -f "$LINE" ]; then
        tr ' ' '\n' < "$LINE" | grep -i 'sz3' > "$WORK/$mode.link"
    fi

    case "$mode:$EXPECT" in
    EXTERNAL:external)
        if [ "$rc" = 0 ]; then ok "EXTERNAL configures"; else bad "EXTERNAL configures" "$(tail -15 "$WORK/$mode.log")"; fi
        want "EXTERNAL reports the installed SZ3 $SZ3_VERSION" \
             "Found external SZ3 library (found version $SZ3_VERSION)" "$WORK/$mode.log"
        if grep -qF "$SZ3_PREFIX" "$WORK/$mode.link"; then
            ok "EXTERNAL links the installed libhdf5sz3 by path"
            sed 's/^/        /' "$WORK/$mode.link"
        else
            bad "EXTERNAL links the installed libhdf5sz3 by path" \
                "the unversioned hdf5sz3 target has to come out of the export, not as a bare -l" \
                "link line: $(cat "$WORK/$mode.link")"
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
        # Only when there is an installed SZ3 to ignore: a check that cannot fail is not a check.
        if [ -n "$SZ3_VERSION" ]; then
            if grep -qF "$SZ3_PREFIX" "$WORK/$mode.link"; then
                bad "INTERNAL ignores the installed SZ3" "link line: $(cat "$WORK/$mode.link")"
            else
                ok "INTERNAL ignores the installed SZ3"
            fi
        fi
        ;;
    OFF:*)
        if [ "$rc" = 0 ]; then ok "OFF configures"; else bad "OFF configures" "$(tail -15 "$WORK/$mode.log")"; fi
        if [ -s "$WORK/$mode.link" ]; then
            bad "OFF links no SZ3 at all" "link line: $(cat "$WORK/$mode.link")"
        else
            ok "OFF links no SZ3 at all"
        fi
        ;;
    *)
        bad "unknown expectation" "$mode:$EXPECT"
        ;;
    esac
done

echo "== $pass passed, $fail failed =="
[ "$fail" = 0 ]
