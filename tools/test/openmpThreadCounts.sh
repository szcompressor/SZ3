#!/bin/bash
# An OpenMP stream carries the thread count it was written with, and the reader asks the runtime
# for that many. It may get fewer, and then it must still decode every chunk to the same bytes.
#
#   tools/test/openmpThreadCounts.sh <sz3 binary> <input file> <abs error bound> <dim args...>
#   tools/test/openmpThreadCounts.sh build/tools/sz3/sz3 input.dat 1 -3 128 8 8
#
# Only OMP_THREAD_LIMIT can force the reader to get fewer: it caps the whole program, so a
# num_threads request cannot exceed it. OMP_NUM_THREADS cannot -- num_threads overrides it -- so
# setting that on the reader would leave the count unchanged and test nothing.
#
# Assert the bytes, never just that the exit was 0. No set -e: every check below is explicit, and
# the total at the bottom catches a section that stopped early.
set -u
SZ3=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
INPUT=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
BOUND=$3
shift 3
DIMS=("$@")
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Required, not skipped: a missing binary must never look like a passing check.
if [ ! -x "$SZ3" ]; then echo "no sz3 binary at $SZ3"; exit 1; fi
if [ ! -f "$INPUT" ]; then echo "no input file at $INPUT"; exit 1; fi
if [ ${#DIMS[@]} -eq 0 ]; then echo "no dimension arguments given"; exit 1; fi

cd "$WORK" || exit 1
printf '[GlobalSettings]\nOpenMP = YES\n' > omp.config

pass=0; fail=0
ok()  { echo "PASS  $1"; pass=$((pass+1)); }
bad() { echo "FAIL  $1"; shift; for l in "$@"; do echo "        $l"; done; fail=$((fail+1)); }

echo "=== OpenMP thread counts: $SZ3 ==="
echo "    input  $INPUT"
echo "    dims   ${DIMS[*]}  bound $BOUND"

WRITERS="1 4 8"
LIMITS="1 2 3"

for w in $WRITERS; do
    if OMP_NUM_THREADS=$w "$SZ3" -f -i "$INPUT" "${DIMS[@]}" -M ABS "$BOUND" \
            -c omp.config -z "omp$w.sz3" > "w$w.log" 2>&1 && [ -s "omp$w.sz3" ]; then
        ok "omp-$w-thread-stream-written"
    else
        bad "omp-$w-thread-stream-written" "$(tail -5 "w$w.log")"
    fi
    # The reference this stream's own reader produces, before any cap is put on it.
    if "$SZ3" -f -s "omp$w.sz3" -o "omp$w.ref.dat" "${DIMS[@]}" > "r$w.log" 2>&1; then
        ok "omp-$w-thread-stream-reads-back"
    else
        bad "omp-$w-thread-stream-reads-back" "$(tail -5 "r$w.log")"
    fi
    for l in $LIMITS; do
        if OMP_THREAD_LIMIT=$l "$SZ3" -f -s "omp$w.sz3" -o "omp$w.$l.dat" "${DIMS[@]}" \
                > "r$w.$l.log" 2>&1; then
            ok "omp-$w-thread-stream-reads-under-limit-$l"
        else
            bad "omp-$w-thread-stream-reads-under-limit-$l" "$(tail -5 "r$w.$l.log")"
        fi
        if [ -s "omp$w.ref.dat" ] && cmp -s "omp$w.ref.dat" "omp$w.$l.dat"; then
            ok "omp-$w-thread-stream-identical-under-limit-$l"
        else
            bad "omp-$w-thread-stream-identical-under-limit-$l" \
                "OMP_THREAD_LIMIT=$l changed the data a $w-thread stream decodes to"
        fi
    done
done

# Without this the suite is 24 checks on a library that was built with no OpenMP in it: the config
# key is ignored, every stream is the plain serial one, and all of them round-trip perfectly.
differ=0
for w in 4 8; do
    cmp -s omp1.sz3 "omp$w.sz3" || differ=1
done
if [ "$differ" = 1 ]; then
    ok "thread-count-changes-the-stream"
else
    bad "thread-count-changes-the-stream" \
        "1, 4 and 8 writer threads produced the same bytes, so this binary has no OpenMP path" \
        "$(ls -l omp1.sz3 omp4.sz3 omp8.sz3 2>&1)"
fi

echo
echo "  $pass passed, $fail failed"

# Raise this with the check it comes with. A loop that stopped early otherwise shows only as a
# smaller number at the bottom that nobody compares.
EXPECTED=25
ran=$((pass + fail))
if [ "$ran" -ne "$EXPECTED" ]; then
    echo "  the suite accounted for $ran checks, not $EXPECTED"
    exit 1
fi
exit $((fail > 0))
