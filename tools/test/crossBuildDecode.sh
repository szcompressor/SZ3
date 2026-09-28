#!/bin/bash
# Data compressed by either of two builds (e.g. one with FMA, one without) must decompress to the same values
# in both.
#
#   tools/test/crossBuildDecode.sh <build dir A> <build dir B> <float32 input> <dim args...>
#   tools/test/crossBuildDecode.sh build build-fma tools/sz3/testfloat_8_8_128.dat -3 128 8 8
set -u
A=$(cd "$1" && pwd)/tools/sz3/sz3
B=$(cd "$2" && pwd)/tools/sz3/sz3
INPUT=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
shift 3
DIMS=("$@")
for f in "$A" "$B"; do [ -x "$f" ] || { echo "no sz3 binary at $f"; exit 1; }; done
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK" || exit 1
python3 -c "import array, sys; array.array('d', array.array('f', open(sys.argv[1], 'rb').read())).tofile(open('in.f64', 'wb'))" "$INPUT"

streams=0; fail=0
while IFS='|' read -r name algo settings; do
    for eb in 1e-2 1e-3 1e-4 1e-5 1e-6; do
        printf '[GlobalSettings]\nCmprAlgo = %s\nErrorBoundMode = ABS\nAbsErrorBound = %s\n[AlgoSettings]\n%b' \
            "$algo" "$eb" "$settings" > c.ini
        for t in f d; do
            in=$INPUT; [ $t = d ] && in=in.f64
            for w in A B; do
                streams=$((streams + 1))
                rm -f s.sz a.out b.out
                "${!w}" -$t -i "$in" -z s.sz "${DIMS[@]}" -c c.ini > /dev/null &&
                    "$A" -$t -s s.sz -o a.out "${DIMS[@]}" > /dev/null &&
                    "$B" -$t -s s.sz -o b.out "${DIMS[@]}" > /dev/null &&
                    cmp -s a.out b.out ||
                    { echo "MISMATCH $name eb=$eb type=$t written by $w"; fail=$((fail + 1)); }
            done
        done
    done
done <<'LIST'
lorenzo_reg|ALGO_LORENZO_REG|
lorenzo2_reg|ALGO_LORENZO_REG|Lorenzo2ndOrder = YES\n
regression|ALGO_LORENZO_REG|Lorenzo = NO\n
interp_linear|ALGO_INTERP|InterpolationAlgo = INTERP_ALGO_LINEAR\n
interp_cubic|ALGO_INTERP|InterpolationAlgo = INTERP_ALGO_CUBIC\n
interp_lorenzo|ALGO_INTERP_LORENZO|
biomd|ALGO_BIOMD|
LIST

echo "$fail of $streams streams decode differently"
[ "$streams" -eq 140 ] && [ "$fail" -eq 0 ]
