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
# ALGO_BIOMD takes coordinates {frames, atoms, 3}: 5 frames of 1000 rigid waters and a chain of 300 bonded atoms, in nm
python3 - <<'PY'
import array, math, random
random.seed(1)
def unit(v):
    n = math.sqrt(sum(c * c for c in v)); return [c / n for c in v]
def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
def rvec():
    return unit([random.gauss(0, 1) for _ in range(3)])
r, th = 0.09572, math.radians(104.52)
waters = [([random.uniform(0, 3) for _ in range(3)], rvec(), rvec()) for _ in range(1000)]
chain = [rvec() for _ in range(300)]
x = array.array('f')
for f in range(5):
    p = [1.0, 1.0, 1.0]
    for k, d in enumerate(chain):
        p = [p[c] + (0.153 if k % 2 else 0.109) * d[c] for c in range(3)]; x.extend(p)
    for o, u, w in waters:
        v = unit(cross(u, w))
        x.extend(o); x.extend([o[c] + r * u[c] for c in range(3)])
        x.extend([o[c] + r * (math.cos(th) * u[c] + math.sin(th) * v[c]) for c in range(3)])
    waters = [([o[c] + 0.003 * random.gauss(0, 1) for c in range(3)], unit([u[c] + 0.05 * random.gauss(0, 1) for c in range(3)]), w)
              for o, u, w in waters]
    chain = [unit([d[c] + 0.03 * random.gauss(0, 1) for c in range(3)]) for d in chain]
x.tofile(open('coords.f32', 'wb'))
array.array('d', x).tofile(open('coords.f64', 'wb'))
PY
COORD_DIMS=(-3 3 3300 5)

streams=0; fail=0
while IFS='|' read -r name algo settings; do
    for eb in 1e-2 1e-3 1e-4 1e-5 1e-6; do
        printf '[GlobalSettings]\nCmprAlgo = %s\nErrorBoundMode = ABS\nAbsErrorBound = %s\n[AlgoSettings]\n%b' \
            "$algo" "$eb" "$settings" > c.ini
        for t in f d; do
            in=$INPUT; [ $t = d ] && in=in.f64
            dims=("${DIMS[@]}")
            if [ "$algo" = ALGO_BIOMD ]; then in=coords.f$([ $t = d ] && echo 64 || echo 32); dims=("${COORD_DIMS[@]}"); fi
            for w in A B; do
                streams=$((streams + 1))
                rm -f s.sz a.out b.out
                "${!w}" -$t -i "$in" -z s.sz "${dims[@]}" -c c.ini > /dev/null &&
                    "$A" -$t -s s.sz -o a.out "${dims[@]}" > /dev/null &&
                    "$B" -$t -s s.sz -o b.out "${dims[@]}" > /dev/null &&
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
