#!/bin/bash
# GROMACS writes an SZ3-compressed H5MD trajectory, and everything that reads it back is checked
# against the coordinates GROMACS handed the compressor.
#
#   tools/test/gromacsH5mdRoundtrip.sh <gmx> <plugin dir> <hdf5 bin dir> <traj> <topology> <workdir> [ndec]
#
# ndec is GROMACS's output precision in decimal places; the error bound it asks SZ3 for is
# 1/(2*10^ndec) nm. Assert the numbers, never just that the exit was 0.
set -u
# Everything is resolved before the cd below, so the caller may pass relative paths.
GMX=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
PLUGIN_DIR=$(cd "$2" && pwd)
H5BIN=$(cd "$3" && pwd)
TRAJ=$(cd "$(dirname "$4")" && pwd)/$(basename "$4")
TOPOL=$(cd "$(dirname "$5")" && pwd)/$(basename "$5")
WORK=$6
NDEC=${7:-2}
# PYTHON only has to import numpy; nothing here needs HDF5 or GROMACS bindings.
PY=${PYTHON:-python3}
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
cd "$WORK" || exit 1

DSET=/particles/system/position/value
NOPLUGIN=$WORK/no-such-plugin-dir

pass=0; fail=0
ok()  { echo "PASS  $1"; pass=$((pass+1)); }
bad() { echo "FAIL  $1"; shift; for l in "$@"; do echo "        $l"; done; fail=$((fail+1)); }
# want <name> <expected substring> <file>
want() { if grep -qF "$2" "$3"; then ok "$1"; else bad "$1" "expected to find: $2" "got:" "$(head -20 "$3")"; fi; }

# Required, not skipped: a missing tool must never look like a passing check.
missing=
for tool in h5dump h5repack; do
    command -v "$H5BIN/$tool" > /dev/null 2>&1 || missing="$missing $tool"
done
command -v "$GMX" > /dev/null 2>&1 || missing="$missing $(basename "$GMX")"
"$PY" -c 'import numpy' > /dev/null 2>&1 || missing="$missing python3-numpy"
if [ -n "$missing" ]; then
    echo "required tools not found:$missing"
    exit 1
fi
# Required, not defaulted: empty, the version assertion below looks for "H5Z-SZ3-" and any version
# satisfies it, which is the one thing this pin exists to rule out.
if [ -z "${SZ3_EXPECT_VERSION:-}" ]; then
    echo "set SZ3_EXPECT_VERSION to the version the filter has to report"
    exit 1
fi

BOUND=$("$PY" -c "print(1.0 / (2.0 * 10 ** $NDEC))")
echo "=== $($GMX --version 2>/dev/null | grep -i 'GROMACS version' | head -1 | sed 's/  */ /g') ==="
echo "    trajectory  $TRAJ"
echo "    error bound $BOUND nm (-ndec $NDEC)"

# ------------------------------------------------------------------ write
# Two passes over the same frames with the same selection: the .trr carries the float32
# coordinates GROMACS gave the compressor, so it is the reference the .h5md is judged against.
echo 0 | HDF5_PLUGIN_PATH=$PLUGIN_DIR "$GMX" trjconv -f "$TRAJ" -s "$TOPOL" -o lossy.h5md \
         -ndec "$NDEC" > write_lossy.log 2>&1
lossy_rc=$?
echo 0 | "$GMX" trjconv -f "$TRAJ" -s "$TOPOL" -o ref.trr > write_ref.log 2>&1
ref_rc=$?
if [ "$lossy_rc" = 0 ] && [ -s lossy.h5md ]; then
    ok "gmx trjconv writes an SZ3-compressed .h5md"
else
    bad "gmx trjconv writes an SZ3-compressed .h5md" "exit $lossy_rc" "$(tail -15 write_lossy.log)"
    echo "== $pass passed, $fail failed =="
    exit 1
fi
if [ "$ref_rc" = 0 ] && [ -s ref.trr ]; then
    ok "gmx trjconv writes the uncompressed .trr reference"
else
    bad "gmx trjconv writes the uncompressed .trr reference" "exit $ref_rc" "$(tail -15 write_ref.log)"
    echo "== $pass passed, $fail failed =="
    exit 1
fi

# ------------------------------------------------------------------ the filter is really on it
# No plugin reachable: h5dump can still report the filter, and cannot silently decode instead.
HDF5_PLUGIN_PATH=$NOPLUGIN "$H5BIN/h5dump" -pH -d "$DSET" lossy.h5md > layout.txt 2>&1
want "the position dataset carries filter 32024" "FILTER_ID 32024" layout.txt
# SZ3_EXPECT_VERSION pins which SZ3 did the compressing, for a consumer that vendors one.
want "the filter names itself as SZ3 ${SZ3_EXPECT_VERSION:-}" \
     "H5Z-SZ3-${SZ3_EXPECT_VERSION:-}" layout.txt
sed -n '/STORAGE_LAYOUT/,/FILTERS/p' layout.txt | sed 's/^/    /'

# h5repack decodes every chunk through the plugin and writes the same trajectory with no filter.
HDF5_PLUGIN_PATH=$PLUGIN_DIR "$H5BIN/h5repack" -f NONE lossy.h5md plain.h5md > repack.log 2>&1
repack_rc=$?
# h5repack exits 0 when it cannot decode a dataset, and drops it with a warning on stdout, so the
# exit status on its own says nothing.
repack_ok=yes
[ "$repack_rc" = 0 ] || repack_ok=no
[ -s plain.h5md ] || repack_ok=no
if grep -q "cannot be read" repack.log; then repack_ok=no; fi
if [ "$repack_ok" = yes ]; then
    ok "h5repack decompresses the whole file through the plugin"
else
    bad "h5repack decompresses the whole file through the plugin" "exit $repack_rc" "$(tail -10 repack.log)"
fi
LOSSY_BYTES=$(stat -c %s lossy.h5md)
PLAIN_BYTES=$(stat -c %s plain.h5md 2>/dev/null || echo 0)
if [ "$PLAIN_BYTES" -gt "$LOSSY_BYTES" ]; then
    ok "the compressed file is smaller than the same trajectory uncompressed"
else
    bad "the compressed file is smaller than the same trajectory uncompressed" \
        "compressed $LOSSY_BYTES bytes, uncompressed $PLAIN_BYTES bytes"
fi

# ------------------------------------------------------------------ a reader with no GROMACS in it
cat > rd.c <<'EOF'
/* Reads the position dataset and writes it as raw float32. Links libhdf5 and nothing else, so the
   SZ3 plugin is its only way to the coordinates. */
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: rd <file> <dataset> <out.f32>\n"); return 2; }
    hid_t f = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0) { printf("OPEN FAILED\n"); return 2; }
    hid_t d = H5Dopen2(f, argv[2], H5P_DEFAULT);
    if (d < 0) { printf("DATASET NOT FOUND\n"); return 2; }
    hid_t s = H5Dget_space(d);
    hsize_t dims[3] = {0, 0, 0};
    int nd = H5Sget_simple_extent_dims(s, dims, NULL);
    if (nd != 3) { printf("UNEXPECTED RANK %d\n", nd); return 2; }
    size_t n = (size_t)dims[0] * dims[1] * dims[2];
    float* buf = (float*)malloc(n * sizeof(float));
    if (!buf) { printf("OUT OF MEMORY\n"); return 2; }
    if (H5Dread(d, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0) {
        printf("READ FAILED\n"); return 1;
    }
    FILE* out = fopen(argv[3], "wb");
    if (!out || fwrite(buf, sizeof(float), n, out) != n) { printf("WRITE FAILED\n"); return 2; }
    fclose(out);
    printf("READ OK frames=%llu atoms=%llu dim=%llu\n",
           (unsigned long long)dims[0], (unsigned long long)dims[1], (unsigned long long)dims[2]);
    return 0;
}
EOF
# Built against the HDF5 whose tools are being driven, so the reader and h5dump agree on a version.
H5ROOT=$(cd "$H5BIN/.." && pwd)
if [ -f "$H5ROOT/include/hdf5.h" ]; then
    cc -O1 -o rd rd.c -I"$H5ROOT/include" -L"$H5ROOT/lib" -lhdf5 -Wl,-rpath,"$H5ROOT/lib" \
       > build_rd.log 2>&1
else
    cc -O1 -o rd rd.c $(pkg-config --cflags --libs hdf5 2>/dev/null || echo "-lhdf5") > build_rd.log 2>&1
fi
if [ -x ./rd ]; then
    ok "the plain HDF5 reader builds"
else
    bad "the plain HDF5 reader builds" "$(tail -10 build_rd.log)"
    echo "== $pass passed, $fail failed =="
    exit 1
fi
ldd ./rd > rd.ldd 2>&1
# The dependency names, not the directories they resolve in: a reader that found its HDF5 under
# a directory called sz3 or gromacs is still a reader with neither of them linked into it.
deps=$(awk '{print $1}' rd.ldd | sed 's|.*/||')
banned=$(echo "$deps" | grep -Ei 'gromacs|sz3' | tr '\n' ' ')
if ! echo "$deps" | grep -q '^libhdf5'; then
    # ldd's own error text matches nothing below, which would read as a pass.
    bad "the reader links no GROMACS and no SZ3" "ldd listed no libhdf5:" "$(cat rd.ldd)"
elif [ -n "$banned" ]; then
    bad "the reader links no GROMACS and no SZ3" "$banned"
else
    ok "the reader links no GROMACS and no SZ3"
fi

HDF5_PLUGIN_PATH=$PLUGIN_DIR ./rd lossy.h5md "$DSET" decoded.f32 > read_ok.log 2>&1
want "the reader decodes the coordinates through the plugin" "READ OK" read_ok.log
sed 's/^/    /' read_ok.log

# The repacked copy carries no filter, so this read goes nowhere near SZ3 and is the second,
# independent decode of the same chunks.
./rd plain.h5md "$DSET" repacked.f32 > read_plain.log 2>&1
want "the reader reads the repacked copy" "READ OK" read_plain.log
if cmp -s decoded.f32 repacked.f32; then
    ok "both decodes of the trajectory agree bit for bit"
else
    bad "both decodes of the trajectory agree bit for bit" \
        "the plugin and h5repack disagree about what the file holds"
fi

# The same read with nowhere to find the plugin: HDF5 cannot decode, and must say so.
HDF5_PLUGIN_PATH=$NOPLUGIN ./rd lossy.h5md "$DSET" /dev/null > read_noplugin.log 2>&1
noplugin_rc=$?
if [ "$noplugin_rc" = 0 ]; then
    bad "without the plugin the read fails" "the read succeeded, so the file was not SZ3-compressed"
else
    ok "without the plugin the read fails"
fi

# ------------------------------------------------------------------ the numbers
"$PY" - "$NDEC" "$LOSSY_BYTES" "$PLAIN_BYTES" <<'PY' > compare.txt 2>&1
import struct, sys
import numpy as np

ndec, lossy_bytes, plain_bytes = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])

# GROMACS computes the bound in its own `real`, so the number SZ3 is given is the float32 one.
bound = float(np.float32(1.0 / (2.0 * 10 ** ndec)))

# TRR is XDR big-endian: magic, two string lengths, the block sizes, natoms, step, nre, t, lambda,
# then the box and the coordinates. Nothing here is GROMACS code.
with open("ref.trr", "rb") as fh:
    blob = fh.read()
if struct.unpack_from(">i", blob, 0)[0] != 1993:
    print("TRR MAGIC MISMATCH"); raise SystemExit(1)
strlen = struct.unpack_from(">i", blob, 8)[0]
off = 12 + ((strlen + 3) // 4) * 4
(box_size, _vir, _pres, _top, _sym, x_size, v_size, f_size, natoms) = struct.unpack_from(">9i", blob, off + 8)
off += 4 * 13
_t, _lam = struct.unpack_from(">2f", blob, off)
header = off + 8
frame_bytes = header + box_size + x_size + v_size + f_size
if len(blob) % frame_bytes != 0:
    print("TRR FRAME SIZE MISMATCH", len(blob), frame_bytes); raise SystemExit(1)
nframes = len(blob) // frame_bytes
ref = np.empty((nframes, natoms, 3), dtype=">f4")
for i in range(nframes):
    start = i * frame_bytes + header + box_size
    ref[i] = np.frombuffer(blob, dtype=">f4", count=natoms * 3, offset=start).reshape(natoms, 3)
ref = ref.astype("<f4")

got = np.fromfile("decoded.f32", dtype="<f4")
if got.size % (natoms * 3) != 0:
    print("DECODED SIZE MISMATCH", got.size, natoms); raise SystemExit(1)
got = got.reshape(-1, natoms, 3)
if got.shape[0] < nframes:
    print("DECODED HAS FEWER FRAMES", got.shape[0], nframes); raise SystemExit(1)

# The dataset is extended a chunk at a time, so it can be longer than the trajectory. Only the
# frames GROMACS wrote are compared; what the padding holds is reported, not asserted.
tail = got[nframes:]
err = np.abs(got[:nframes].astype("f8") - ref.astype("f8"))
worst = float(err.max())
rms = float(np.sqrt((err ** 2).mean()))
raw = ref.size * 4

print("frames            %d" % nframes)
print("atoms             %d" % natoms)
print("coordinates       %d (%.1f MiB raw float32)" % (ref.size, raw / 1048576.0))
print("padding frames    %d (all fill value: %s)" % (tail.shape[0], bool((tail == -1).all()) if tail.size else "n/a"))
# The file stores float32, so a reconstruction exactly on the bound can land one ulp beyond it.
ulp = float(np.spacing(np.float32(np.abs(ref).max())))
print("max abs error     %.17g nm" % worst)
print("rms error         %.8g nm" % rms)
print("requested bound   %.17g nm (float32, plus one ulp %.3g)" % (bound, ulp))
print("file compressed   %d bytes" % lossy_bytes)
print("file uncompressed %d bytes" % plain_bytes)
print("ratio vs raw x    %.3f" % (raw / float(lossy_bytes)))
print("ratio vs file     %.3f" % (plain_bytes / float(lossy_bytes)))
if worst <= bound + ulp:
    print("VERDICT within-bound")
else:
    print("VERDICT over-bound")
PY
cat compare.txt | sed 's/^/    /'
want "every coordinate is reconstructed within the requested bound" "VERDICT within-bound" compare.txt

echo "== $pass passed, $fail failed =="
[ "$fail" = 0 ]
