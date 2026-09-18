#!/bin/bash
# Every way an application or a user reaches the SZ3 HDF5 filter, asserted against an install tree.
#
#   tools/H5Z-SZ3/test/filter_access_modes.sh <install-prefix> [<hdf5 bin dir>]
#
# CMAKE_PREFIX_PATH is passed through for anything else the probes need. Each check asserts which
# path was taken, not only that the command exited zero.
set -u
PREFIX=$(cd "$1" && pwd)
H5BIN=${2:-}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK" || exit 1

h5() { if [ -n "$H5BIN" ]; then echo "$H5BIN/$1"; else echo "$1"; fi; }
PLUGIN_DIR=$PREFIX/lib/plugin
NOPLUGIN=$WORK/no-such-plugin-dir
LIBDIR=$PREFIX/lib
[ -d "$LIBDIR" ] || LIBDIR=$PREFIX/lib64

pass=0; fail=0; skip=0
ok()   { echo "PASS  $1"; pass=$((pass+1)); }
bad()  { echo "FAIL  $1"; shift; for l in "$@"; do echo "        $l"; done; fail=$((fail+1)); }
# Never counted as a pass. A suite that reports more checks than it ran is the failure this
# whole file exists to catch.
skipped() { echo "SKIP  $1"; skip=$((skip+1)); }
# want <name> <expected substring> <file> -- the diagnostic a user sees, not just an exit code
want() { if grep -qF "$2" "$3"; then ok "$1"; else bad "$1" "expected to find: $2" "got:" "$(head -5 "$3")"; fi; }
notwant() { if grep -qF "$2" "$3"; then bad "$1" "did not expect: $2"; else ok "$1"; fi; }

# Required, not optional: they cover how a user reaches this filter without writing code. Refused
# up front so a missing tool can never be mistaken for a passing check.
missing=
for tool in h5repack h5dump h5ls; do
    command -v "$(h5 $tool)" > /dev/null 2>&1 || missing="$missing $tool"
done
if [ -n "$missing" ]; then
    echo "HDF5 command-line tools not found:$missing"
    echo "Install them (Debian/Ubuntu: hdf5-tools) or pass their directory as the second argument."
    exit 1
fi

echo "=== prefix $PREFIX ==="
"$(h5 h5dump)" --version 2>&1 | head -1 | sed 's/^/    /'

# ---------------------------------------------------------------- fixtures
cat > gen.c <<'EOF'
/* Writes an uncompressed float dataset; no SZ3 involved. */
#include <hdf5.h>
#include <math.h>
int main(int argc, char **argv) {
    hsize_t dims[2] = {64, 64}, ch[2] = {16, 64};
    float d[64 * 64];
    for (int i = 0; i < 64 * 64; i++) d[i] = (float)sin(0.01 * i) * 100.0f;
    hid_t f = H5Fcreate(argv[1], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t s = H5Screate_simple(2, dims, NULL);
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(p, 2, ch);
    hid_t d2 = H5Dcreate2(f, "ds", H5T_NATIVE_FLOAT, s, H5P_DEFAULT, p, H5P_DEFAULT);
    herr_t r = H5Dwrite(d2, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, d);
    H5Dclose(d2); H5Pclose(p); H5Sclose(s); H5Fclose(f);
    return r < 0;
}
EOF
# A reader that links no SZ3 at all: the plugin is its only way in.
cat > read.c <<'EOF'
#include <hdf5.h>
#include <math.h>
#include <stdio.h>
int main(int argc, char **argv) {
    float d[64 * 64];
    hid_t f = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0) { printf("OPEN FAILED\n"); return 2; }
    hid_t ds = H5Dopen2(f, "ds", H5P_DEFAULT);
    if (H5Dread(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, d) < 0) {
        printf("READ FAILED\n"); return 1;
    }
    double m = 0;
    for (int i = 0; i < 64 * 64; i++) {
        double e = fabs((double)d[i] - sin(0.01 * i) * 100.0);
        if (e > m) m = e;
    }
    printf("READ OK maxerr %.6f\n", m);
    return m <= 1e-3 * 1.000001 ? 0 : 1;
}
EOF
# Three link/registration shapes an application can take. argv[1] picks one.
cat > app.c <<'EOF'
#include <H5Z_SZ3.hpp>
#include <math.h>
#include <stdio.h>
int main(int argc, char **argv) {
    int use_init = argv[1][0] == 'i' || argv[1][0] == 'b';
    if (use_init) printf("INIT RETURNED %d\n", (int)H5Z_SZ3_initialize());
    printf("AVAIL %d\n", (int)H5Zfilter_avail(H5Z_FILTER_SZ3));
    float d[64 * 64];
    for (int i = 0; i < 64 * 64; i++) d[i] = (float)sin(0.01 * i) * 100.0f;
    hsize_t dims[2] = {64, 64}, ch[2] = {16, 64};
    hid_t f = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t s = H5Screate_simple(2, dims, NULL);
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(p, 2, ch);
    SZ3::Config conf(16, 64);
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-3;
    set_SZ3_conf_to_H5(p, conf);
    hid_t d2 = H5Dcreate2(f, "ds", H5T_NATIVE_FLOAT, s, H5P_DEFAULT, p, H5P_DEFAULT);
    if (d2 < 0 || H5Dwrite(d2, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, d) < 0) {
        printf("WRITE FAILED\n"); return 1;
    }
    H5Dclose(d2); H5Pclose(p); H5Sclose(s); H5Fclose(f);
    printf("WRITE OK\n");
    if (use_init) printf("FINI RETURNED %d\n", (int)H5Z_SZ3_finalize());
    return 0;
}
EOF

printf '#include <hdf5.h>\nint main(void){return H5open() < 0;}\n' > noref.c
cat > CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(modes C CXX)
find_package(SZ3 REQUIRED)
find_package(HDF5 COMPONENTS C REQUIRED)
find_library(MATH_LIB m)
if (NOT MATH_LIB)
    set(MATH_LIB "")
endif ()
add_executable(gen gen.c)
add_executable(read read.c)
target_link_libraries(gen PRIVATE HDF5::HDF5 ${MATH_LIB})
target_link_libraries(read PRIVATE HDF5::HDF5 ${MATH_LIB})
add_executable(app app.c)
set_source_files_properties(app.c PROPERTIES LANGUAGE CXX)
target_link_libraries(app PRIVATE SZ3::SZ3 SZ3::hdf5sz3)
# The shape GROMACS had: linked, referencing nothing. It is the control for the DT_NEEDED
# check below -- without it that check passes vacuously on a toolchain that never drops.
add_executable(noref noref.c)
set_source_files_properties(noref.c PROPERTIES LANGUAGE CXX)
target_link_libraries(noref PRIVATE SZ3::SZ3 SZ3::hdf5sz3)
EOF
# When the tools were given explicitly, that HDF5 is the one to build against as well. Anything
# else the consumer needs goes in CMAKE_PREFIX_PATH -- on Homebrew that is where zstd.h lives,
# since the export deliberately does not put a whole prefix on a consumer's include path.
CMPFX=$PREFIX
[ -n "$H5BIN" ] && CMPFX="$PREFIX;$(cd "$H5BIN/.." && pwd)"
[ -n "${CMAKE_PREFIX_PATH:-}" ] && CMPFX="$CMPFX;$CMAKE_PREFIX_PATH"
cmake -S . -B b -DCMAKE_PREFIX_PATH="$CMPFX" -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 \
  && cmake --build b -j 4 >> cmake.log 2>&1 \
  || { echo "FATAL: could not build the probes"; tail -20 cmake.log; exit 1; }
./b/gen plain.h5 || { echo "FATAL: could not write the fixture"; exit 1; }

# The cd_values a user has to type. Mirrors SZ3::Config::save(); see cdvalueHelper.py.
CD_ABS="UD=32024,0,8,32,0,16777216,4054449152,1348619730,41023,256,0"

# ---------------------------------------------------------------- 1. h5repack
export HDF5_PLUGIN_PATH=$PLUGIN_DIR
"$(h5 h5repack)" -f "$CD_ABS" plain.h5 rp.h5 > rp.log 2>&1
"$(h5 h5dump)" -pH rp.h5 > rp.head 2>&1
want "h5repack-applies-sz3"        "FILTER_ID 32024" rp.head
want "h5repack-records-version"    "H5Z-SZ3-" rp.head
./b/read rp.h5 > rp.read 2>&1
want "h5repack-output-reads-back"  "READ OK" rp.read
# A COMMENT line only appears when the filter was registered as the dataset was created, so it
# is also the tell that the filter really ran rather than being skipped as optional.
"$(h5 h5repack)" -f NONE rp.h5 rp_none.h5 > strip.log 2>&1
"$(h5 h5dump)" -pH rp_none.h5 > strip.head 2>&1
notwant "h5repack-none-strips-sz3" "32024" strip.head
"$(h5 h5ls)" rp_none.h5 > strip.ls 2>&1
want "h5repack-none-keeps-dataset" "ds" strip.ls

# The trap: with the filter unreachable h5repack exits 0 and quietly writes an UNFILTERED copy.
export HDF5_PLUGIN_PATH=$NOPLUGIN
"$(h5 h5repack)" -f "$CD_ABS" plain.h5 rp_noplug.h5 > rpn.log 2>&1
"$(h5 h5dump)" -pH rp_noplug.h5 > rpn.head 2>&1
notwant "h5repack-noplugin-drops-filter-silently" "32024" rpn.head
# and stripping a filtered file loses the dataset outright, also at exit 0
"$(h5 h5repack)" -f NONE rp.h5 rp_lost.h5 > lost.log 2>&1
want "h5repack-noplugin-warns-on-read" "filter is not available" lost.log
"$(h5 h5ls)" rp_lost.h5 > lost.ls 2>&1
notwant "h5repack-noplugin-loses-dataset" "ds" lost.ls

# A cd_values array this filter did not write. The element count cannot be what identifies it:
# it is dataset-dependent and has been 9, 11, 13, 14 and 15 across released versions, so
# H5Z_sz3_set_local has to refuse on content.
# h5repack turns that refusal into the same quiet unfiltered copy as a missing plugin, so what
# is asserted is that no dataset was written claiming a filter that never parsed its parameters.
export HDF5_PLUGIN_PATH=$PLUGIN_DIR
"$(h5 h5repack)" -f "UD=32024,0,9,3,0,3,3341,20,0,1062232653,3539053052,0" plain.h5 rp_old.h5 > old.log 2>&1
"$(h5 h5dump)" -pH rp_old.h5 > old.head 2>&1
notwant "foreign-cdvalues-refused" "32024" old.head
./b/read rp_old.h5 > old.read 2>&1
want "foreign-cdvalues-leaves-data-intact" "READ OK" old.read

# ---------------------------------------------------------------- 2. h5dump / h5ls
export HDF5_PLUGIN_PATH=$NOPLUGIN
"$(h5 h5dump)" -pH rp.h5 > d_meta.head 2>&1
want "h5dump-header-needs-no-plugin"   "FILTER_ID 32024" d_meta.head
want "h5dump-header-shows-version"     "H5Z-SZ3-" d_meta.head
"$(h5 h5ls)" -v rp.h5 > d_meta.ls 2>&1
want "h5ls-verbose-shows-version"      "H5Z-SZ3-" d_meta.ls
if "$(h5 h5dump)" -d /ds rp.h5 > d_data.out 2> d_data.err; then
    bad "h5dump-data-noplugin-fails" "h5dump exited 0 with the filter unreachable"
else
    ok "h5dump-data-noplugin-fails"
fi
want "h5dump-data-noplugin-message"    "unable to print data" d_data.err
"$(h5 h5dump)" --enable-error-stack -d /ds rp.h5 > d_stack.out 2> d_stack.err
want "h5dump-names-the-filter"         "is not registered" d_stack.err
want "h5dump-names-our-version"        "H5Z-SZ3-" d_stack.err
export HDF5_PLUGIN_PATH=$PLUGIN_DIR
"$(h5 h5dump)" -d /ds rp.h5 > d_ok.out 2> d_ok.err
want "h5dump-data-with-plugin"         "DATA {" d_ok.out
notwant "h5dump-data-with-plugin-clean" "unable to print data" d_ok.err

# ---------------------------------------------------------------- 3. application shapes
# (a) links and calls H5Z_SZ3_initialize(), nothing on the plugin path
export HDF5_PLUGIN_PATH=$NOPLUGIN
./b/app init a_init.h5 > a_init.log 2>&1
want "app-init-registers"        "INIT RETURNED 1" a_init.log
want "app-init-writes"           "WRITE OK" a_init.log
want "app-init-finalizes"        "FINI RETURNED 1" a_init.log
./b/read a_init.h5 > /dev/null 2>&1 && bad "app-init-file-needs-filter" "read without the filter succeeded" \
  || ok "app-init-file-needs-filter"
# (b) links but never calls it, reaching the filter only through HDF5_PLUGIN_PATH
export HDF5_PLUGIN_PATH=$PLUGIN_DIR
./b/app plugin a_plug.h5 > a_plug.log 2>&1
want "app-plugin-only-writes"    "WRITE OK" a_plug.log
want "app-plugin-only-avail"     "AVAIL 1" a_plug.log
# (c) does both: the guard must defer to what the plugin already registered, and say so
./b/app both a_both.h5 > a_both.log 2>&1
want "app-both-defers"           "INIT RETURNED 0" a_both.log
want "app-both-writes"           "WRITE OK" a_both.log
# finalize must not unregister a filter it did not register
want "app-both-finalize-declines" "FINI RETURNED 1" a_both.log
# (d) the plugin-only reader, the mode every third-party tool uses
./b/read a_both.h5 > r_plug.log 2>&1
want "plugin-only-reader-works"  "READ OK" r_plug.log
export HDF5_PLUGIN_PATH=$NOPLUGIN
./b/read a_both.h5 > r_noplug.log 2>&1
want "plugin-only-reader-fails-without" "READ FAILED" r_noplug.log
# Whether an unreferenced library is dropped is up to the backend compiler -- Debian's gcc has
# --as-needed in its link spec, the conda-forge driver does not -- so without the unreferencing
# control below this check passes by itself on any toolchain that keeps everything.
if ! command -v readelf > /dev/null 2>&1; then
    skipped "app-keeps-dt-needed (no readelf)"
elif [ ! -f ./b/noref ] || [ ! -f ./b/app ]; then
    bad "app-keeps-dt-needed" "one of the two probes is missing, so there is nothing to compare"
else
    n_app=$(readelf -d ./b/app   | grep -c 'NEEDED.*hdf5sz3')
    n_ref=$(readelf -d ./b/noref | grep -c 'NEEDED.*hdf5sz3')
    if [ "$n_ref" -gt 0 ]; then
        skipped "app-keeps-dt-needed (this toolchain keeps unreferenced libraries too)"
    elif [ "$n_app" -gt 0 ]; then
        ok "app-keeps-dt-needed"
    else
        bad "app-keeps-dt-needed" "calling H5Z_SZ3_initialize() did not keep libhdf5sz3 as a DT_NEEDED"
    fi
fi

echo
echo "  $pass passed, $fail failed, $skip skipped"
exit $((fail > 0))
