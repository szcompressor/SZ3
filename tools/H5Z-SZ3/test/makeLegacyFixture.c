/*
 * How legacy-v3.3.2.h5 was made. Not built by CMake, and not meant to be re-run in place: the
 * fixture has to come out of the released filter, not out of this tree.
 *
 *   git -C /tmp worktree add sz3-v3.3.2 v3.3.2   # or: git clone -b v3.3.2 <SZ3>
 *   cmake -S /tmp/sz3-v3.3.2 -B /tmp/sz3-v3.3.2/b -DCMAKE_BUILD_TYPE=Release -DBUILD_H5Z_FILTER=ON
 *   cmake --build /tmp/sz3-v3.3.2/b -j
 *   cc makeLegacyFixture.c -o /tmp/mk -lhdf5
 *   HDF5_PLUGIN_PATH=/tmp/sz3-v3.3.2/b/tools/H5Z-SZ3 /tmp/mk legacy-v3.3.2.h5
 *
 * Links no SZ3 itself, so the only SZ3 involved is the one HDF5_PLUGIN_PATH points at.
 */
#include <hdf5.h>
#include <stdio.h>
#define H5Z_FILTER_SZ3 32024
/* Exactly representable in float on every platform, so the reader needs no libm to know the values. */
static float val(long i) { return (float)i * 0.25f; }

static int mk(hid_t f, const char *name, hsize_t d0, hsize_t d1, hsize_t c0, hsize_t c1) {
    hsize_t dims[2] = {d0, d1}, ch[2] = {c0, c1};
    static float d[64 * 64];
    for (long i = 0; i < (long)(d0 * d1); i++) d[i] = val(i);
    hid_t s = H5Screate_simple(2, dims, NULL);
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(p, 2, ch);
    /* cd_nelmts 0: set_local fills cd_values from the chunk shape and the built-in defaults. */
    if (H5Pset_filter(p, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, 0, NULL) < 0) return 1;
    hid_t ds = H5Dcreate2(f, name, H5T_IEEE_F32LE, s, H5P_DEFAULT, p, H5P_DEFAULT);
    if (ds < 0) return 1;
    if (H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, d) < 0) return 1;
    H5Dclose(ds);
    H5Pclose(p);
    H5Sclose(s);
    return 0;
}

int main(int argc, char **argv) {
    (void)argc;
    /* Oldest on-disk format, so the fixture stays readable by every HDF5 SZ3 supports. */
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_libver_bounds(fapl, H5F_LIBVER_EARLIEST, H5F_LIBVER_V18);
    hid_t f = H5Fcreate(argv[1], H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    if (f < 0) {
        printf("CREATE FAILED\n");
        return 1;
    }
    /* 16 elements per chunk: under the 20 the released filter skipped, so stored raw. */
    if (mk(f, "small_chunk", 4, 4, 4, 4)) {
        printf("SMALL FAILED\n");
        return 1;
    }
    /* 4096 per chunk: over it, so genuinely SZ3-compressed. v3.3.2 sized its compression buffer at
       twice the raw chunk, which the bound SZ_compress checks only clears from about here up. */
    if (mk(f, "normal_chunk", 64, 64, 64, 64)) {
        printf("NORMAL FAILED\n");
        return 1;
    }
    H5Fclose(f);
    H5Pclose(fapl);
    printf("FIXTURE OK\n");
    return 0;
}
