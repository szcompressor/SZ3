/**
 * @file readLegacyChunks.cpp
 * @brief The filter must still read files earlier releases wrote.
 *
 * The fixture, tools/H5Z-SZ3/test/legacy-v3.3.2.h5, was written by the released SZ3 v3.3.2 and is
 * committed as it came out of that build -- see makeLegacyFixture.c for the exact recipe. It has to
 * come from the old filter: a fixture this tree produced could not detect new code that fails to
 * read old files, which is the whole of what this test is for.
 *
 * Up to v3.3.2 the filter returned a chunk untouched whenever cd_values said fewer than 20
 * elements, on write as much as on read, so `small_chunk` sits in that file raw and headerless
 * while `normal_chunk` is SZ3-compressed. Both have to read.
 */

#include <hdf5.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "H5Z_SZ3.hpp"

namespace {

int passes = 0;
int failures = 0;

void ok(const std::string& name) {
    printf("PASS  %s\n", name.c_str());
    passes++;
}

void bad(const std::string& name, const std::string& detail) {
    printf("FAIL  %s\n        %s\n", name.c_str(), detail.c_str());
    failures++;
}

void check(const std::string& name, bool cond, const std::string& detail = "") { cond ? ok(name) : bad(name, detail); }

float val(long i) { return static_cast<float>(i) * 0.25f; }

// Collects the messages HDF5 stacked up, so a failure can be asserted on what it told the user and
// not merely on the fact that it failed.
herr_t collect(unsigned, const H5E_error2_t* err, void* out) {
    static_cast<std::string*>(out)->append(err->desc).append("\n");
    return 0;
}

std::string error_stack() {
    std::string s;
    H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, collect, &s);
    return s;
}

// Reads one dataset whole. Returns false and fills `why` with the HDF5 error stack on failure.
bool read_dataset(const char* path, const char* name, std::vector<float>& out, std::string& why) {
    hid_t f = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0) {
        why = "cannot open " + std::string(path);
        return false;
    }
    hid_t ds = H5Dopen2(f, name, H5P_DEFAULT);
    if (ds < 0) {
        H5Fclose(f);
        why = "cannot open dataset " + std::string(name);
        return false;
    }
    hid_t sp = H5Dget_space(ds);
    out.resize(static_cast<size_t>(H5Sget_simple_extent_npoints(sp)));
    bool okread = H5Dread(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) >= 0;
    if (!okread) why = error_stack();
    H5Sclose(sp);
    H5Dclose(ds);
    H5Fclose(f);
    return okread;
}

const char* kProbeMessage = "H5Z probe filter declines";

// Returns 0 on a read after pushing kProbeMessage, so the caller can see whether HDF5 kept it.
size_t probe_filter(unsigned flags, size_t, const unsigned int*, size_t nbytes, size_t*, void**) {
    if (flags & H5Z_FLAG_REVERSE) {
        H5Epush2(H5E_DEFAULT, __FILE__, "probe_filter", __LINE__, H5E_ERR_CLS, H5E_PLINE, H5E_CALLBACK, "%s",
                 kProbeMessage);
        return 0;
    }
    return nbytes;
}

// HDF5 1.14 discards whatever a filter pushes onto the error stack; 1.10 and 2.x pass it through.
// Probed with a filter of our own rather than read off the version number, so the assertion on the
// filter's wording is relaxed only where the library really does throw it away.
bool hdf5_keeps_filter_messages() {
    const H5Z_filter_t id = 32099;  // unassigned, and only ever registered inside this process
    const H5Z_class2_t probe[1] = {
        {H5Z_CLASS_T_VERS, id, 1, 1, "H5Z probe", NULL, NULL, static_cast<H5Z_func_t>(probe_filter)}};
    if (H5Zregister(probe) < 0) return false;
    const char* path = "h5zProbe.h5";
    hsize_t dims[1] = {64}, ch[1] = {64};
    float v[64] = {0};
    hid_t f = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t s = H5Screate_simple(1, dims, NULL);
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(p, 1, ch);
    H5Pset_filter(p, id, H5Z_FLAG_MANDATORY, 0, NULL);
    hid_t ds = H5Dcreate2(f, "ds", H5T_IEEE_F32LE, s, H5P_DEFAULT, p, H5P_DEFAULT);
    H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
    H5Dclose(ds);
    H5Pclose(p);
    H5Sclose(s);
    H5Fclose(f);
    std::vector<float> out;
    std::string why;
    bool kept = !read_dataset(path, "ds", out, why) && why.find(kProbeMessage) != std::string::npos;
    H5Zunregister(id);
    remove(path);
    return kept;
}

// A chunk that is neither an SZ3 payload nor the legacy shape. H5Dwrite_chunk puts the bytes in
// with the filter marked as applied, so the read pipeline runs over something SZ3 never wrote.
bool write_foreign_chunk(const char* path) {
    hsize_t dims[1] = {64}, ch[1] = {64}, off[1] = {0};
    hid_t f = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t s = H5Screate_simple(1, dims, NULL);
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(p, 1, ch);
    H5Pset_filter(p, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, 0, NULL);
    hid_t ds = H5Dcreate2(f, "ds", H5T_IEEE_F32LE, s, H5P_DEFAULT, p, H5P_DEFAULT);
    std::vector<unsigned char> junk(96, 0xAB);
    bool okwrite = ds >= 0 && H5Dwrite_chunk(ds, H5P_DEFAULT, 0, off, junk.size(), junk.data()) >= 0;
    if (ds >= 0) H5Dclose(ds);
    H5Pclose(p);
    H5Sclose(s);
    H5Fclose(f);
    return okwrite;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("usage: readLegacyChunks <legacy-v3.3.2.h5>\n");
        return 2;
    }
    const char* fixture = argv[1];
    if (H5Z_SZ3_initialize() < 0) {
        printf("FATAL: could not register the filter\n");
        return 1;
    }
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    // The regression: v3.3.2 stored this chunk raw, so it carries no SZ3 header to recognize.
    std::vector<float> d;
    std::string why;
    if (!read_dataset(fixture, "small_chunk", d, why)) {
        bad("v3.3.2 chunk under 20 elements reads", why);
    } else {
        size_t wrong = 0;
        for (size_t i = 0; i < d.size(); i++) {
            if (d[i] != val(static_cast<long>(i))) wrong++;
        }
        check("v3.3.2 chunk under 20 elements reads", d.size() == 16 && wrong == 0,
              "got " + std::to_string(d.size()) + " elements, " + std::to_string(wrong) + " of them wrong");
    }

    // The control. Without it the check above would also pass on a filter that accepts anything.
    if (!read_dataset(fixture, "normal_chunk", d, why)) {
        bad("v3.3.2 compressed chunk reads", why);
    } else {
        double worst = 0;
        for (size_t i = 0; i < d.size(); i++) {
            worst = std::max(worst, std::fabs(static_cast<double>(d[i]) - val(static_cast<long>(i))));
        }
        check("v3.3.2 compressed chunk reads", d.size() == 4096 && worst <= 1e-3,
              "got " + std::to_string(d.size()) + " elements, worst error " + std::to_string(worst));
    }

    // And the fallback stays narrow: something SZ3 did not write is still refused, and says so.
    const bool messages_survive = hdf5_keeps_filter_messages();
    const std::string name = messages_survive ? "a chunk SZ3 did not write is refused"
                                              : "a chunk SZ3 did not write is refused (this HDF5 drops the "
                                                "filter's own message, so only the refusal is checked)";
    const char* foreign = "foreignChunk.h5";
    if (!write_foreign_chunk(foreign)) {
        bad(name, "could not build the foreign-chunk file");
    } else if (read_dataset(foreign, "ds", d, why)) {
        bad(name, "the read succeeded");
    } else {
        check(name, !messages_survive || why.find("carries no SZ3 header") != std::string::npos,
              "the message did not say the chunk has no SZ3 header:\n" + why);
    }
    remove(foreign);

    H5Z_SZ3_finalize();
    printf("  %d passed, %d failed\n", passes, failures);
    return failures == 0 ? 0 : 1;
}
