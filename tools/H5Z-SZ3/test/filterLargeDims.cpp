// A chunk with a side past 4096 stores a Config longer than a default-constructed one estimates,
// which is the case where reading cd_values back used to run off the end of its own buffer.
#include <cmath>
#include <cstdio>
#include <vector>

#include "H5Z_SZ3.hpp"

static const hsize_t DIMS[3] = {8192, 4, 4};
static const size_t NUM = 8192 * 4 * 4;
static const double ABS_BOUND = 1e-3;
// The quantizer lands values a hair under the bound -- 0.99999x measured on clang and on gcc -- so
// leave room for the float rounding another compiler's reconstruction could add on top.
static const double BOUND_SLACK = 1.001;

static int failures = 0;

static void check(const char* what, bool ok) {
    printf("  %-56s %s\n", what, ok ? "ok" : "FAILED");
    if (!ok) failures++;
}

static std::vector<unsigned char> serialize(const SZ3::Config& conf) {
    std::vector<unsigned char> bytes(conf.size_est());
    auto pos = bytes.data();
    conf.save(pos);
    return bytes;
}

static void report_config(const char* what, const SZ3::Config& expect, const SZ3::Config& got) {
    auto a = serialize(expect), b = serialize(got);
    check(what, a == b);
    if (a == b) return;
    printf("    stored   :");
    for (auto v : a) printf(" %02x", v);
    printf("\n    read back:");
    for (auto v : b) printf(" %02x", v);
    printf("\n    read back dims =");
    for (auto d : got.dims) printf(" %zu", d);
    printf(", dataType %d, quantbinCnt %d, blockSize %d, predDim %d\n", static_cast<int>(got.dataType), got.quantbinCnt,
           got.blockSize, static_cast<int>(got.predDim));
}

static SZ3::Config make_config() {
    SZ3::Config conf(DIMS[0], DIMS[1], DIMS[2]);
    // Two bounds rather than one: the second double pushes the stored config further past the
    // capacity a default-constructed Config asks for, so the tail fields are squarely outside it.
    conf.errorBoundMode = SZ3::EB_ABS_AND_REL;
    conf.absErrorBound = ABS_BOUND;
    conf.relErrorBound = 1.0;
    conf.dataType = SZ_FLOAT;
    return conf;
}

static std::vector<float> make_data() {
    std::vector<float> data(NUM);
    for (size_t i = 0; i < NUM; i++) {
        data[i] = static_cast<float>(std::sin(static_cast<double>(i) * 0.001));
    }
    return data;
}

static void round_trip(const char* path, const char* name, hid_t dcpl, const std::vector<float>& in) {
    hid_t file = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t space = H5Screate_simple(3, DIMS, NULL);
    hid_t dset = H5Dcreate2(file, name, H5T_IEEE_F32LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    check("dataset created", dset >= 0);
    if (dset < 0) {
        H5Sclose(space);
        H5Fclose(file);
        return;
    }
    check("dataset written", H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, in.data()) >= 0);

    hid_t stored_dcpl = H5Dget_create_plist(dset);
    SZ3::Config stored;
    check("config readable from the dataset", get_SZ3_conf_from_H5(stored_dcpl, stored) == 1);
    check("stored config kept all three dimensions", stored.dims == std::vector<size_t>(DIMS, DIMS + 3));
    check("stored config kept its float datatype", stored.dataType == SZ_FLOAT);
    check("stored config kept its quantiser bin count", stored.quantbinCnt == SZ3::Config().quantbinCnt);
    check("stored config kept the block size for 3D", stored.blockSize == 6);
    check("stored config kept the prediction dimension", stored.predDim == 3);
    H5Pclose(stored_dcpl);

    H5Dclose(dset);
    H5Sclose(space);
    H5Fclose(file);

    std::vector<float> out(NUM, 0);
    file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen2(file, name, H5P_DEFAULT);
    check("dataset reopened", dset >= 0);
    check("dataset read back",
          dset >= 0 && H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) >= 0);
    double worst = 0;
    for (size_t i = 0; i < NUM; i++) {
        worst = std::max(worst, std::fabs(static_cast<double>(out[i]) - in[i]));
    }
    printf("    worst absolute error %.6e against a bound of %.6e\n", worst, ABS_BOUND);
    check("every value came back inside the bound", worst <= ABS_BOUND * BOUND_SLACK);
    if (dset >= 0) H5Dclose(dset);
    H5Fclose(file);
}

int main(int argc, char** argv) {
    const char* path = argc > 1 ? argv[1] : "filterLargeDims.h5";
    if (H5Z_SZ3_initialize() < 0) {
        printf("H5Z_SZ3_initialize failed\n");
        return 2;
    }

    SZ3::Config conf = make_config();
    printf("config for %llux%llux%llu serialises to %zu bytes; a default Config estimates %zu\n",
           static_cast<unsigned long long>(DIMS[0]), static_cast<unsigned long long>(DIMS[1]),
           static_cast<unsigned long long>(DIMS[2]), conf.size_est(), SZ3::Config().size_est());

    printf("\nget_SZ3_conf_from_H5 on a property list carrying it\n");
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 3, DIMS);
    check("config stored on the property list", set_SZ3_conf_to_H5(dcpl, conf) == 1);
    SZ3::Config got;
    check("config read back", get_SZ3_conf_from_H5(dcpl, got) == 1);
    report_config("config survives the round trip unchanged", conf, got);

    std::vector<float> in = make_data();
    hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    check("file created", file >= 0);
    if (file >= 0) H5Fclose(file);

    printf("\nthrough the filter, from a property list the caller filled\n");
    round_trip(path, "fresh", dcpl, in);

    // What h5repack does: create a dataset from a property list that already carries a stored
    // config, so set_local reads back the longer one rather than the caller's.
    printf("\nthrough the filter, from a property list the filter itself filled\n");
    file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t first = H5Dopen2(file, "fresh", H5P_DEFAULT);
    hid_t reused = H5Dget_create_plist(first);
    H5Dclose(first);
    H5Fclose(file);
    round_trip(path, "repacked", reused, in);
    H5Pclose(reused);

    H5Pclose(dcpl);
    H5Z_SZ3_finalize();
    printf("\n%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
