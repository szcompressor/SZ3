// The compressor predicts from its own reconstructed values, so the decompressor must reconstruct exactly the
// same values. Built with FMA (-march=x86-64-v3), this fails if only one of the two is protected by nofma().
// ALGO_BIOMD leaves the input as it was and predicts on its integer lattice: the values it decodes, compressed again,
// must decode to themselves.

#include <cmath>
#include <cstring>
#include <vector>

#include "SZ3/api/sz.hpp"
#include "gtest/gtest.h"

namespace {

template <class T>
size_t differing_values(SZ3::ALGO algo, SZ3::INTERP_ALGO interp, double eb) {
    SZ3::Config conf = algo == SZ3::ALGO_BIOMD ? SZ3::Config(32, 341, 3) : SZ3::Config(32, 32, 32);
    conf.cmprAlgo = algo;
    conf.interpAlgo = interp;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
    std::vector<T> data(conf.num), out(conf.num);
    for (size_t i = 0; i < conf.num; i++) {
        double x = i / 1024, y = i / 32 % 32, z = i % 32;
        data[i] = static_cast<T>(std::sin(0.3 * x) * std::cos(0.2 * y) + 0.5 * std::sin(0.1 * z + 0.05 * x));
    }

    std::vector<SZ3::uchar> cmp(SZ3::SZ_compress_size_bound<T>(conf));
    if (algo == SZ3::ALGO_INTERP) {
        size_t size = SZ3::SZ_compress_Interp<T, 3>(conf, data.data(), cmp.data(), cmp.size());
        SZ3::SZ_decompress_Interp<T, 3>(conf, cmp.data(), size, out.data());
    } else {
        size_t size = SZ3::SZ_compress_bioMD<T, 3>(conf, data.data(), cmp.data(), cmp.size());
        SZ3::SZ_decompress_bioMD<T, 3>(conf, cmp.data(), size, data.data());
        size = SZ3::SZ_compress_bioMD<T, 3>(conf, data.data(), cmp.data(), cmp.size());
        SZ3::SZ_decompress_bioMD<T, 3>(conf, cmp.data(), size, out.data());
    }
    size_t differ = 0;
    for (size_t i = 0; i < conf.num; i++) {
        differ += std::memcmp(&data[i], &out[i], sizeof(T)) != 0;
    }
    return differ;
}

template <class T>
void expect_identical(SZ3::ALGO algo, SZ3::INTERP_ALGO interp) {
    for (double eb : {1e-2, 1e-3, 1e-4, 1e-5, 1e-6}) {
        EXPECT_EQ(differing_values<T>(algo, interp, eb), 0u) << "eb = " << eb;
    }
}

}  // namespace

TEST(SZ3_Reconstruction, InterpLinear) {
    expect_identical<float>(SZ3::ALGO_INTERP, SZ3::INTERP_ALGO_LINEAR);
    expect_identical<double>(SZ3::ALGO_INTERP, SZ3::INTERP_ALGO_LINEAR);
}

TEST(SZ3_Reconstruction, InterpCubic) {
    expect_identical<float>(SZ3::ALGO_INTERP, SZ3::INTERP_ALGO_CUBIC);
    expect_identical<double>(SZ3::ALGO_INTERP, SZ3::INTERP_ALGO_CUBIC);
}

TEST(SZ3_Reconstruction, BioMD) {
    expect_identical<float>(SZ3::ALGO_BIOMD, SZ3::INTERP_ALGO_CUBIC);
    expect_identical<double>(SZ3::ALGO_BIOMD, SZ3::INTERP_ALGO_CUBIC);
}

// SZ_compress copies the input only for the algorithms that write their reconstruction into it.
TEST(SZ3_Reconstruction, SZCompressLeavesTheInputAsItWas) {
    for (auto algo : {SZ3::ALGO_LORENZO_REG, SZ3::ALGO_INTERP_LORENZO, SZ3::ALGO_INTERP, SZ3::ALGO_NOPRED,
                      SZ3::ALGO_BIOMD, SZ3::ALGO_BIOMDXTC}) {
        SZ3::Config conf(32, 341, 3);
        conf.cmprAlgo = algo;
        conf.errorBoundMode = SZ3::EB_ABS;
        conf.absErrorBound = 1e-3;
        std::vector<float> data(conf.num);
        for (size_t i = 0; i < conf.num; i++) data[i] = static_cast<float>(std::sin(0.01 * i) + 0.1 * (i % 3));
        const std::vector<float> before = data;
        size_t size = 0;
        delete[] SZ_compress(conf, data.data(), size);
        EXPECT_EQ(std::memcmp(before.data(), data.data(), data.size() * sizeof(float)), 0) << algo;
    }
}
