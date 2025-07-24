#include <cmath>
#include <cstdint>
#include <random>

#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "gtest/gtest.h"

template <class Lossless>
void runFunctionalTest() {
    Lossless lossless;
    size_t N = 1000;
    std::vector<SZ3::uchar> src(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 255);
    for (size_t i = 0; i < N; i++) {
        src[i] = static_cast<SZ3::uchar>(dis(gen));
    }
    std::vector<SZ3::uchar> dst(ZSTD_compressBound(src.size())+ sizeof(size_t));
    size_t compressedSize = lossless.compress(src.data(), src.size(), dst.data(), dst.size());

    std::vector<SZ3::uchar> decompressed(N);
    SZ3::uchar* decompressed_pos = decompressed.data();
    size_t decompressedSize;
    lossless.decompress(dst.data(), compressedSize, decompressed_pos, decompressedSize);

    EXPECT_EQ(decompressedSize, src.size());
    EXPECT_EQ(std::vector<SZ3::uchar>(decompressed.data(), decompressed.data() + decompressedSize), src);
}

template <typename Lossless>
void runAllTest() {
    runFunctionalTest<Lossless>();
}

TEST(LosslessTest, LosslessZstd) { runAllTest<SZ3::Lossless_zstd>(); }
TEST(LosslessTest, LosslessBypass) { runAllTest<SZ3::Lossless_bypass>(); }
