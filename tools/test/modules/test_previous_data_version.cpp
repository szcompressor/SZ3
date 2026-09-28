// Streams SZ3 wrote in the previous data version, SZ3_DATA_VER_PREV, which this build still reads.
//
// tools/test/streams/<version>/ holds tools/sz3/testfloat_8_8_128.dat compressed by that release's sz3 CLI
// (-3 8 8 128, ABS 1e-3 unless the name says otherwise; the _omp ones with OpenMP on and 4 threads), and each
// decode must match that release's decode bit for bit.

#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "SZ3/api/sz.hpp"
#include "gtest/gtest.h"

namespace {

std::vector<char> slurp(const std::string &name) {
    std::ifstream in(std::string(SZ3_TEST_STREAMS) + "/" SZ3_DATA_VER_PREV "/" + name + ".sz", std::ios::binary);
    return std::vector<char>(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// FNV-1a over the decoded bytes.
uint64_t fnv1a(const float *data, size_t n) {
    uint64_t h = 0xcbf29ce484222325ull;
    const auto *b = reinterpret_cast<const unsigned char *>(data);
    for (size_t i = 0; i < n * sizeof(float); i++) h = (h ^ b[i]) * 0x100000001b3ull;
    return h;
}

TEST(SZ3_PreviousDataVersion, DecodesAsItsReleaseDid) {
    const std::pair<const char *, uint64_t> streams[] = {
        {"interp", 0xe49ee9e883e32415ull},
        {"interp_omp", 0x1172c536c4c77855ull},
        {"lorenzo_reg", 0x3f7750a23094223dull},
        {"lorenzo_reg_omp", 0xe32514629ac5fd3dull},
        {"lorenzo_reg_rel1e-5", 0xe0885e9322ea80ccull},
        {"lossless", 0x51c93bc7f98c12a5ull},
        {"nopred", 0xe686c4ff7e13dca5ull},
    };
    for (const auto &s : streams) {
        const auto cmp = slurp(s.first);
        ASSERT_FALSE(cmp.empty()) << s.first;
        SZ3::Config conf;
        std::vector<float> dec(8 * 8 * 128);
        float *p = dec.data();
        SZ_decompress(conf, cmp.data(), cmp.size(), p);
        EXPECT_EQ(conf.sz3DataVer, versionInt(SZ3_DATA_VER_PREV)) << s.first;
        EXPECT_EQ(fnv1a(dec.data(), dec.size()), s.second) << s.first;
    }
}

// That version's ALGO_BIOMD streams used a Huffman format this build no longer has.
TEST(SZ3_PreviousDataVersion, BioMDIsRefused) {
    const auto cmp = slurp("biomd");
    ASSERT_FALSE(cmp.empty());
    SZ3::Config conf;
    std::vector<float> dec(8 * 8 * 128);
    float *p = dec.data();
    EXPECT_THROW(SZ_decompress(conf, cmp.data(), cmp.size(), p), std::invalid_argument);
}

}  // namespace
