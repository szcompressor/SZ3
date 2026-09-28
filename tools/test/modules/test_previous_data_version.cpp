// Streams SZ3 wrote in data version 3.3.2, which this build still reads.
//
// tools/test/streams/<version>/ holds tools/sz3/testfloat_8_8_128.dat compressed by that release's sz3 CLI
// (-3 8 8 128, ABS 1e-3 unless the name says otherwise; the _omp ones with OpenMP on and 4 threads), and each
// decode must match that release's decode bit for bit.

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "SZ3/api/sz.hpp"
#include "gtest/gtest.h"

namespace {

std::vector<char> slurp(const std::string &name) {
    std::ifstream in(std::string(SZ3_TEST_STREAMS) + "/3.3.2/" + name + ".sz", std::ios::binary);
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
        {"biomdxtc", 0xe686c4ff7e13dca5ull},        {"interp", 0xe49ee9e883e32415ull},
        {"interp_omp", 0x1172c536c4c77855ull},      {"lorenzo_reg", 0x3f7750a23094223dull},
        {"lorenzo_reg_omp", 0xe32514629ac5fd3dull}, {"lorenzo_reg_rel1e-5", 0xe0885e9322ea80ccull},
        {"lossless", 0x51c93bc7f98c12a5ull},
    };
    for (const auto &s : streams) {
        const auto cmp = slurp(s.first);
        ASSERT_FALSE(cmp.empty()) << s.first;
        SZ3::Config conf;
        std::vector<float> dec(8 * 8 * 128);
        float *p = dec.data();
        SZ_decompress(conf, cmp.data(), cmp.size(), p);
        EXPECT_EQ(conf.sz3DataVer, versionInt("3.3.2")) << s.first;
        EXPECT_EQ(fnv1a(dec.data(), dec.size()), s.second) << s.first;
    }
}

// That version's ALGO_NOPRED and ALGO_BIOMD streams are not supported.
TEST(SZ3_PreviousDataVersion, NopredAndBioMDAreRefused) {
    for (const char *name : {"nopred", "biomd"}) {
        const auto cmp = slurp(name);
        ASSERT_FALSE(cmp.empty()) << name;
        SZ3::Config conf;
        std::vector<float> dec(8 * 8 * 128);
        float *p = dec.data();
        EXPECT_THROW(SZ_decompress(conf, cmp.data(), cmp.size(), p), std::invalid_argument) << name;
    }
}

// Outside SZ3_DATA_VER_OLDEST to SZ3_DATA_VER, a stream is refused with the way to read it.
TEST(SZ3_PreviousDataVersion, VersionsOutsideTheRangeAreRefused) {
    auto cmp = slurp("lorenzo_reg");
    ASSERT_FALSE(cmp.empty());
    for (const auto &v : {std::make_pair(versionInt("3.3.0"), "Use SZ3 v3.3.0"),
                          std::make_pair(versionInt(SZ3_DATA_VER) + (1u << 8), "Upgrade SZ3")}) {
        auto pos = reinterpret_cast<SZ3::uchar *>(cmp.data() + 4);  // after the magic number
        SZ3::write(v.first, pos);
        SZ3::Config conf;
        std::vector<float> dec(8 * 8 * 128);
        float *p = dec.data();
        try {
            SZ_decompress(conf, cmp.data(), cmp.size(), p);
            ADD_FAILURE() << versionStr(v.first) << " was read";
        } catch (const std::invalid_argument &e) {
            EXPECT_NE(std::string(e.what()).find(v.second), std::string::npos) << e.what();
        }
    }
}

// Compressing with a Config loaded from 3.3.2 data writes this build's data version, which is what the payload is.
TEST(SZ3_PreviousDataVersion, RecompressingWritesThisVersion) {
    const auto cmp = slurp("lorenzo_reg");
    ASSERT_FALSE(cmp.empty());
    SZ3::Config conf;
    std::unique_ptr<float[]> dec(SZ_decompress<float>(conf, cmp.data(), cmp.size()));
    size_t n = 0;
    std::unique_ptr<char[]> again(SZ_compress(conf, dec.get(), n));
    SZ3::Config back;
    std::unique_ptr<float[]> out(SZ_decompress<float>(back, again.get(), n));
    EXPECT_EQ(back.sz3DataVer, versionInt(SZ3_DATA_VER));
}

}  // namespace
