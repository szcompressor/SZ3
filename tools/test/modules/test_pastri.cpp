// Unit tests for PaSTRIDecomposition (pattern-scaling decomposition for two-electron integrals).

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/PaSTRIDecomposition.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Config.hpp"
#include "gtest/gtest.h"

using SZ3::PaSTRIDecomposition;

// ----- Helpers -------------------------------------------------------------

static size_t idxRange(int b) { return static_cast<size_t>((b + 1) * (b + 2) / 2); }

struct Geometry {
    size_t sbSize;
    size_t sbNum;
    size_t bSize;
};

static Geometry geometryOf(const std::array<int, 4>& bf) {
    Geometry g{};
    g.sbSize = idxRange(bf[2]) * idxRange(bf[3]);
    g.sbNum = idxRange(bf[0]) * idxRange(bf[1]);
    g.bSize = g.sbSize * g.sbNum;
    return g;
}

static void setAbsBound(SZ3::Config& conf, double eb) {
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
}

/**
 * Synthetic ERI-like data: within each block, every sub-block is a scaled copy of one pattern
 * (scale in [-1,1]) plus a small deviation. That is precisely the structure PaSTRI exploits.
 */
template <typename T>
static std::vector<T> makeERILike(const std::array<int, 4>& bf, size_t numBlocks, double deviation,
                                  unsigned seed = 12345) {
    const Geometry g = geometryOf(bf);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> patternDist(-1.0, 1.0);
    std::uniform_real_distribution<double> scaleDist(-1.0, 1.0);
    std::uniform_real_distribution<double> noiseDist(-1.0, 1.0);

    std::vector<T> out(numBlocks * g.bSize);
    for (size_t b = 0; b < numBlocks; b++) {
        // Pattern values span a few orders of magnitude, like real integral shells.
        std::vector<double> pattern(g.sbSize);
        for (size_t i = 0; i < g.sbSize; i++) {
            pattern[i] = patternDist(rng) * std::pow(10.0, -static_cast<double>(i % 4));
        }
        for (size_t s = 0; s < g.sbNum; s++) {
            const double scale = scaleDist(rng);
            for (size_t i = 0; i < g.sbSize; i++) {
                const double v = scale * pattern[i] + deviation * noiseDist(rng);
                out[b * g.bSize + s * g.sbSize + i] = static_cast<T>(v);
            }
        }
    }
    return out;
}

/// compress -> save -> load into a fresh instance -> decompress. Returns the max abs error.
template <typename T, SZ3::uint N>
static double roundtripSaveLoad(const SZ3::Config& conf, const std::vector<T>& original, double eb,
                                const std::array<int, 4>& bf, int radius = 32768) {
    std::vector<T> data = original;
    PaSTRIDecomposition<T, N> decomp(eb, bf, radius);
    auto bins = decomp.compress(conf, data.data());
    EXPECT_EQ(bins.size(), conf.num);
    EXPECT_EQ(decomp.get_out_range().first, 0);
    for (int bin : bins) {
        EXPECT_GE(bin, decomp.get_out_range().first);
        EXPECT_LE(bin, decomp.get_out_range().second);
    }

    std::vector<unsigned char> buf(64 + 32 * conf.num * sizeof(T) + (1u << 16));
    unsigned char* sp = buf.data();
    decomp.save(sp);
    const size_t saved = static_cast<size_t>(sp - buf.data());
    EXPECT_LE(saved, decomp.size_est()) << "size_est() must upper-bound what save() writes";

    // Construct with deliberately different parameters: load() must restore everything.
    PaSTRIDecomposition<T, N> decomp2(eb * 7.0, {0, 0, 0, 0}, 4096);
    const unsigned char* lp = buf.data();
    size_t remaining = saved;
    decomp2.load(lp, remaining);
    EXPECT_EQ(remaining, 0u) << "load() must consume exactly what save() wrote";
    EXPECT_EQ(static_cast<size_t>(lp - buf.data()), saved);

    std::vector<T> dec(conf.num, static_cast<T>(0));
    decomp2.decompress(conf, bins, dec.data());

    double maxErr = 0;
    for (size_t i = 0; i < conf.num; i++) {
        const double e = std::fabs(static_cast<double>(original[i]) - static_cast<double>(dec[i]));
        if (e > maxErr) maxErr = e;
    }
    // The input array must not be touched by compress().
    for (size_t i = 0; i < conf.num; i++) {
        EXPECT_EQ(data[i], original[i]) << "compress() modified the input at " << i;
    }
    return maxErr;
}

// ----- The central guarantee: the error bound holds for every point --------

TEST(SZ3_PaSTRITest, ErrorBoundHoldsForEveryPoint) {
    constexpr SZ3::uint N = 1;
    const std::array<std::array<int, 4>, 4> bfs = {{{1, 1, 1, 1}, {0, 1, 2, 1}, {2, 2, 1, 1}, {0, 0, 1, 3}}};
    const std::vector<double> ebs = {1e-1, 1e-2, 1e-3, 1e-5, 1e-7};

    for (const auto& bf : bfs) {
        const Geometry g = geometryOf(bf);
        for (double eb : ebs) {
            auto original = makeERILike<double>(bf, 7, 1e-4);
            SZ3::Config conf(original.size());
            setAbsBound(conf, eb);
            const double maxErr = roundtripSaveLoad<double, N>(conf, original, eb, bf);
            EXPECT_LE(maxErr, eb) << "bf={" << bf[0] << "," << bf[1] << "," << bf[2] << "," << bf[3]
                                  << "} bSize=" << g.bSize << " eb=" << eb;
        }
    }
}

TEST(SZ3_PaSTRITest, ErrorBoundHoldsForFloat) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {1, 1, 1, 1};
    for (double eb : {1e-2, 1e-4, 1e-6}) {
        auto original = makeERILike<float>(bf, 5, 1e-4);
        SZ3::Config conf(original.size());
        setAbsBound(conf, eb);
        EXPECT_LE((roundtripSaveLoad<float, N>(conf, original, eb, bf)), eb) << "eb=" << eb;
    }
}

// ----- Round trip / save-load ---------------------------------------------

TEST(SZ3_PaSTRITest, RoundTrip3D) {
    // N is only a template parameter for PaSTRI; the data is treated as a flat block sequence.
    constexpr SZ3::uint N = 3;
    const std::array<int, 4> bf = {1, 1, 1, 1};  // sbSize = sbNum = 9, bSize = 81
    const double eb = 1e-4;
    auto original = makeERILike<double>(bf, 8, 1e-5);
    ASSERT_EQ(original.size(), 8u * 81u);
    SZ3::Config conf(9, 9, 8);  // 648 points, same total as 8 blocks of 81
    ASSERT_EQ(conf.num, original.size());
    setAbsBound(conf, eb);
    EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb);
}

TEST(SZ3_PaSTRITest, SaveLoadPreservesModelExactly) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {2, 1, 1, 2};
    const double eb = 1e-4;
    auto original = makeERILike<double>(bf, 3, 1e-5);
    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);

    std::vector<double> data = original;
    PaSTRIDecomposition<double, N> a(eb, bf);
    auto bins = a.compress(conf, data.data());

    std::vector<unsigned char> buf(1u << 20);
    unsigned char* sp = buf.data();
    a.save(sp);
    const size_t saved = static_cast<size_t>(sp - buf.data());

    // Decompressing through the original instance and through a reloaded one must agree bit for bit.
    std::vector<double> decA(conf.num), decB(conf.num);
    a.decompress(conf, bins, decA.data());

    PaSTRIDecomposition<double, N> b(eb, bf);
    const unsigned char* lp = buf.data();
    size_t remaining = saved;
    b.load(lp, remaining);
    b.decompress(conf, bins, decB.data());

    for (size_t i = 0; i < conf.num; i++) {
        EXPECT_EQ(decA[i], decB[i]) << "mismatch at " << i;
    }
}

TEST(SZ3_PaSTRITest, LoadRejectsBadUid) {
    const std::array<int, 4> bf = {1, 1, 1, 1};
    PaSTRIDecomposition<double, 1> a(1e-3, bf);
    std::vector<unsigned char> buf(1u << 12, 0);
    unsigned char* sp = buf.data();
    a.save(sp);
    buf[0] ^= 0xFF;

    PaSTRIDecomposition<double, 1> b(1e-3, bf);
    const unsigned char* lp = buf.data();
    size_t remaining = static_cast<size_t>(sp - buf.data());
    EXPECT_THROW(b.load(lp, remaining), std::invalid_argument);
}

// ----- Contracts -----------------------------------------------------------

TEST(SZ3_PaSTRITest, OutRangeStartsAtZero) {
    for (int radius : {2, 1024, 32768}) {
        PaSTRIDecomposition<float, 1> decomp(1e-3, {1, 1, 1, 1}, radius);
        EXPECT_EQ(decomp.get_out_range().first, 0);
        EXPECT_EQ(decomp.get_out_range().second, radius * 2);
    }
}

TEST(SZ3_PaSTRITest, GeometryFromBasisFunctions) {
    PaSTRIDecomposition<double, 1> d(1e-3, {1, 1, 1, 1});
    EXPECT_EQ(d.get_sub_block_size(), 9u);
    EXPECT_EQ(d.get_sub_block_num(), 9u);
    EXPECT_EQ(d.get_block_size(), 81u);

    PaSTRIDecomposition<double, 1> d2(1e-3, {0, 0, 2, 3});  // 1*1 sub-blocks of 6*10
    EXPECT_EQ(d2.get_sub_block_size(), 60u);
    EXPECT_EQ(d2.get_sub_block_num(), 1u);
    EXPECT_EQ(d2.get_block_size(), 60u);
}

TEST(SZ3_PaSTRITest, RejectsInvalidConstructorArguments) {
    EXPECT_THROW((PaSTRIDecomposition<double, 1>(0.0, {1, 1, 1, 1})), std::invalid_argument);
    EXPECT_THROW((PaSTRIDecomposition<double, 1>(-1e-3, {1, 1, 1, 1})), std::invalid_argument);
    EXPECT_THROW((PaSTRIDecomposition<double, 1>(1e-3, {1, 1, 1, 1}, 1)), std::invalid_argument);
    EXPECT_THROW((PaSTRIDecomposition<double, 1>(1e-3, {-1, 1, 1, 1})), std::invalid_argument);
}

// ----- Degenerate inputs ---------------------------------------------------

TEST(SZ3_PaSTRITest, AllZeroBlockIsExact) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {1, 1, 1, 1};
    const Geometry g = geometryOf(bf);
    const double eb = 1e-3;
    std::vector<double> original(4 * g.bSize, 0.0);
    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);

    std::vector<double> data = original;
    PaSTRIDecomposition<double, N> decomp(eb, bf);
    auto bins = decomp.compress(conf, data.data());
    std::vector<double> dec(conf.num, 1.0);
    decomp.decompress(conf, bins, dec.data());
    for (size_t i = 0; i < conf.num; i++) {
        EXPECT_EQ(dec[i], 0.0) << "all-zero block must reconstruct exactly at " << i;
    }
}

TEST(SZ3_PaSTRITest, ZeroExtremumBlockMixedWithRealBlocks) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {1, 1, 1, 1};
    const Geometry g = geometryOf(bf);
    const double eb = 1e-4;
    // Block 0 real, block 1 all zeros (extremum is zero), block 2 real.
    auto original = makeERILike<double>(bf, 3, 1e-5);
    for (size_t i = 0; i < g.bSize; i++) original[g.bSize + i] = 0.0;

    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);
    EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb);
}

TEST(SZ3_PaSTRITest, SingleSubBlock) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {0, 0, 2, 2};  // sbNum == 1, sbSize == 36
    const Geometry g = geometryOf(bf);
    ASSERT_EQ(g.sbNum, 1u);
    const double eb = 1e-5;
    auto original = makeERILike<double>(bf, 6, 1e-6);
    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);
    EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb);
}

TEST(SZ3_PaSTRITest, SmallestBlockGeometry) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {0, 0, 0, 0};  // bSize == 1
    const double eb = 1e-3;
    PaSTRIDecomposition<double, N> probe(eb, bf);
    ASSERT_EQ(probe.get_block_size(), 1u);
    auto original = makeERILike<double>(bf, 40, 1e-4);
    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);
    EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb);
}

TEST(SZ3_PaSTRITest, ConstantAndHugeValues) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {1, 1, 1, 1};
    const Geometry g = geometryOf(bf);
    const double eb = 1e-3;
    std::vector<double> original(2 * g.bSize);
    for (size_t i = 0; i < g.bSize; i++) original[i] = 42.0;            // constant block
    for (size_t i = 0; i < g.bSize; i++) original[g.bSize + i] = 1e12;  // far outside the radius
    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);
    EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb);
}

TEST(SZ3_PaSTRITest, WideBitWidthsFromVeryTightErrorBound) {
    // eb = 1e-15 on O(1) data drives patternBits to ~50, so scalesQ*patternQ needs ~100 bits.
    // The reference implementation forms that product in int64 and would overflow; this module
    // multiplies in the double domain instead. The error bound must still hold either way.
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {1, 1, 1, 1};
    for (double eb : {1e-12, 1e-15}) {
        auto original = makeERILike<double>(bf, 6, eb * 10);
        SZ3::Config conf(original.size());
        setAbsBound(conf, eb);
        EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb) << "eb=" << eb;
    }
}

// ----- Data without the pattern structure ----------------------------------

TEST(SZ3_PaSTRITest, RandomNoiseStillRespectsErrorBound) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {1, 1, 1, 1};
    const Geometry g = geometryOf(bf);
    const double eb = 1e-3;
    std::mt19937 rng(999);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<double> original(10 * g.bSize);
    for (auto& v : original) v = dist(rng);

    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);
    EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb);
}

// ----- Tail handling (conf.num not a multiple of bSize) --------------------

TEST(SZ3_PaSTRITest, TailIsHandledAndBounded) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {1, 1, 1, 1};
    const Geometry g = geometryOf(bf);
    const double eb = 1e-4;
    for (size_t extra : {size_t(1), size_t(7), g.bSize - 1}) {
        auto original = makeERILike<double>(bf, 4, 1e-5);
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> dist(-0.5, 0.5);
        for (size_t i = 0; i < extra; i++) original.push_back(dist(rng));
        ASSERT_NE(original.size() % g.bSize, 0u);

        SZ3::Config conf(original.size());
        setAbsBound(conf, eb);
        EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb) << "extra=" << extra;
    }
}

TEST(SZ3_PaSTRITest, ShorterThanOneBlock) {
    constexpr SZ3::uint N = 1;
    const std::array<int, 4> bf = {1, 1, 1, 1};
    const double eb = 1e-4;
    std::vector<double> original = {0.1, -0.25, 3.5, 0.0, -1e-6};
    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);
    EXPECT_LE((roundtripSaveLoad<double, N>(conf, original, eb, bf)), eb);
}

// ----- Full pipeline: decomposition + Huffman + Zstd ------------------------

template <typename T>
static double pipelineRatio(const std::vector<T>& original, double eb, const std::array<int, 4>& bf,
                            const char* label) {
    constexpr SZ3::uint N = 1;
    SZ3::Config conf(original.size());
    setAbsBound(conf, eb);

    std::vector<T> data = original;
    const size_t cmpCap = original.size() * sizeof(T) * 2 + (1u << 16);
    std::vector<unsigned char> cmp(cmpCap);

    auto sz = SZ3::make_compressor_sz_generic<T, N>(PaSTRIDecomposition<T, N>(eb, bf), SZ3::HuffmanEncoder<int>(),
                                                    SZ3::Lossless_zstd());
    const size_t cmpSize = sz->compress(conf, data.data(), cmp.data(), cmpCap);

    std::vector<T> dec(conf.num, static_cast<T>(0));
    auto sz2 = SZ3::make_compressor_sz_generic<T, N>(PaSTRIDecomposition<T, N>(eb, bf), SZ3::HuffmanEncoder<int>(),
                                                     SZ3::Lossless_zstd());
    sz2->decompress(conf, cmp.data(), cmpSize, dec.data());

    double maxErr = 0;
    for (size_t i = 0; i < conf.num; i++) {
        maxErr = std::max(maxErr, std::fabs(static_cast<double>(original[i]) - static_cast<double>(dec[i])));
    }
    EXPECT_LE(maxErr, eb) << label;

    const double ratio = static_cast<double>(original.size() * sizeof(T)) / static_cast<double>(cmpSize);
    printf("[PaSTRI pipeline] %-22s eb=%-8.1e raw=%zuB cmp=%zuB CR=%.2fx maxErr=%.3e\n", label, eb,
           original.size() * sizeof(T), cmpSize, ratio, maxErr);
    return ratio;
}

TEST(SZ3_PaSTRITest, PipelineCompressionRatioOnPatternedData) {
    const std::array<int, 4> bf = {1, 1, 1, 1};
    auto patterned = makeERILike<double>(bf, 2000, 1e-6);

    const double crTight = pipelineRatio<double>(patterned, 1e-6, bf, "patterned");
    const double crLoose = pipelineRatio<double>(patterned, 1e-4, bf, "patterned");
    EXPECT_GT(crTight, 1.0);
    EXPECT_GT(crLoose, crTight) << "a looser bound must not compress worse";

    // Noise has no pattern structure: the bound must still hold (checked inside), CR is just poor.
    std::mt19937 rng(4242);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<double> noise(patterned.size());
    for (auto& v : noise) v = dist(rng);
    const double crNoise = pipelineRatio<double>(noise, 1e-4, bf, "random noise");
    EXPECT_LT(crNoise, crLoose) << "structured data must compress better than noise";
}
