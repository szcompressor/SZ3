// Tests for the decomposed MGARD path: MGARDTransform (preprocessor) +
// geometric_level_ebs (error budget policy) + multilevel_quantize/recover
// (level walk) + MultiLevelDecomposition (the composition of the three).
//
// The load-bearing claim is that `MultiLevelDecomposition` wired with its
// defaults is numerically identical to the monolithic `MGARDFusedDecomposition`,
// which stays the supported path.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <vector>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/MGARDFusedDecomposition.hpp"
#include "SZ3/decomposition/MultiLevelDecomposition.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/preprocessor/MGARDTransform.hpp"
#include "SZ3/quantizer/LevelQuantizer.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MultiLevelErrorBound.hpp"
#include "SZ3/utils/MultiLevelQuantization.hpp"
#include "gtest/gtest.h"

// ----- Shared helpers ------------------------------------------------------

template <typename T, typename F>
static std::vector<T> make1D(size_t dim, F f) {
    std::vector<T> v(dim);
    for (size_t i = 0; i < dim; i++) v[i] = static_cast<T>(f(i));
    return v;
}
template <typename T, typename F>
static std::vector<T> make2D(size_t H, size_t W, F f) {
    std::vector<T> v(H * W);
    for (size_t y = 0; y < H; y++)
        for (size_t x = 0; x < W; x++) v[y * W + x] = static_cast<T>(f(y, x));
    return v;
}
template <typename T, typename F>
static std::vector<T> make3D(size_t Z, size_t Y, size_t X, F f) {
    std::vector<T> v(Z * Y * X);
    for (size_t z = 0; z < Z; z++)
        for (size_t y = 0; y < Y; y++)
            for (size_t x = 0; x < X; x++) v[(z * Y + y) * X + x] = static_cast<T>(f(z, y, x));
    return v;
}

static void setAbsBound(SZ3::Config& conf, double eb) {
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
}

static double maxAbsErr(const std::vector<float>& a, const std::vector<float>& b) {
    double m = 0;
    for (size_t i = 0; i < a.size(); i++) m = std::max(m, std::fabs(static_cast<double>(a[i] - b[i])));
    return m;
}

// The default composition: MGARD multigrid + geometric schedule + LinearQuantizer.
using DefaultSplit1D = SZ3::MultiLevelDecomposition<float, 1>;
using DefaultSplit2D = SZ3::MultiLevelDecomposition<float, 2>;
using DefaultSplit3D = SZ3::MultiLevelDecomposition<float, 3>;

static std::function<SZ3::LinearQuantizer<float>(double)> linearFactory(int radius) {
    return [radius](double level_eb) { return SZ3::LinearQuantizer<float>(level_eb, radius); };
}

// ----- 1. MGARDTransform round-trip ----------------------------------------

TEST(SZ3_MGARDComposable, TransformRoundTrip1D) {
    constexpr size_t dim = 257;
    const std::vector<size_t> dims{dim};
    auto original = make1D<float>(dim, [](size_t i) { return std::sin(0.05 * i) + 0.3 * std::cos(0.13 * i); });
    std::vector<float> data = original;

    SZ3::MGARDTransform<float, 1> transform;
    const size_t level = SZ3::MGARDTransform<float, 1>::default_target_level(dims);
    EXPECT_EQ(level, 7u);  // floor(log2(257)) - 1

    transform.preprocess(data.data(), dims, level);
    EXPECT_GT(maxAbsErr(original, data), 1e-3) << "forward transform should actually change the buffer";
    transform.postprocess(data.data(), dims, level);
    EXPECT_LT(maxAbsErr(original, data), 1e-4) << "inverse transform must undo the forward transform";
}

TEST(SZ3_MGARDComposable, TransformRoundTrip2D) {
    constexpr size_t dim = 33;
    const std::vector<size_t> dims{dim, dim};
    auto original = make2D<float>(dim, dim, [](size_t y, size_t x) { return std::sin(0.05 * x) * std::cos(0.03 * y); });
    std::vector<float> data = original;

    SZ3::MGARDTransform<float, 2> transform;
    const size_t level = SZ3::MGARDTransform<float, 2>::default_target_level(dims);
    transform.preprocess(data.data(), dims, level);
    transform.postprocess(data.data(), dims, level);
    EXPECT_LT(maxAbsErr(original, data), 1e-4);
}

TEST(SZ3_MGARDComposable, TransformRoundTrip3D) {
    constexpr size_t dim = 17;
    const std::vector<size_t> dims{dim, dim, dim};
    auto original = make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1 * z) + 0.5 * std::cos(0.07 * y) + 0.25 * std::sin(0.03 * x);
    });
    std::vector<float> data = original;

    SZ3::MGARDTransform<float, 3> transform;
    const size_t level = SZ3::MGARDTransform<float, 3>::default_target_level(dims);
    EXPECT_EQ(level, 3u);
    transform.preprocess(data.data(), dims, level);
    transform.postprocess(data.data(), dims, level);
    EXPECT_LT(maxAbsErr(original, data), 1e-4);
}

TEST(SZ3_MGARDComposable, TransformRejectsBadDims) {
    SZ3::MGARDTransform<float, 2> transform;
    std::vector<float> data(16, 0.0f);
    EXPECT_THROW(transform.preprocess(data.data(), std::vector<size_t>{16}, 1), std::invalid_argument);
    EXPECT_THROW(transform.preprocess(nullptr, std::vector<size_t>{4, 4}, 1), std::invalid_argument);
}

// ----- 2. Level walk covers every element exactly once ----------------------

template <SZ3::uint N>
static void expectLevelWalkPartitions(const std::vector<size_t>& dims) {
    const size_t level = SZ3::MGARDTransform<float, N>::default_target_level(dims);
    const auto levels = SZ3::MGARDTransform<float, N>::level_dims(dims, level);
    size_t num = 1;
    for (size_t d : dims) num *= d;

    std::vector<int> visits(num, 0);
    const std::vector<size_t> no_coarse(dims.size(), 0);
    for (size_t l = 0; l < levels.size(); l++) {
        SZ3::for_each_level_coefficient<N>(dims, (l == 0) ? no_coarse : levels[l - 1], levels[l], [&](size_t offset) {
            ASSERT_LT(offset, num);
            visits[offset]++;
        });
    }
    for (size_t i = 0; i < num; i++) {
        EXPECT_EQ(visits[i], 1) << "offset " << i << " visited " << visits[i] << " times";
    }
}

TEST(SZ3_MGARDComposable, LevelWalkPartitions1D) { expectLevelWalkPartitions<1>({257}); }
TEST(SZ3_MGARDComposable, LevelWalkPartitions2D) { expectLevelWalkPartitions<2>({33, 33}); }
TEST(SZ3_MGARDComposable, LevelWalkPartitions3D) { expectLevelWalkPartitions<3>({17, 17, 17}); }
TEST(SZ3_MGARDComposable, LevelWalkPartitionsAnisotropic3D) { expectLevelWalkPartitions<3>({9, 17, 12}); }

// ----- 3. The level-eb schedule in isolation --------------------------------

TEST(SZ3_MGARDComposable, GeometricLevelEbsBudget) {
    const double eb = 1e-2;
    const double growth = SZ3::MGARDTransform<float, 3>::level_eb_growth();
    const double amplification = SZ3::MGARDTransform<float, 3>::level_eb_amplification();

    for (size_t L : {size_t(0), size_t(1), size_t(3), size_t(8)}) {
        auto ebs = SZ3::geometric_level_ebs(eb, L, growth, amplification);
        ASSERT_EQ(ebs.size(), L + 1);

        // Coarse levels are tighter; the ratio between neighbours is exactly `growth`.
        double sum = ebs[0];
        for (size_t l = 1; l <= L; l++) {
            EXPECT_NEAR(ebs[l] / ebs[l - 1], growth, 1e-12) << "L=" << L << " l=" << l;
            EXPECT_GT(ebs[l], ebs[l - 1]);
            sum += ebs[l];
        }
        // The budget invariant: the per-level bounds sum to eb / amplification.
        EXPECT_NEAR(sum, eb / amplification, 1e-12 * eb / amplification) << "L=" << L;
        EXPECT_GT(ebs[0], 0.0);
    }
}

TEST(SZ3_MGARDComposable, GeometricLevelEbsMatchesMonolithicSchedule) {
    // Same closed form the monolithic MGARDFusedDecomposition::build_level_ebs uses.
    const double eb = 3.7e-3;
    const double c = std::sqrt(8.0);
    const double C2 = 1.0 + 3.0 * std::sqrt(3.0) / 4.0;
    for (size_t L = 0; L <= 8; L++) {
        const double cc = (1.0 - c) / (1.0 - std::pow(c, static_cast<double>(L + 1)));
        std::vector<double> expected(L + 1);
        expected[0] = cc * eb / C2;
        for (size_t l = 1; l <= L; l++) expected[l] = expected[l - 1] * c;

        auto actual = SZ3::geometric_level_ebs(eb, L, c, C2);
        ASSERT_EQ(actual.size(), expected.size());
        for (size_t l = 0; l <= L; l++) {
            EXPECT_DOUBLE_EQ(actual[l], expected[l]) << "L=" << L << " l=" << l;
        }
    }
}

TEST(SZ3_MGARDComposable, GeometricLevelEbsDegenerateGrowth) {
    // growth == 1 is the limit of the geometric series: an even split.
    auto ebs = SZ3::geometric_level_ebs(1.0, 3, 1.0, 2.0);
    ASSERT_EQ(ebs.size(), 4u);
    for (double e : ebs) EXPECT_DOUBLE_EQ(e, 1.0 / (2.0 * 4.0));
}

TEST(SZ3_MGARDComposable, GeometricLevelEbsRejectsBadArgs) {
    EXPECT_THROW(SZ3::geometric_level_ebs(0.0, 3, 2.0, 2.0), std::invalid_argument);
    EXPECT_THROW(SZ3::geometric_level_ebs(1.0, 3, -1.0, 2.0), std::invalid_argument);
    EXPECT_THROW(SZ3::geometric_level_ebs(1.0, 3, 2.0, 0.0), std::invalid_argument);
}

// ----- 4. Numerical equivalence with the monolithic MGARDFusedDecomposition ------

// Runs both modules over the same input/eb and requires bit-identical bins and
// bit-identical reconstructions. `out_bins`, when given, receives the bin stream.
template <SZ3::uint N, class Split>
static void expectEquivalence(SZ3::Config& conf, const std::vector<float>& original, double eb, int radius,
                              std::vector<int>* out_bins = nullptr) {
    std::vector<float> data_mono = original;
    SZ3::MGARDFusedDecomposition<float, N> mono(eb, radius);
    auto bins_mono = mono.compress(conf, data_mono.data());

    std::vector<float> data_split = original;
    Split split(eb, linearFactory(radius));
    auto bins_split = split.compress(conf, data_split.data());

    ASSERT_EQ(bins_mono.size(), conf.num);
    ASSERT_EQ(bins_split.size(), bins_mono.size());
    for (size_t i = 0; i < bins_mono.size(); i++) {
        ASSERT_EQ(bins_split[i], bins_mono[i]) << "bin " << i << " differs";
    }
    EXPECT_EQ(split.get_out_range(), mono.get_out_range());
    EXPECT_EQ(split.get_out_range().first, 0);

    std::vector<float> dec_mono(conf.num), dec_split(conf.num);
    mono.decompress(conf, bins_mono, dec_mono.data());
    split.decompress(conf, bins_split, dec_split.data());
    for (size_t i = 0; i < conf.num; i++) {
        ASSERT_EQ(dec_split[i], dec_mono[i]) << "reconstruction differs at " << i;
    }
    EXPECT_LE(maxAbsErr(original, dec_split), eb);
    if (out_bins != nullptr) {
        *out_bins = bins_split;
    }
}

TEST(SZ3_MGARDComposable, EquivalentToMonolithic1D) {
    constexpr size_t dim = 257;
    const double eb = 1e-3;
    SZ3::Config conf(dim);
    setAbsBound(conf, eb);
    auto original = make1D<float>(dim, [](size_t i) { return std::sin(0.05 * i) + 0.3 * std::cos(0.13 * i); });
    expectEquivalence<1, DefaultSplit1D>(conf, original, eb, 32768);
}

TEST(SZ3_MGARDComposable, EquivalentToMonolithic2D) {
    constexpr size_t dim = 33;
    const double eb = 1e-3;
    SZ3::Config conf(dim, dim);
    setAbsBound(conf, eb);
    auto original = make2D<float>(dim, dim, [](size_t y, size_t x) { return std::sin(0.05 * x) * std::cos(0.03 * y); });
    expectEquivalence<2, DefaultSplit2D>(conf, original, eb, 32768);
}

TEST(SZ3_MGARDComposable, EquivalentToMonolithic3D) {
    constexpr size_t dim = 17;
    const double eb = 1e-2;
    SZ3::Config conf(dim, dim, dim);
    setAbsBound(conf, eb);
    auto original = make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1 * z) + 0.5 * std::cos(0.07 * y) + 0.25 * std::sin(0.03 * x);
    });
    expectEquivalence<3, DefaultSplit3D>(conf, original, eb, 32768);
}

TEST(SZ3_MGARDComposable, EquivalentToMonolithicAnisotropic3D) {
    const double eb = 5e-3;
    SZ3::Config conf(9, 17, 12);
    setAbsBound(conf, eb);
    auto original = make3D<float>(9, 17, 12, [](size_t z, size_t y, size_t x) {
        return std::exp(-0.01 * (z * z + y * y)) + 0.2 * std::sin(0.31 * x);
    });
    expectEquivalence<3, DefaultSplit3D>(conf, original, eb, 32768);
}

TEST(SZ3_MGARDComposable, EquivalentToMonolithicNoisyData3D) {
    // Noise forces LinearQuantizer down its unpredictable-value path on the
    // coarse levels, so this also checks the unpred lists line up.
    constexpr size_t dim = 17;
    const double eb = 1e-4;
    SZ3::Config conf(dim, dim, dim);
    setAbsBound(conf, eb);
    uint32_t seed = 12345;
    auto rnd = [&seed]() {
        seed = seed * 1664525u + 1013904223u;
        return static_cast<double>(seed >> 8) / static_cast<double>(1u << 24);
    };
    auto original = make3D<float>(dim, dim, dim, [&](size_t z, size_t y, size_t x) {
        return std::sin(0.1 * z) * std::cos(0.13 * y) + 0.5 * rnd() + 0.01 * x;
    });
    std::vector<int> bins;
    expectEquivalence<3, DefaultSplit3D>(conf, original, eb, 32768, &bins);
    EXPECT_GT(std::count(bins.begin(), bins.end(), 0), 0)
        << "expected LinearQuantizer to fall back to unpredictable values here";
}

// ----- 5. save/load round-trip through the composed decomposition -----------

TEST(SZ3_MGARDComposable, SaveLoadRoundTrip3D) {
    constexpr size_t dim = 17;
    const double eb = 1e-2;
    const int radius = 32768;
    SZ3::Config conf(dim, dim, dim);
    setAbsBound(conf, eb);
    auto original = make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1 * z) + 0.5 * std::cos(0.07 * y) + 0.25 * std::sin(0.03 * x);
    });

    std::vector<float> data = original;
    DefaultSplit3D writer(eb, linearFactory(radius));
    auto bins = writer.compress(conf, data.data());

    std::vector<unsigned char> buf(1u << 22);
    unsigned char* sp = buf.data();
    writer.save(sp);
    const size_t saved = static_cast<size_t>(sp - buf.data());
    EXPECT_LE(saved, writer.size_est()) << "size_est must not under-report the serialised size";

    DefaultSplit3D reader(eb, linearFactory(radius));
    const unsigned char* lp = buf.data();
    size_t remaining = saved;
    reader.load(lp, remaining);
    EXPECT_EQ(static_cast<size_t>(lp - buf.data()), saved) << "load must consume exactly what save wrote";
    EXPECT_EQ(remaining, 0u);

    std::vector<float> dec(conf.num);
    reader.decompress(conf, bins, dec.data());
    EXPECT_LE(maxAbsErr(original, dec), eb);
}

// ----- 6. Full pipeline equivalence (decomposition + Huffman + Zstd) --------

TEST(SZ3_MGARDComposable, PipelineEquivalence3D) {
    constexpr SZ3::uint N = 3;
    constexpr size_t dim = 17;
    const double eb = 1e-2;
    const int radius = 32768;
    SZ3::Config conf(dim, dim, dim);
    setAbsBound(conf, eb);
    auto original = make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1 * z) + 0.5 * std::cos(0.07 * y) + 0.25 * std::sin(0.03 * x);
    });

    const size_t cap = conf.num * sizeof(float) * 2 + 4096;

    std::vector<unsigned char> cmp_mono(cap);
    std::vector<float> in_mono = original;
    auto mono = SZ3::make_compressor_sz_generic<float, N>(SZ3::MGARDFusedDecomposition<float, N>(eb, radius),
                                                          SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
    const size_t size_mono = mono->compress(conf, in_mono.data(), cmp_mono.data(), cap);

    std::vector<unsigned char> cmp_split(cap);
    std::vector<float> in_split = original;
    auto split = SZ3::make_compressor_sz_generic<float, N>(DefaultSplit3D(eb, linearFactory(radius)),
                                                           SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
    const size_t size_split = split->compress(conf, in_split.data(), cmp_split.data(), cap);

    // The two streams differ only by the monolithic module's extra serialised
    // `radius` field (its quantizers already carry their own radius).
    EXPECT_NEAR(static_cast<double>(size_split), static_cast<double>(size_mono), 16.0);

    std::vector<float> dec_mono(conf.num), dec_split(conf.num);
    auto mono_dec = SZ3::make_compressor_sz_generic<float, N>(SZ3::MGARDFusedDecomposition<float, N>(eb, radius),
                                                              SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
    mono_dec->decompress(conf, cmp_mono.data(), size_mono, dec_mono.data());
    auto split_dec = SZ3::make_compressor_sz_generic<float, N>(DefaultSplit3D(eb, linearFactory(radius)),
                                                               SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
    split_dec->decompress(conf, cmp_split.data(), size_split, dec_split.data());

    for (size_t i = 0; i < conf.num; i++) {
        ASSERT_EQ(dec_split[i], dec_mono[i]) << "pipeline reconstruction differs at " << i;
    }
    EXPECT_LE(maxAbsErr(original, dec_split), eb);
}

// ----- 7. The point of the split: a quantizer other than LinearQuantizer ----

TEST(SZ3_MGARDComposable, PairsWithLevelQuantizer3D) {
    constexpr SZ3::uint N = 3;
    constexpr size_t dim = 17;
    const double eb = 1e-2;
    const int radius = 32768;
    SZ3::Config conf(dim, dim, dim);
    setAbsBound(conf, eb);
    auto original = make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1 * z) + 0.5 * std::cos(0.07 * y) + 0.25 * std::sin(0.03 * x);
    });

    using Quant = SZ3::LevelQuantizer<float>;
    using Split = SZ3::MultiLevelDecomposition<float, N, Quant>;
    // Both sides build the quantizer from the same factory, so the curve matches on load.
    auto factory = [radius](double level_eb) { return Quant(level_eb, radius, SZ3::LevelCurve::Quadratic); };

    std::vector<float> data = original;
    Split split(eb, factory);
    auto bins = split.compress(conf, data.data());
    ASSERT_EQ(bins.size(), conf.num);
    EXPECT_EQ(split.get_out_range().first, 0);

    // save/load so the non-default quantizer's serialisation is exercised too.
    std::vector<unsigned char> buf(1u << 22);
    unsigned char* sp = buf.data();
    split.save(sp);
    const size_t saved = static_cast<size_t>(sp - buf.data());

    Split reader(eb, factory);
    const unsigned char* lp = buf.data();
    size_t remaining = saved;
    reader.load(lp, remaining);
    EXPECT_EQ(remaining, 0u);

    std::vector<float> dec(conf.num);
    reader.decompress(conf, bins, dec.data());
    EXPECT_LE(maxAbsErr(original, dec), eb) << "MGARD transform + LevelQuantizer must honour the eb";

    // It really is a different quantizer: the bin stream must not match the
    // LinearQuantizer one.
    std::vector<float> data_lin = original;
    DefaultSplit3D linear(eb, linearFactory(radius));
    auto bins_linear = linear.compress(conf, data_lin.data());
    ASSERT_EQ(bins_linear.size(), bins.size());
    EXPECT_NE(bins, bins_linear);
}

TEST(SZ3_MGARDComposable, PairsWithLevelQuantizerInFullPipeline2D) {
    constexpr SZ3::uint N = 2;
    constexpr size_t dim = 33;
    const double eb = 1e-2;
    const int radius = 32768;
    SZ3::Config conf(dim, dim);
    setAbsBound(conf, eb);
    auto original = make2D<float>(dim, dim, [](size_t y, size_t x) { return std::sin(0.05 * x) * std::cos(0.03 * y); });

    using Quant = SZ3::LevelQuantizer<float>;
    using Split = SZ3::MultiLevelDecomposition<float, N, Quant>;
    // Both sides build the quantizer from the same factory, so the curve matches on load.
    auto factory = [radius](double level_eb) { return Quant(level_eb, radius, SZ3::LevelCurve::Quadratic); };

    const size_t cap = conf.num * sizeof(float) * 2 + 4096;
    std::vector<unsigned char> cmp(cap);
    std::vector<float> in = original;
    auto comp =
        SZ3::make_compressor_sz_generic<float, N>(Split(eb, factory), SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
    const size_t cmp_size = comp->compress(conf, in.data(), cmp.data(), cap);
    EXPECT_GT(cmp_size, 0u);

    std::vector<float> dec(conf.num);
    auto decomp =
        SZ3::make_compressor_sz_generic<float, N>(Split(eb, factory), SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
    decomp->decompress(conf, cmp.data(), cmp_size, dec.data());
    EXPECT_LE(maxAbsErr(original, dec), eb);
}

// ----- 8. Level-count override ---------------------------------------------

TEST(SZ3_MGARDComposable, TargetLevelOverrideIsHonoured) {
    constexpr size_t dim = 17;
    const double eb = 1e-2;
    const int radius = 32768;
    SZ3::Config conf(dim, dim, dim);
    setAbsBound(conf, eb);
    auto original = make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1 * z) + 0.5 * std::cos(0.07 * y) + 0.25 * std::sin(0.03 * x);
    });

    for (size_t level : {size_t(0), size_t(1), size_t(2)}) {
        std::vector<float> data = original;
        DefaultSplit3D split(eb, linearFactory(radius), level);
        auto bins = split.compress(conf, data.data());
        ASSERT_EQ(bins.size(), conf.num);

        std::vector<float> dec(conf.num);
        split.decompress(conf, bins, dec.data());
        EXPECT_LE(maxAbsErr(original, dec), eb) << "target_level=" << level;
    }
}

TEST(SZ3_MGARDComposable, RejectsQuantizerWithNonZeroRangeStart) {
    // ScalarQuantizer-style signed bins would break SZGenericCompressor's
    // "out range starts at 0" contract, so the composition refuses them up front.
    struct SignedQuantizer : SZ3::LinearQuantizer<float> {
        using SZ3::LinearQuantizer<float>::LinearQuantizer;
        std::pair<int, int> get_out_range() const override { return std::make_pair(-1, 1); }
    };
    using Split = SZ3::MultiLevelDecomposition<float, 3, SignedQuantizer>;
    EXPECT_THROW(Split(1e-2, [](double eb) { return SignedQuantizer(eb, 32768); }), std::invalid_argument);
}
