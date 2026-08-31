/**
 * Tests for GranularBitRoundQuantizer and ClusterQuantizer.
 *
 * Both are pure scalar quantizers (they ignore `pred`), so every test drives them
 * directly rather than through a compressor.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

#include "SZ3/quantizer/ClusterQuantizer.hpp"
#include "SZ3/quantizer/GranularBitRoundQuantizer.hpp"
#include "gtest/gtest.h"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Independent restatement of the guarantee: 0.5 * 10^(floor(log10|x|) + 1 - nsd).
/// Deliberately does NOT reuse the quantizer's own arithmetic.
static double granular_bound(double x, int nsd) {
    const double a = std::fabs(x);
    if (a == 0.0 || !std::isfinite(a)) return 0.0;
    const int d10 = static_cast<int>(std::floor(std::log10(a))) + 1;
    return 0.5 * std::pow(10.0, static_cast<double>(d10 - nsd));
}

/// Slack for the comparison itself: the bound is computed with a different rounding
/// path than the quantizer uses, and ties-to-even can land exactly on the bound.
static constexpr double kBoundSlack = 1.0 + 1e-9;

template <typename T>
static void expect_granular_guarantee(T original, T reconstructed, int nsd, const char *where) {
    if (!std::isfinite(static_cast<double>(original))) {
        // Non-finite input must be passed through bit-exactly.
        if (std::isnan(static_cast<double>(original))) {
            EXPECT_TRUE(std::isnan(static_cast<double>(reconstructed))) << where;
        } else {
            EXPECT_EQ(reconstructed, original) << where;
        }
        return;
    }
    if (!std::isnormal(original)) {
        // Zero and subnormals are passed through bit-exactly.
        EXPECT_EQ(reconstructed, original) << where;
        return;
    }
    const double err = std::fabs(static_cast<double>(reconstructed) - static_cast<double>(original));
    const double bound = granular_bound(static_cast<double>(original), nsd);
    EXPECT_LE(err, bound * kBoundSlack) << where << " value=" << static_cast<double>(original) << " nsd=" << nsd;
    // Relative form of the same guarantee.
    EXPECT_LE(err / std::fabs(static_cast<double>(original)), 0.5 * std::pow(10.0, 1.0 - nsd) * kBoundSlack) << where;
    // A normal value must stay normal and keep its sign (the product would underflow, so
    // compare sign bits instead).
    EXPECT_TRUE(std::isnormal(reconstructed)) << where;
    EXPECT_EQ(std::signbit(reconstructed), std::signbit(original)) << where;
}

// ---------------------------------------------------------------------------
// GranularBitRoundQuantizer
// ---------------------------------------------------------------------------

TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundOutRangeStartsAtZero) {
    SZ3::GranularBitRoundQuantizer<float> qf(4);
    SZ3::GranularBitRoundQuantizer<double> qd(4);
    EXPECT_EQ(qf.get_out_range().first, 0u);
    EXPECT_EQ(qd.get_out_range().first, 0u);
    // No usable bin range, so the encoder derives one from the bins.
    EXPECT_EQ(qf.get_out_range().second, 0u);
    EXPECT_EQ(qd.get_out_range().second, 0u);
    // A signed bin type would go negative here; check the sign-bit-set case explicitly.
    float negf = -1.5f;
    double negd = -1.5;
    EXPECT_GE(qf.quantize_and_overwrite(negf, 0.f), qf.get_out_range().first);
    EXPECT_GE(qd.quantize_and_overwrite(negd, 0.0), qd.get_out_range().first);
}

TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundRejectsBadNsd) {
    EXPECT_THROW(SZ3::GranularBitRoundQuantizer<float>(0), std::invalid_argument);
    EXPECT_THROW(SZ3::GranularBitRoundQuantizer<float>(19), std::invalid_argument);
    EXPECT_NO_THROW(SZ3::GranularBitRoundQuantizer<float>(1));
    EXPECT_NO_THROW(SZ3::GranularBitRoundQuantizer<double>(18));
}

template <typename T>
static void runGranularSpecialValues(int nsd) {
    SZ3::GranularBitRoundQuantizer<T> q(nsd);
    const std::vector<T> vals = {
        static_cast<T>(0.0),
        static_cast<T>(-0.0),
        static_cast<T>(1.0),
        static_cast<T>(-1.0),
        static_cast<T>(0.1),
        static_cast<T>(-123.456),
        static_cast<T>(999.9999),
        static_cast<T>(1000.0),
        static_cast<T>(1e-30),
        static_cast<T>(-1e-30),
        static_cast<T>(1e30),
        static_cast<T>(-1e30),
        std::numeric_limits<T>::max(),
        std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::min(),         // smallest normal
        std::numeric_limits<T>::denorm_min(),  // subnormal
        std::numeric_limits<T>::infinity(),
        -std::numeric_limits<T>::infinity(),
        std::numeric_limits<T>::quiet_NaN(),
    };

    for (size_t i = 0; i < vals.size(); i++) {
        T data = vals[i];
        const uint64_t bin = q.quantize_and_overwrite(data, T(7));  // pred must be ignored
        const T rec = q.recover(T(-11), bin);                       // pred must be ignored
        EXPECT_TRUE(std::memcmp(&rec, &data, sizeof(T)) == 0) << "recover != overwrite at i=" << i;
        expect_granular_guarantee<T>(vals[i], data, nsd, "special values");

        // The bin is a rounded IEEE-754 pattern, so it must fit in T's width and stay
        // non-negative even for a negative-signed `double`.
        EXPECT_GE(bin, q.get_out_range().first) << "bin below range at i=" << i;
        if (sizeof(T) < sizeof(uint64_t)) {
            EXPECT_LT(bin, uint64_t{1} << (sizeof(T) * 8)) << "bin too wide at i=" << i;
        }
    }

    // Zero must keep its sign bit.
    T zero = static_cast<T>(-0.0);
    q.quantize_and_overwrite(zero, T(0));
    EXPECT_TRUE(std::signbit(zero));
}

TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundSpecialValues) {
    for (int nsd = 1; nsd <= 7; nsd++) {
        runGranularSpecialValues<float>(nsd);
        runGranularSpecialValues<double>(nsd);
    }
}

/// The core claim: over a wide random sample the stated bound is never violated.
template <typename T>
static void runGranularRandomizedGuarantee(int nsd, int min_exp10, int max_exp10) {
    SZ3::GranularBitRoundQuantizer<T> q(nsd);
    std::mt19937 gen(20240613u + nsd);
    std::uniform_real_distribution<double> mant(1.0, 10.0);
    std::uniform_int_distribution<int> ex(min_exp10, max_exp10);
    std::uniform_int_distribution<int> sign(0, 1);

    for (int i = 0; i < 20000; i++) {
        const double raw = mant(gen) * std::pow(10.0, ex(gen)) * (sign(gen) ? 1.0 : -1.0);
        const T original = static_cast<T>(raw);
        if (!std::isnormal(original)) continue;
        T data = original;
        const uint64_t bin = q.quantize_and_overwrite(data, T(3.5));
        expect_granular_guarantee<T>(original, data, nsd, "randomized");
        EXPECT_EQ(q.recover(T(0), bin), data) << "recover mismatch at i=" << i;

        // Granular BitRound is deliberately NOT idempotent (rounding can move a value into
        // a coarser decade), but a second pass must still honour the guarantee relative to
        // its own input.
        T second = data;
        q.quantize_and_overwrite(second, T(0));
        expect_granular_guarantee<T>(data, second, nsd, "second pass");

        // The retained-bit count must actually be granular, i.e. the dropped bits
        // really are zero in the reconstruction.
        const int keep = q.keep_mantissa_bits(original);
        const int drop = SZ3::GranularBitRoundQuantizer<T>::kMantBits - keep;
        if (drop > 0) {
            uint64_t bits = 0;
            std::memcpy(&bits, &data, sizeof(T));
            EXPECT_EQ(bits & ((uint64_t{1} << drop) - 1), uint64_t{0}) << "low bits not cleared at i=" << i;
        }
    }
}

TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundRandomizedGuarantee) {
    for (int nsd = 1; nsd <= 7; nsd++) {
        runGranularRandomizedGuarantee<float>(nsd, -30, 30);
        runGranularRandomizedGuarantee<double>(nsd, -280, 280);
    }
}

/// Values sitting right at decade and binade boundaries, where the derivation of the
/// retained-bit count is most fragile.
template <typename T>
static void runGranularBoundaryValues(int nsd) {
    SZ3::GranularBitRoundQuantizer<T> q(nsd);
    std::vector<T> vals;
    for (int e = -20; e <= 20; e++) {
        const double p = std::pow(10.0, e);
        for (double f : {1.0, 0.9999999, 1.0000001, 9.999999, 5.0, 2.0}) {
            const T v = static_cast<T>(f * p);
            if (std::isnormal(v)) vals.push_back(v);
        }
        const double b = std::pow(2.0, e);
        for (double f : {1.0, 0.99999994, 1.9999999}) {
            const T v = static_cast<T>(f * b);
            if (std::isnormal(v)) vals.push_back(v);
        }
    }
    for (T v : vals) {
        T data = v;
        q.quantize_and_overwrite(data, T(0));
        expect_granular_guarantee<T>(v, data, nsd, "boundary");
        T neg = -v;
        q.quantize_and_overwrite(neg, T(0));
        expect_granular_guarantee<T>(static_cast<T>(-v), neg, nsd, "boundary-negative");
    }
}

TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundBoundaryValues) {
    for (int nsd = 1; nsd <= 8; nsd++) {
        runGranularBoundaryValues<float>(nsd);
        runGranularBoundaryValues<double>(nsd);
    }
}

/// A big enough nsd cannot be met by dropping any bits, so the quantizer degrades to
/// a bit-exact pass-through: 8 digits for float, 17 for double.
TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundLosslessAtFullPrecision) {
    SZ3::GranularBitRoundQuantizer<float> qf(8);
    SZ3::GranularBitRoundQuantizer<double> qd(17);
    std::mt19937 gen(99u);
    std::uniform_real_distribution<double> mant(1.0, 10.0);
    std::uniform_int_distribution<int> ex(-20, 20);
    for (int i = 0; i < 4000; i++) {
        const double raw = mant(gen) * std::pow(10.0, ex(gen));
        const float f = static_cast<float>(raw);
        const double d = raw;
        float fdata = f;
        double ddata = d;
        qf.quantize_and_overwrite(fdata, 0.f);
        qd.quantize_and_overwrite(ddata, 0.0);
        EXPECT_EQ(fdata, f) << "float not lossless at nsd=8, i=" << i;
        EXPECT_EQ(ddata, d) << "double not lossless at nsd=17, i=" << i;
    }
}

/// Lower nsd must never keep more bits than higher nsd (monotonicity of the model).
TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundKeepBitsMonotone) {
    std::mt19937 gen(7u);
    std::uniform_real_distribution<double> mant(1.0, 10.0);
    std::uniform_int_distribution<int> ex(-25, 25);
    for (int i = 0; i < 2000; i++) {
        const float v = static_cast<float>(mant(gen) * std::pow(10.0, ex(gen)));
        if (!std::isnormal(v)) continue;
        int prev = -1;
        for (int nsd = 1; nsd <= 8; nsd++) {
            SZ3::GranularBitRoundQuantizer<float> q(nsd);
            const int keep = q.keep_mantissa_bits(v);
            EXPECT_GE(keep, prev) << "keep bits not monotone in nsd, v=" << v;
            prev = keep;
        }
    }
}

TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundSaveLoad) {
    constexpr int nsd = 4;
    constexpr int N = 256;
    SZ3::GranularBitRoundQuantizer<double> q(nsd);

    std::vector<double> originals(N);
    std::vector<double> reconstructed(N);
    std::vector<uint64_t> bins(N);
    std::mt19937 gen(4242u);
    std::uniform_real_distribution<double> mant(1.0, 10.0);
    std::uniform_int_distribution<int> ex(-12, 12);
    for (int i = 0; i < N; i++) {
        originals[i] = mant(gen) * std::pow(10.0, ex(gen)) * ((i % 3) ? 1.0 : -1.0);
        double data = originals[i];
        bins[i] = q.quantize_and_overwrite(data, 0.0);
        reconstructed[i] = data;
    }

    std::vector<unsigned char> buffer(64);
    unsigned char *save_ptr = buffer.data();
    q.save(save_ptr);
    const size_t saved_size = static_cast<size_t>(save_ptr - buffer.data());
    EXPECT_GT(saved_size, 0u);

    // A quantizer configured differently must be fully overwritten by load().
    SZ3::GranularBitRoundQuantizer<double> loaded(9);
    const unsigned char *load_ptr = buffer.data();
    size_t remaining = saved_size;
    loaded.load(load_ptr, remaining);
    EXPECT_EQ(remaining, 0u) << "load() must consume exactly what save() wrote";
    EXPECT_EQ(load_ptr, buffer.data() + saved_size);
    EXPECT_EQ(loaded.nsd(), nsd);

    for (int i = 0; i < N; i++) {
        EXPECT_EQ(loaded.recover(0.0, bins[i]), reconstructed[i]) << "i=" << i;
        // The reloaded quantizer must also make the same decisions on fresh data.
        double data = originals[i];
        EXPECT_EQ(loaded.quantize_and_overwrite(data, 0.0), bins[i]) << "i=" << i;
        expect_granular_guarantee<double>(originals[i], data, nsd, "after load");
    }
}

TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundLoadRejectsForeignStream) {
    SZ3::GranularBitRoundQuantizer<float> q(3);
    unsigned char buffer[16] = {0};
    buffer[0] = 0xAB;  // wrong uid
    const unsigned char *p = buffer;
    size_t remaining = sizeof(buffer);
    EXPECT_THROW(q.load(p, remaining), std::invalid_argument);
}

TEST(SZ3_GranularBitRoundQuantizerTest, GranularBitRoundForceSaveUnpredIsExact) {
    SZ3::GranularBitRoundQuantizer<float> q(2);
    const float vals[] = {3.14159265f, -2.71828f, 1e-20f, 1e20f};
    for (float v : vals) {
        const uint64_t bin = q.force_save_unpred(v);
        EXPECT_LT(bin, uint64_t{1} << (sizeof(float) * 8));
        EXPECT_EQ(q.recover(0.f, bin), v);
    }
}

// ---------------------------------------------------------------------------
// ClusterQuantizer
// ---------------------------------------------------------------------------

/// Build MD-like data: a uniform lattice of levels plus small noise.
static std::vector<float> makeClusteredData(size_t n, float start, float offset, int nlevels, float noise_sigma,
                                            uint32_t seed) {
    std::mt19937 gen(seed);
    std::normal_distribution<float> noise(0.f, noise_sigma);
    std::uniform_int_distribution<int> lvl(0, nlevels - 1);
    std::vector<float> data(n);
    for (auto &v : data) {
        v = start + offset * static_cast<float>(lvl(gen)) + noise(gen);
    }
    return data;
}

TEST(SZ3_ClusterQuantizerTest, ClusterOutRangeStartsAtZero) {
    SZ3::ClusterQuantizer<float> q(10.0f, 2.5f, 12);
    EXPECT_EQ(q.get_out_range().first, 0);
    EXPECT_EQ(q.get_out_range().second, 13);  // bin 0 (unpred) + 12 levels
}

TEST(SZ3_ClusterQuantizerTest, ClusterRejectsBadLevels) {
    EXPECT_THROW(SZ3::ClusterQuantizer<float>(0.f, 0.f, 4), std::invalid_argument);
    EXPECT_THROW(SZ3::ClusterQuantizer<float>(0.f, -1.f, 4), std::invalid_argument);
    EXPECT_THROW(SZ3::ClusterQuantizer<float>(0.f, 1.f, 0), std::invalid_argument);
    EXPECT_THROW(SZ3::ClusterQuantizer<float>(SZ3::ClusterLevels{}), std::invalid_argument);
    EXPECT_NO_THROW(SZ3::ClusterQuantizer<float>(0.f, 1.f, 4));
}

/// End-to-end: derive levels from a sample, quantize, verify the bound and the codebook hit rate.
TEST(SZ3_ClusterQuantizerTest, ClusterDeriveLevelsAndRoundTrip) {
    const float start = 10.0f, offset = 2.5f;
    const int nlevels = 12;
    auto data = makeClusteredData(4000, start, offset, nlevels, 0.05f, 1234u);

    const SZ3::ClusterLevels levels = SZ3::derive_cluster_levels(data.data(), data.size(), 2000);
    ASSERT_TRUE(levels.valid()) << "k-means failed to find the planted lattice";
    EXPECT_NEAR(levels.level_offset, offset, 0.05f);
    EXPECT_GE(levels.level_num, nlevels);

    SZ3::ClusterQuantizer<float> q(levels);
    const double eb = q.get_eb();
    EXPECT_NEAR(eb, 0.5 * levels.level_offset, 1e-6);

    std::vector<int> bins(data.size());
    std::vector<float> recon(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        float scratch = data[i];
        bins[i] = q.quantize_and_overwrite(scratch, 0.f);
        recon[i] = scratch;
        EXPECT_LE(std::fabs(static_cast<double>(scratch) - data[i]), eb) << "i=" << i;
        EXPECT_GE(bins[i], q.get_out_range().first);
        EXPECT_LT(bins[i], q.get_out_range().second);
    }
    // With noise this small essentially everything should land on the codebook.
    EXPECT_LT(q.get_unpred_num(), data.size() / 100) << "too many samples missed the codebook";

    SZ3::ClusterQuantizer<float> q2(levels);
    // Feed the same unpredictable values back for decoding.
    std::vector<unsigned char> buf(64 + q.get_unpred_num() * sizeof(float) + 64);
    unsigned char *sp = buf.data();
    q.save(sp);
    const unsigned char *lp = buf.data();
    size_t remaining = static_cast<size_t>(sp - buf.data());
    q2.load(lp, remaining);
    q2.predecompress_data();
    for (size_t i = 0; i < data.size(); i++) {
        EXPECT_EQ(q2.recover(0.f, bins[i]), recon[i]) << "i=" << i;
    }
}

TEST(SZ3_ClusterQuantizerTest, ClusterDeriveLevelsRejectsUnclusterableData) {
    std::mt19937 gen(5u);
    std::uniform_real_distribution<float> u(0.f, 1.f);
    std::vector<float> uniform(4000);
    for (auto &v : uniform) v = u(gen);
    EXPECT_FALSE(SZ3::derive_cluster_levels(uniform.data(), uniform.size(), 2000).valid());

    std::vector<float> constant(1000, 7.5f);
    EXPECT_FALSE(SZ3::derive_cluster_levels(constant.data(), constant.size()).valid());

    std::vector<float> tiny(3, 1.f);
    EXPECT_FALSE(SZ3::derive_cluster_levels(tiny.data(), tiny.size()).valid());
    EXPECT_FALSE(SZ3::derive_cluster_levels<float>(nullptr, 0).valid());
}

/// `derive_cluster_levels` must also work on double input (KmeansUtil itself only
/// instantiates for float, so the wrapper funnels the sample through float).
TEST(SZ3_ClusterQuantizerTest, ClusterDeriveLevelsDoubleInput) {
    auto f = makeClusteredData(3000, -3.0f, 1.75f, 6, 0.02f, 777u);
    std::vector<double> d(f.begin(), f.end());
    const SZ3::ClusterLevels levels = SZ3::derive_cluster_levels(d.data(), d.size(), 1500);
    ASSERT_TRUE(levels.valid());
    EXPECT_NEAR(levels.level_offset, 1.75f, 0.05f);

    SZ3::ClusterQuantizer<double> q(levels);
    for (double v : d) {
        double scratch = v;
        q.quantize_and_overwrite(scratch, 0.0);
        EXPECT_LE(std::fabs(scratch - v), q.get_eb());
    }
}

/// The derived lattice must span the full data range, not just the sampled subset,
/// and values outside that range must still round-trip exactly via the unpred path.
TEST(SZ3_ClusterQuantizerTest, ClusterOutsideSampledRange) {
    auto data = makeClusteredData(3000, 100.0f, 8.0f, 5, 0.02f, 31337u);
    const SZ3::ClusterLevels levels = SZ3::derive_cluster_levels(data.data(), data.size(), 1500);
    ASSERT_TRUE(levels.valid());

    SZ3::ClusterQuantizer<float> q(levels);
    const double eb = q.get_eb();

    // Values far outside [level(0), level(level_num-1)] -- including zero, negatives and
    // extreme magnitudes -- must be stored verbatim and recovered bit-exactly.
    const std::vector<float> outsiders = {0.f,
                                          -0.f,
                                          -1e6f,
                                          1e6f,
                                          -273.15f,
                                          std::numeric_limits<float>::max(),
                                          std::numeric_limits<float>::lowest(),
                                          1e-30f,
                                          std::numeric_limits<float>::infinity(),
                                          -std::numeric_limits<float>::infinity(),
                                          std::numeric_limits<float>::quiet_NaN()};
    std::vector<int> bins;
    for (float v : outsiders) {
        float scratch = v;
        const int bin = q.quantize_and_overwrite(scratch, 0.f);
        EXPECT_EQ(bin, 0) << "value " << v << " should not have hit the codebook";
        if (std::isnan(v)) {
            EXPECT_TRUE(std::isnan(scratch)) << "data must be left untouched on the unpred path";
        } else {
            EXPECT_EQ(scratch, v) << "data must be left untouched on the unpred path";
        }
        bins.push_back(bin);
    }
    EXPECT_EQ(q.get_unpred_num(), outsiders.size());

    q.predecompress_data();
    for (size_t i = 0; i < outsiders.size(); i++) {
        const float rec = q.recover(0.f, bins[i]);
        if (std::isnan(outsiders[i])) {
            EXPECT_TRUE(std::isnan(rec));
        } else {
            EXPECT_EQ(rec, outsiders[i]) << "i=" << i;
        }
    }

    // A value just past the top level but within eb of it is still allowed on the codebook
    // path -- the guarantee is about the error, not about being inside the range.
    const float top = levels.level_start + levels.level_offset * static_cast<float>(levels.level_num - 1);
    float near_top = top + static_cast<float>(eb) * 0.5f;
    const float near_top_ori = near_top;
    const int bin = q.quantize_and_overwrite(near_top, 0.f);
    EXPECT_GT(bin, 0);
    EXPECT_LE(std::fabs(static_cast<double>(near_top) - near_top_ori), eb);
}

TEST(SZ3_ClusterQuantizerTest, ClusterSaveLoad) {
    const float start = -3.0f, offset = 1.25f;
    const int nlevels = 8;
    SZ3::ClusterQuantizer<float> q(start, offset, nlevels, /*eb=*/0.2);

    std::mt19937 gen(2468u);
    std::uniform_int_distribution<int> lvl(0, nlevels - 1);
    std::normal_distribution<float> noise(0.f, 0.4f);  // large enough to force some unpreds
    const int N = 500;
    std::vector<float> originals(N), recon(N);
    std::vector<int> bins(N);
    for (int i = 0; i < N; i++) {
        originals[i] = start + offset * static_cast<float>(lvl(gen)) + noise(gen);
        float scratch = originals[i];
        bins[i] = q.quantize_and_overwrite(scratch, 0.f);
        recon[i] = scratch;
    }
    ASSERT_GT(q.get_unpred_num(), 0u) << "test needs some unpredictable values to be meaningful";
    ASSERT_LT(q.get_unpred_num(), static_cast<size_t>(N)) << "test needs some codebook hits too";

    std::vector<unsigned char> buffer(1024 + q.get_unpred_num() * sizeof(float));
    unsigned char *sp = buffer.data();
    q.save(sp);
    const size_t saved = static_cast<size_t>(sp - buffer.data());
    EXPECT_GT(saved, 0u);

    SZ3::ClusterQuantizer<float> q2;  // default-constructed: everything comes from the stream
    const unsigned char *lp = buffer.data();
    size_t remaining = saved;
    q2.load(lp, remaining);
    EXPECT_EQ(remaining, 0u) << "load() must consume exactly what save() wrote";
    EXPECT_EQ(q2.get_level_start(), start);
    EXPECT_EQ(q2.get_level_offset(), offset);
    EXPECT_EQ(q2.get_level_num(), nlevels);
    EXPECT_EQ(q2.get_eb(), 0.2);
    EXPECT_EQ(q2.get_unpred_num(), q.get_unpred_num());
    EXPECT_EQ(q2.get_out_range(), q.get_out_range());

    for (int i = 0; i < N; i++) {
        const float rec = q2.recover(0.f, bins[i]);
        EXPECT_EQ(rec, recon[i]) << "i=" << i;
        EXPECT_LE(std::fabs(static_cast<double>(rec) - originals[i]), q2.get_eb()) << "i=" << i;
    }
}

TEST(SZ3_ClusterQuantizerTest, ClusterLoadRejectsForeignStream) {
    SZ3::ClusterQuantizer<float> q;
    unsigned char buffer[64] = {0};
    buffer[0] = 0xAB;  // wrong uid
    const unsigned char *p = buffer;
    size_t remaining = sizeof(buffer);
    EXPECT_THROW(q.load(p, remaining), std::invalid_argument);
}

/// The core claim: whatever the data and whatever the codebook, the reconstruction error
/// never exceeds eb, and the unpredictable path is bit-exact.
TEST(SZ3_ClusterQuantizerTest, ClusterRandomizedGuarantee) {
    std::mt19937 gen(864213u);
    std::uniform_real_distribution<float> start_dist(-1000.f, 1000.f);
    std::uniform_real_distribution<float> offset_dist(1e-3f, 50.f);
    std::uniform_int_distribution<int> nlevel_dist(1, 64);
    std::uniform_int_distribution<int> ebmode(0, 2);

    for (int trial = 0; trial < 200; trial++) {
        const float start = start_dist(gen);
        const float offset = offset_dist(gen);
        const int nlevels = nlevel_dist(gen);
        const int mode = ebmode(gen);
        const double eb_arg = (mode == 0) ? 0.0 : (mode == 1 ? 0.1 * offset : 0.9 * offset);
        SZ3::ClusterQuantizer<float> q(start, offset, nlevels, eb_arg);
        const double eb = q.get_eb();
        EXPECT_GT(eb, 0.0);

        // Data spanning far more than the codebook range, plus specials.
        std::uniform_real_distribution<float> data_dist(start - 5.f * offset * nlevels, start + 5.f * offset * nlevels);
        std::vector<float> originals;
        std::vector<int> bins;
        std::vector<float> recon;
        for (int i = 0; i < 400; i++) {
            float v;
            switch (i % 40) {
                case 0:
                    v = 0.f;
                    break;
                case 1:
                    v = std::numeric_limits<float>::infinity();
                    break;
                case 2:
                    v = -std::numeric_limits<float>::infinity();
                    break;
                case 3:
                    v = std::numeric_limits<float>::quiet_NaN();
                    break;
                case 4:
                    v = std::numeric_limits<float>::max();
                    break;
                case 5:
                    v = std::numeric_limits<float>::lowest();
                    break;
                case 6:
                    v = 1e-38f;
                    break;
                default:
                    v = data_dist(gen);
                    break;
            }
            originals.push_back(v);
            float scratch = v;
            const int bin = q.quantize_and_overwrite(scratch, 12345.f);  // pred ignored
            bins.push_back(bin);
            recon.push_back(scratch);

            EXPECT_GE(bin, q.get_out_range().first) << "trial=" << trial << " i=" << i;
            EXPECT_LT(bin, q.get_out_range().second) << "trial=" << trial << " i=" << i;
            if (bin == 0) {
                // Exact path: value untouched.
                if (std::isnan(v)) {
                    EXPECT_TRUE(std::isnan(scratch));
                } else {
                    EXPECT_EQ(scratch, v) << "trial=" << trial << " i=" << i;
                }
            } else {
                ASSERT_FALSE(std::isnan(v)) << "NaN must never hit the codebook";
                EXPECT_LE(std::fabs(static_cast<double>(scratch) - static_cast<double>(v)), eb)
                    << "trial=" << trial << " i=" << i << " v=" << v;
            }
        }

        // Decode in encode order and confirm the reconstruction matches exactly.
        q.predecompress_data();
        for (size_t i = 0; i < bins.size(); i++) {
            const float rec = q.recover(0.f, bins[i]);
            if (std::isnan(recon[i])) {
                EXPECT_TRUE(std::isnan(rec)) << "trial=" << trial << " i=" << i;
            } else {
                EXPECT_EQ(rec, recon[i]) << "trial=" << trial << " i=" << i;
            }
        }
    }
}

TEST(SZ3_ClusterQuantizerTest, ClusterForceSaveUnpredIsExact) {
    SZ3::ClusterQuantizer<double> q(0.0f, 1.0f, 4);
    const double vals[] = {1.0 / 3.0, -1e300, 1e-300, 0.0};
    std::vector<int> bins;
    for (double v : vals) {
        const int bin = q.force_save_unpred(v);
        EXPECT_EQ(bin, 0);
        bins.push_back(bin);
    }
    q.predecompress_data();
    for (size_t i = 0; i < bins.size(); i++) {
        EXPECT_EQ(q.recover(0.0, bins[i]), vals[i]) << "i=" << i;
    }
}
