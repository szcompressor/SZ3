/**
 * Tests for GranularBitRoundQuantizer: the "keep N significant decimal digits"
 * guarantee, its special-value handling, and its serialized state.
 *
 * It is a pure scalar quantizer (it ignores `pred`), so every test drives it
 * directly rather than through a compressor.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

#include "SZ3/quantizer/GranularBitRoundQuantizer.hpp"
#include "gtest/gtest.h"

static constexpr double kBoundSlack = 1.0 + 1e-9;

/// Independent restatement of the guarantee: 0.5 * 10^(floor(log10|x|) + 1 - nsd).
/// Deliberately does NOT reuse the quantizer's own arithmetic.
static double granular_bound(double x, int nsd) {
    const double a = std::fabs(x);
    if (a == 0.0 || !std::isfinite(a)) return 0.0;
    const int d10 = static_cast<int>(std::floor(std::log10(a))) + 1;
    return 0.5 * std::pow(10.0, static_cast<double>(d10 - nsd));
}

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
