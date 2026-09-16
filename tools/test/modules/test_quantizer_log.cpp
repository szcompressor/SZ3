// Tests for LogDomainQuantizer -- quantizes log|x| for a pointwise RELATIVE error bound.
//
// The log-spaced *level* LUT that used to be tested here (LogLevelQuantizer) was merged into
// LevelQuantizer; see test_quantizer_level.cpp.
//
// The load-bearing test is `BoundHoldsRandomized*`: it asserts the advertised bound over
// randomized inputs, including the inputs that must fall back to verbatim storage.

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "SZ3/quantizer/LogDomainQuantizer.hpp"
#include "gtest/gtest.h"

namespace {

template <typename Q>
std::vector<unsigned char> saveToBuffer(const Q &q, size_t extra_bytes) {
    std::vector<unsigned char> buffer(256 + extra_bytes);
    unsigned char *save_ptr = buffer.data();
    q.save(save_ptr);
    const size_t saved_size = static_cast<size_t>(save_ptr - buffer.data());
    EXPECT_GT(saved_size, 0u) << "Saved state must be non-empty.";
    EXPECT_LE(saved_size, buffer.size()) << "save() overran the buffer.";
    buffer.resize(saved_size);
    return buffer;
}

template <typename Q>
void loadFromBuffer(Q &q, const std::vector<unsigned char> &buffer) {
    const unsigned char *load_ptr = buffer.data();
    size_t remaining = buffer.size();
    q.load(load_ptr, remaining);
    EXPECT_EQ(remaining, 0u) << "load() must consume exactly what save() wrote.";
    EXPECT_EQ(static_cast<size_t>(load_ptr - buffer.data()), buffer.size()) << "load() must advance the pointer.";
}

// ---------------------------------------------------------------------------
// LogDomainQuantizer
// ---------------------------------------------------------------------------

TEST(SZ3_LogDomainQuantizer, OutRangeStartsAtZero) {
    SZ3::LogDomainQuantizer<float> q(1e-3, 1e-12, 1024);
    EXPECT_EQ(q.get_out_range().first, 0);
    EXPECT_EQ(q.get_out_range().second, 2 * q.get_radius() + 2);

    SZ3::LogDomainQuantizer<double> qd;
    EXPECT_EQ(qd.get_out_range().first, 0);
}

TEST(SZ3_LogDomainQuantizer, ConstructorRejectsBadArgs) {
    // A relative bound finer than the type's own resolution cannot be honored.
    EXPECT_THROW(SZ3::LogDomainQuantizer<float>(1e-9), std::invalid_argument);
    EXPECT_THROW(SZ3::LogDomainQuantizer<double>(0.0), std::invalid_argument);
    EXPECT_THROW(SZ3::LogDomainQuantizer<double>(1.5), std::invalid_argument);
    EXPECT_THROW(SZ3::LogDomainQuantizer<double>(1e-3, 0.0), std::invalid_argument);
    EXPECT_THROW(SZ3::LogDomainQuantizer<double>(1e-3, 1e-20, 0), std::invalid_argument);
    // A relative bound coarser than eps is fine for double.
    EXPECT_NO_THROW(SZ3::LogDomainQuantizer<double>(1e-9));
}

TEST(SZ3_LogDomainQuantizer, ZeroIsExact) {
    SZ3::LogDomainQuantizer<double> q(1e-3, 1e-12, 4096);

    double zero = 0.0;
    const int bin = q.quantize_and_overwrite(zero, 123.0);
    EXPECT_EQ(bin, 1);
    EXPECT_EQ(zero, 0.0);
    EXPECT_EQ(q.recover(123.0, bin), 0.0);
    EXPECT_EQ(q.unpred_count(), 0u) << "exact zero has its own bin, it is not unpredictable";

    // -0 is normalized to +0 (documented behavior).
    double neg_zero = -0.0;
    EXPECT_EQ(q.quantize_and_overwrite(neg_zero, 0.0), 1);
    EXPECT_EQ(neg_zero, 0.0);
    EXPECT_FALSE(std::signbit(neg_zero));
}

TEST(SZ3_LogDomainQuantizer, SignIsCarried) {
    const double rel_eb = 1e-3;
    SZ3::LogDomainQuantizer<double> q(rel_eb, 1e-12, 32768);

    for (double mag : {1e-8, 1e-3, 1.0, 7.25, 1234.5, 9.87e6}) {
        for (double sign : {1.0, -1.0}) {
            const double original = sign * mag;
            double data = original;
            const int bin = q.quantize_and_overwrite(data, 0.0);
            ASSERT_NE(bin, 0) << "mag=" << mag;
            EXPECT_EQ(std::signbit(data), std::signbit(original)) << "sign lost for " << original;
            EXPECT_LE(std::fabs(data - original) / std::fabs(original), rel_eb);
            EXPECT_EQ(q.recover(0.0, bin), data);
        }
    }
    EXPECT_EQ(q.unpred_count(), 0u);
}

TEST(SZ3_LogDomainQuantizer, OutOfWindowAndNonFiniteGoUnpredictable) {
    const double rel_eb = 1e-3;
    const double min_abs = 1e-6;
    SZ3::LogDomainQuantizer<double> q(rel_eb, min_abs, 1024);
    // Window is narrow on purpose so both edges are easy to cross.
    EXPECT_DOUBLE_EQ(q.min_magnitude(), min_abs);
    EXPECT_GT(q.max_magnitude(), min_abs);

    const double too_small = min_abs * 1e-3;
    double d1 = too_small;
    EXPECT_EQ(q.quantize_and_overwrite(d1, 0.0), 0);
    EXPECT_EQ(d1, too_small) << "unpredictable values must not be overwritten";

    const double too_big = q.max_magnitude() * 1e3;
    double d2 = too_big;
    EXPECT_EQ(q.quantize_and_overwrite(d2, 0.0), 0);
    EXPECT_EQ(d2, too_big);

    double nan_data = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(q.quantize_and_overwrite(nan_data, 0.0), 0);
    double inf_data = -std::numeric_limits<double>::infinity();
    EXPECT_EQ(q.quantize_and_overwrite(inf_data, 0.0), 0);

    EXPECT_EQ(q.force_save_unpred(42.0), 0);
    EXPECT_EQ(q.unpred_count(), 5u);

    EXPECT_EQ(q.recover(0.0, 0), too_small);
    EXPECT_EQ(q.recover(0.0, 0), too_big);
    EXPECT_TRUE(std::isnan(q.recover(0.0, 0)));
    EXPECT_TRUE(std::isinf(q.recover(0.0, 0)));
    EXPECT_EQ(q.recover(0.0, 0), 42.0);
}

TEST(SZ3_LogDomainQuantizer, Calibrate) {
    SZ3::LogDomainQuantizer<double> q(1e-3, 1e-20, 4096);
    const double max_abs = 5000.0;
    q.calibrate(max_abs);
    // The top rung lands at max_abs, so the whole window sits just under the data maximum.
    EXPECT_NEAR(q.max_magnitude(), max_abs, max_abs * 1e-9);
    EXPECT_LT(q.min_magnitude(), max_abs);

    double data = max_abs;
    EXPECT_NE(q.quantize_and_overwrite(data, 0.0), 0);
    EXPECT_LE(std::fabs(data - max_abs) / max_abs, 1e-3);

    EXPECT_THROW(q.calibrate(0.0), std::invalid_argument);
    EXPECT_THROW(q.calibrate(std::numeric_limits<double>::infinity()), std::invalid_argument);
}

// The important one: the advertised RELATIVE bound must hold for every input.
template <typename T>
static void runLogDomainBoundHoldsRandomized(double rel_eb, double min_abs, int radius, double log10_lo,
                                             double log10_hi) {
    SZ3::LogDomainQuantizer<T> q(rel_eb, min_abs, radius);

    std::mt19937 rng(20260830u);
    std::uniform_real_distribution<double> exponent(log10_lo, log10_hi);
    std::uniform_int_distribution<int> sign_bit(0, 1);
    std::uniform_int_distribution<int> kind(0, 9);

    constexpr int N = 20000;
    std::vector<T> originals(N);
    std::vector<int> bins(N);
    std::vector<T> in_place(N);

    for (int i = 0; i < N; i++) {
        T v;
        const int k = kind(rng);
        if (k == 0) {
            v = static_cast<T>(0);  // exact zeros
        } else if (k == 1) {
            v = static_cast<T>(min_abs * 1e-4);  // below the window
        } else {
            const double mag = std::pow(10.0, exponent(rng));
            v = static_cast<T>(sign_bit(rng) ? -mag : mag);
        }
        originals[i] = v;
        T data = v;
        // `pred` is ignored by this quantizer; feed junk to prove it.
        bins[i] = q.quantize_and_overwrite(data, static_cast<T>(-999.5));
        in_place[i] = data;

        EXPECT_GE(bins[i], q.get_out_range().first) << "i=" << i;
        EXPECT_LT(bins[i], q.get_out_range().second) << "i=" << i;
    }

    const auto buffer = saveToBuffer(q, q.size_est());
    SZ3::LogDomainQuantizer<T> loaded(0.5, 1.0, 2);  // deliberately wrong config
    loadFromBuffer(loaded, buffer);
    EXPECT_DOUBLE_EQ(loaded.get_eb(), rel_eb);
    EXPECT_EQ(loaded.get_radius(), q.get_radius());
    EXPECT_EQ(loaded.get_out_range(), q.get_out_range());

    int quantized = 0;
    for (int i = 0; i < N; i++) {
        const T recovered = loaded.recover(static_cast<T>(7.0), bins[i]);
        EXPECT_EQ(recovered, in_place[i]) << "compress/decompress mismatch at i=" << i;

        const double orig = static_cast<double>(originals[i]);
        const double rec = static_cast<double>(recovered);
        if (bins[i] == 0) {
            EXPECT_EQ(recovered, originals[i]) << "verbatim value must round-trip exactly, i=" << i;
        } else if (bins[i] == 1) {
            EXPECT_EQ(orig, 0.0);
            EXPECT_EQ(rec, 0.0);
        } else {
            ASSERT_NE(orig, 0.0);
            EXPECT_EQ(std::signbit(rec), std::signbit(orig)) << "sign flipped at i=" << i;
            EXPECT_LE(std::fabs(rec - orig) / std::fabs(orig), rel_eb)
                << "relative error bound violated at i=" << i << " (orig=" << orig << ", rec=" << rec << ")";
            quantized++;
        }
    }
    EXPECT_GT(quantized, N / 2) << "suspiciously few values were actually quantized";
}

TEST(SZ3_LogDomainQuantizer, BoundHoldsRandomizedFloat) {
    // Window [1e-12, ~2.9e16] comfortably covers the sampled magnitudes.
    runLogDomainBoundHoldsRandomized<float>(1e-3, 1e-12, 32768, -10.0, 6.0);
    // Coarse bound: the requested radius is clamped down so the top rung stays finite in float.
    runLogDomainBoundHoldsRandomized<float>(1e-1, 1e-8, 32768, -7.0, 5.0);
    // A bound close to float's guard band. The log-domain step shrinks with rel_eb, so covering
    // six decades needs a correspondingly large level budget -- this is the module's central
    // trade-off between relative precision and alphabet size.
    runLogDomainBoundHoldsRandomized<float>(1e-5, 1e-6, 1048576, -5.0, 1.0);
}

TEST(SZ3_LogDomainQuantizer, BoundHoldsRandomizedDouble) {
    runLogDomainBoundHoldsRandomized<double>(1e-3, 1e-12, 32768, -10.0, 6.0);
    runLogDomainBoundHoldsRandomized<double>(1e-4, 1e-9, 262144, -8.0, 2.0);
}

// Even when the representable window is far too narrow for the data, the bound must hold:
// everything that cannot be represented is stored verbatim rather than approximated badly.
TEST(SZ3_LogDomainQuantizer, NarrowWindowStillHonorsBound) {
    const double rel_eb = 1e-4;
    SZ3::LogDomainQuantizer<double> q(rel_eb, 1.0, 4);  // window spans a factor of ~1.0006

    std::mt19937 rng(7u);
    std::uniform_real_distribution<double> exponent(-6.0, 6.0);

    constexpr int N = 5000;
    int unpred = 0;
    for (int i = 0; i < N; i++) {
        const double original = std::pow(10.0, exponent(rng));
        double data = original;
        const int bin = q.quantize_and_overwrite(data, 0.0);
        if (bin == 0) {
            unpred++;
            EXPECT_EQ(data, original);
        } else {
            EXPECT_LE(std::fabs(data - original) / original, rel_eb) << "i=" << i;
        }
    }
    EXPECT_GT(unpred, N / 2) << "a window this narrow must reject most values";
}

TEST(SZ3_LogDomainQuantizer, RejectsForeignSerialization) {
    SZ3::LogDomainQuantizer<double> q(1e-3, 1e-12, 64);
    auto buffer = saveToBuffer(q, 0);
    buffer[0] = 0x7f;  // corrupt the uid
    SZ3::LogDomainQuantizer<double> loaded;
    const unsigned char *load_ptr = buffer.data();
    size_t remaining = buffer.size();
    EXPECT_THROW(loaded.load(load_ptr, remaining), std::invalid_argument);
}

}  // namespace
