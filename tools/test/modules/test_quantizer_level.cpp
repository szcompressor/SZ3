// Tests for LevelQuantizer -- a LUT of non-uniformly spaced reconstruction levels with an
// ABSOLUTE error bound. The level curve is a runtime constructor parameter, so every behavioural
// test below runs against both `LevelCurve::Quadratic` and `LevelCurve::Log`.
//
// This file consolidates what used to be split across QuadraticLevelQuantizer (tested via the
// generic harness in test_quantizer.cpp) and LogLevelQuantizer (test_quantizer_log.cpp).
//
// The load-bearing test is `BoundHoldsRandomized*`: it asserts the advertised bound over
// randomized residuals swept *well beyond* the level table's range. That is exactly the regime in
// which the old QuadraticLevelQuantizer silently clamped to its extreme level and violated its own
// bound by ~150000x eb; see `QuadraticRespectsBoundBeyondTableRange` for the pinned regression.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

#include "SZ3/quantizer/LevelQuantizer.hpp"
#include "gtest/gtest.h"

namespace {

using SZ3::LevelCurve;

/// Both curves, so every test can loop over them.
constexpr LevelCurve kCurves[] = {LevelCurve::Quadratic, LevelCurve::Log};

const char *curveName(LevelCurve c) { return SZ3::LevelQuantizer<float>::curve_name(c); }

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
// Shape and configuration
// ---------------------------------------------------------------------------

TEST(SZ3_LevelQuantizer, OutRangeStartsAtZero) {
    for (LevelCurve curve : kCurves) {
        SZ3::LevelQuantizer<float> qf(1e-3, 1024, curve);
        EXPECT_EQ(qf.get_out_range().first, 0) << curveName(curve);
        EXPECT_EQ(qf.get_out_range().second, 2 * 1024 + 1) << curveName(curve);

        SZ3::LevelQuantizer<double> qd(1e-3, 64, curve);
        EXPECT_EQ(qd.get_out_range().first, 0) << curveName(curve);
    }
    // The default-constructed quantizer must satisfy the contract too.
    SZ3::LevelQuantizer<double> def;
    EXPECT_EQ(def.get_out_range().first, 0);
    EXPECT_EQ(def.get_curve(), LevelCurve::Quadratic) << "default curve is quadratic";
}

TEST(SZ3_LevelQuantizer, ConstructorStoresCurve) {
    for (LevelCurve curve : kCurves) {
        SZ3::LevelQuantizer<float> q(1e-2, 512, curve);
        EXPECT_EQ(q.get_curve(), curve);
        EXPECT_DOUBLE_EQ(q.get_eb(), 1e-2);
        EXPECT_EQ(q.get_radius(), 512);
    }
}

TEST(SZ3_LevelQuantizer, TableShape) {
    const double eb = 1e-3;
    const int radius = 4096;

    for (LevelCurve curve : kCurves) {
        SZ3::LevelQuantizer<double> q(eb, radius, curve);

        // Documented weakness: the smallest positive level is exactly eb, so there is no zero
        // level and a residual of exactly 0 reconstructs to -/+eb -- the worst case the bound
        // allows. Pinned here so the behaviour cannot drift silently.
        double data = 0.0;
        const int bin = q.quantize_and_overwrite(data, 0.0);
        EXPECT_NE(bin, 0) << curveName(curve);
        EXPECT_DOUBLE_EQ(std::fabs(data), eb) << curveName(curve);

        // Post-clamping keeps the top level strictly below the unclamped nominal top eb*(2R-1),
        // which both curves share.
        EXPECT_GT(q.max_level(), eb) << curveName(curve);
        EXPECT_LT(q.max_level(), eb * (2 * radius - 1)) << curveName(curve);
    }
}

// The two curves must really be different tables, not one curve wearing two names.
TEST(SZ3_LevelQuantizer, CurvesAreDistinct) {
    const double eb = 1e-3;
    const int radius = 4096;
    SZ3::LevelQuantizer<double> quad(eb, radius, LevelCurve::Quadratic);
    SZ3::LevelQuantizer<double> logq(eb, radius, LevelCurve::Log);

    // The log curve grows geometrically and so hits the 2*eb clamp earlier, spending more of its
    // level budget on the uniform tail: its top level ends up lower than the quadratic curve's.
    EXPECT_LT(logq.max_level(), quad.max_level());

    // And they assign different bins to the same residual somewhere in the mid range.
    int differing = 0;
    for (int i = 1; i <= 200; i++) {
        const double residual = 0.25 * i;
        double a = residual;
        double b = residual;
        if (quad.quantize_and_overwrite(a, 0.0) != logq.quantize_and_overwrite(b, 0.0)) {
            differing++;
        }
    }
    EXPECT_GT(differing, 0) << "the two curves produced identical bins for every probe";
}

TEST(SZ3_LevelQuantizer, ConstructorRejectsBadArgs) {
    for (LevelCurve curve : kCurves) {
        EXPECT_THROW(SZ3::LevelQuantizer<float>(0.0, 32768, curve), std::invalid_argument);
        EXPECT_THROW(SZ3::LevelQuantizer<float>(-1.0, 32768, curve), std::invalid_argument);
        EXPECT_THROW(SZ3::LevelQuantizer<float>(1e-3, 1, curve), std::invalid_argument);
        EXPECT_NO_THROW(SZ3::LevelQuantizer<float>(1e-3, 2, curve));
    }
    // An out-of-range curve value must be refused rather than quietly falling back to a curve.
    EXPECT_THROW(SZ3::LevelQuantizer<float>(1e-3, 64, static_cast<LevelCurve>(0x7f)), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Unpredictable path
// ---------------------------------------------------------------------------

TEST(SZ3_LevelQuantizer, UnpredictablePath) {
    const double eb = 1e-3;
    for (LevelCurve curve : kCurves) {
        SZ3::LevelQuantizer<double> q(eb, 512, curve);

        // A residual far beyond the table range cannot meet the bound -> stored verbatim.
        const double huge = 1e12;
        double data = huge;
        EXPECT_EQ(q.quantize_and_overwrite(data, 0.0), 0) << curveName(curve);
        EXPECT_EQ(data, huge) << "unpredictable values must not be overwritten, " << curveName(curve);
        EXPECT_EQ(q.unpred_count(), 1u) << curveName(curve);

        // Non-finite inputs must also take the unpredictable path rather than reaching the search.
        double nan_data = std::numeric_limits<double>::quiet_NaN();
        EXPECT_EQ(q.quantize_and_overwrite(nan_data, 0.0), 0) << curveName(curve);
        double inf_data = std::numeric_limits<double>::infinity();
        EXPECT_EQ(q.quantize_and_overwrite(inf_data, 0.0), 0) << curveName(curve);

        // force_save_unpred is the explicit route.
        EXPECT_EQ(q.force_save_unpred(-7.5), 0) << curveName(curve);
        EXPECT_EQ(q.unpred_count(), 4u) << curveName(curve);

        EXPECT_EQ(q.recover(0.0, 0), huge) << curveName(curve);
        EXPECT_TRUE(std::isnan(q.recover(0.0, 0))) << curveName(curve);
        EXPECT_TRUE(std::isinf(q.recover(0.0, 0))) << curveName(curve);
        EXPECT_EQ(q.recover(0.0, 0), -7.5) << curveName(curve);
    }
}

// Pinned regression for the bug that motivated the merge. The old QuadraticLevelQuantizer had no
// unpredictable path: a residual past the end of its table was clamped to the extreme level and
// returned as a valid bin. At eb=1e-3 and the default radius the table tops out near 49.15, so a
// residual of -200 came back with an error of ~150.8 -- about 150000x the advertised bound.
TEST(SZ3_LevelQuantizer, QuadraticRespectsBoundBeyondTableRange) {
    const double eb = 1e-3;
    SZ3::LevelQuantizer<double> q(eb, 32768, LevelCurve::Quadratic);
    ASSERT_LT(q.max_level(), 200.0) << "test premise: -200 must be outside the table";

    const double original = -200.0;
    double data = original;
    const int bin = q.quantize_and_overwrite(data, 0.0);
    EXPECT_EQ(bin, 0) << "a residual past the table must go to the unpredictable path";
    EXPECT_EQ(data, original) << "unpredictable values must not be overwritten";
    // recover() consumes one entry of the unpredictable buffer per call, so call it exactly once.
    const double recovered = q.recover(0.0, bin);
    EXPECT_EQ(recovered, original) << "verbatim storage must reconstruct exactly";
    EXPECT_LE(std::fabs(recovered - original), eb) << "the old quadratic module reported ~150848x eb here";
}

// ---------------------------------------------------------------------------
// The important one: the advertised ABSOLUTE bound must hold for every input, whether it is
// representable by a table level or falls back to verbatim storage. Residuals are deliberately
// swept far past the end of the level table.
// ---------------------------------------------------------------------------

template <typename T>
static void runBoundHoldsRandomized(double eb, int radius, LevelCurve curve) {
    SZ3::LevelQuantizer<T> q(eb, radius, curve);
    const double table_top = q.max_level();

    std::mt19937 rng(20260830u);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    std::uniform_real_distribution<double> value(-10.0, 10.0);
    std::uniform_real_distribution<double> tiny(-1e-4, 1e-4);
    std::uniform_real_distribution<double> big(-1e6, 1e6);

    constexpr int N = 24000;
    std::vector<T> originals(N);
    std::vector<T> preds(N);
    std::vector<int> bins(N);
    std::vector<T> in_place(N);

    int beyond_table = 0;
    for (int i = 0; i < N; i++) {
        // Mix scales so residuals land inside the dense part of the table, right at its clamped
        // tail, just past the end, and vastly past the end.
        T v;
        T p;
        switch (i % 8) {
            case 0:
                v = static_cast<T>(value(rng));
                p = static_cast<T>(value(rng));
                break;
            case 1:
                p = static_cast<T>(value(rng));
                v = static_cast<T>(p + tiny(rng));  // tiny residual
                break;
            case 2:
                v = static_cast<T>(big(rng));
                p = static_cast<T>(big(rng));
                break;
            case 3:
                v = static_cast<T>(0);
                p = static_cast<T>(tiny(rng));  // residual ~0: the documented worst case
                break;
            case 4:
                // Straddle the top of the table: 0.5x .. 1.5x table_top.
                p = static_cast<T>(value(rng));
                v = static_cast<T>(p + table_top * (1.0 + 0.5 * unit(rng)));
                break;
            case 5:
                // Just past the table: 1x .. 10x table_top.
                p = static_cast<T>(value(rng));
                v = static_cast<T>(p + table_top * (5.5 + 4.5 * unit(rng)));
                break;
            case 6:
                // Far past the table.
                p = static_cast<T>(value(rng));
                v = static_cast<T>(p + table_top * 1e3 * unit(rng));
                break;
            default:
                // Absurdly past the table -- the exact regime the old quadratic module broke in.
                p = static_cast<T>(value(rng));
                v = static_cast<T>(p + table_top * 1e6 * unit(rng));
                break;
        }
        originals[i] = v;
        preds[i] = p;
        if (std::fabs(static_cast<double>(v) - static_cast<double>(p)) > table_top) {
            beyond_table++;
        }
        T data = v;
        bins[i] = q.quantize_and_overwrite(data, p);
        in_place[i] = data;

        EXPECT_GE(bins[i], q.get_out_range().first) << "i=" << i;
        EXPECT_LT(bins[i], q.get_out_range().second) << "i=" << i;
    }
    // The sweep must actually have left the table, or it is not testing what it claims to.
    EXPECT_GT(beyond_table, N / 4) << "sweep never went beyond the level table";

    // Exercise save/load in the same test so the bound is checked against a *reloaded* quantizer.
    const auto buffer = saveToBuffer(q, q.size_est());
    SZ3::LevelQuantizer<T> loaded(1.0, 8, curve);  // deliberately wrong eb/radius; load() overwrites
    loadFromBuffer(loaded, buffer);
    EXPECT_DOUBLE_EQ(loaded.get_eb(), eb);
    EXPECT_EQ(loaded.get_radius(), radius);
    EXPECT_EQ(loaded.get_curve(), curve);

    int predictable = 0;
    double worst_ratio = 0.0;
    for (int i = 0; i < N; i++) {
        const T recovered = loaded.recover(preds[i], bins[i]);
        // Decompression must reproduce exactly what compression wrote back in place.
        EXPECT_EQ(recovered, in_place[i]) << "compress/decompress mismatch at i=" << i;
        // The advertised bound.
        const double err = std::fabs(static_cast<double>(recovered) - static_cast<double>(originals[i]));
        EXPECT_LE(err, eb) << "absolute error bound violated at i=" << i << " (bin " << bins[i] << ", curve "
                           << curveName(curve) << ")";
        worst_ratio = std::max(worst_ratio, err / eb);
        if (bins[i] != 0) {
            predictable++;
        }
    }
    // Sanity: the bound must not be met by sending everything to the unpredictable path.
    EXPECT_GT(predictable, N / 8) << "suspiciously few values were actually quantized";

    std::printf("[bound] curve=%-9s T=%-6s eb=%-8.1E radius=%-7d worst |err|/eb = %.6f (%d/%d quantized)\n",
                curveName(curve), sizeof(T) == 4 ? "float" : "double", eb, radius, worst_ratio, predictable, N);
}

TEST(SZ3_LevelQuantizer, BoundHoldsRandomizedFloat) {
    for (LevelCurve curve : kCurves) {
        runBoundHoldsRandomized<float>(1e-3, 4096, curve);
        runBoundHoldsRandomized<float>(1e-2, 512, curve);
    }
}

TEST(SZ3_LevelQuantizer, BoundHoldsRandomizedDouble) {
    for (LevelCurve curve : kCurves) {
        runBoundHoldsRandomized<double>(1e-4, 8192, curve);
        runBoundHoldsRandomized<double>(1.0, 64, curve);
    }
}

// ---------------------------------------------------------------------------
// Serialization
// ---------------------------------------------------------------------------

TEST(SZ3_LevelQuantizer, SaveLoadRoundTripIncludesUnpred) {
    const double eb = 1e-3;
    for (LevelCurve curve : kCurves) {
        SZ3::LevelQuantizer<float> q(eb, 1024, curve);

        constexpr int N = 200;
        std::vector<float> originals(N);
        std::vector<int> bins(N);
        for (int i = 0; i < N; i++) {
            // Every 10th value is far out of range and must be stored verbatim.
            const float v = (i % 10 == 0) ? static_cast<float>(1e7 + i) : static_cast<float>(0.5 + 0.001 * i);
            originals[i] = v;
            float data = v;
            bins[i] = q.quantize_and_overwrite(data, 0.5f);
        }
        EXPECT_EQ(q.unpred_count(), static_cast<size_t>(N / 10)) << curveName(curve);

        const auto buffer = saveToBuffer(q, q.size_est());
        SZ3::LevelQuantizer<float> loaded(1e-9, 4, curve);
        loadFromBuffer(loaded, buffer);

        EXPECT_DOUBLE_EQ(loaded.get_eb(), eb) << curveName(curve);
        EXPECT_EQ(loaded.get_radius(), 1024) << curveName(curve);
        EXPECT_EQ(loaded.get_curve(), curve) << curveName(curve);
        EXPECT_EQ(loaded.get_out_range(), q.get_out_range()) << curveName(curve);
        EXPECT_EQ(loaded.unpred_count(), q.unpred_count()) << curveName(curve);
        EXPECT_DOUBLE_EQ(loaded.max_level(), q.max_level()) << "the reloaded table must be identical";

        for (int i = 0; i < N; i++) {
            const float recovered = loaded.recover(0.5f, bins[i]);
            if (bins[i] == 0) {
                EXPECT_EQ(recovered, originals[i]) << "verbatim value must round-trip exactly, i=" << i;
            } else {
                EXPECT_LE(std::fabs(static_cast<double>(recovered) - originals[i]), eb) << "i=" << i;
            }
        }
    }
}

// The curve is part of the stream's identity: loading a stream written with one curve into a
// quantizer configured for the other must be detected, not silently accepted with the wrong table.
TEST(SZ3_LevelQuantizer, CurveMismatchIsDetected) {
    for (LevelCurve written : kCurves) {
        const LevelCurve other = (written == LevelCurve::Quadratic) ? LevelCurve::Log : LevelCurve::Quadratic;

        SZ3::LevelQuantizer<float> q(1e-3, 256, written);
        const auto buffer = saveToBuffer(q, 0);

        SZ3::LevelQuantizer<float> wrong(1e-3, 256, other);
        const unsigned char *load_ptr = buffer.data();
        size_t remaining = buffer.size();
        EXPECT_THROW(wrong.load(load_ptr, remaining), std::invalid_argument)
            << "stream written as " << curveName(written) << " was accepted by a " << curveName(other) << " quantizer";

        // The matching curve must of course still load.
        SZ3::LevelQuantizer<float> right(1.0, 8, written);
        loadFromBuffer(right, buffer);
        EXPECT_EQ(right.get_curve(), written);
        EXPECT_DOUBLE_EQ(right.max_level(), q.max_level());
    }
}

TEST(SZ3_LevelQuantizer, RejectsForeignSerialization) {
    for (LevelCurve curve : kCurves) {
        SZ3::LevelQuantizer<float> q(1e-3, 64, curve);

        {  // Corrupt uid.
            auto buffer = saveToBuffer(q, 0);
            buffer[0] = 0x7f;
            SZ3::LevelQuantizer<float> loaded(1e-3, 64, curve);
            const unsigned char *load_ptr = buffer.data();
            size_t remaining = buffer.size();
            EXPECT_THROW(loaded.load(load_ptr, remaining), std::invalid_argument) << curveName(curve);
        }
        {  // Corrupt curve byte to a value no build knows.
            auto buffer = saveToBuffer(q, 0);
            buffer[1] = 0x5a;
            SZ3::LevelQuantizer<float> loaded(1e-3, 64, curve);
            const unsigned char *load_ptr = buffer.data();
            size_t remaining = buffer.size();
            EXPECT_THROW(loaded.load(load_ptr, remaining), std::invalid_argument) << curveName(curve);
        }
    }
}

}  // namespace
