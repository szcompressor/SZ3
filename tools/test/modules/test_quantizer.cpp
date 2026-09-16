#include <bitset>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "SZ3/quantizer/BitTruncationQuantizer.hpp"
#include "SZ3/quantizer/ClusterQuantizer.hpp"
#include "SZ3/quantizer/FixedPointQuantizer.hpp"
#include "SZ3/quantizer/GranularBitRoundQuantizer.hpp"
#include "SZ3/quantizer/LevelQuantizer.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/quantizer/LogDomainQuantizer.hpp"
#include "SZ3/quantizer/OutlierQuantizer.hpp"
#include "SZ3/quantizer/ScalarQuantizer.hpp"
#include "gtest/gtest.h"

template <typename Quantizer, typename T>
void runQuantizeRecoverTest() {
    T eb = 12.1973;
    Quantizer quantizer(eb);
    T data = static_cast<T>(100.11212);
    T data_ori = data;

    int quant_index = quantizer.quantize_and_overwrite(data, static_cast<T>(0));

    T recovered = quantizer.recover(static_cast<T>(0), quant_index);

    EXPECT_NEAR(recovered, data_ori, eb);
}

template <typename Quantizer, typename T>
void runFunctionalTest() {
    T eb = 12.1973;
    Quantizer quantizer(eb);

    const int N = 100;
    std::vector<T> originals(N);
    std::vector<int> quant_indices(N);

    // Generate 100 different input values.
    // We choose values that vary slightly so that quantization succeeds.
    for (int i = 0; i < N; i++) {
        T data = static_cast<T>(10.723 + (i * 0.0005));
        originals[i] = data;
        quant_indices[i] = quantizer.quantize_and_overwrite(data, static_cast<T>(0));
    }

    // Save the quantizer's state to a temporary buffer.
    std::vector<unsigned char> buffer(N * sizeof(T) * 4);
    unsigned char* save_ptr = buffer.data();
    quantizer.save(save_ptr);
    size_t saved_size = save_ptr - buffer.data();
    EXPECT_GT(saved_size, 0u) << "Saved state must be non-empty.";

    // Create a new quantizer instance and load the saved state.
    Quantizer loadedQuantizer(eb);  // The error bound will be overwritten by load().
    const unsigned char* load_ptr = buffer.data();
    size_t remaining_size = saved_size;
    loadedQuantizer.load(load_ptr, remaining_size);

    // Now, recover each value using the stored quantization indices.
    // Verify that the recovered value is within the error bound of the original.
    for (int i = 0; i < N; i++) {
        T recovered = loadedQuantizer.recover(static_cast<T>(0), quant_indices[i]);
        EXPECT_NEAR(recovered, originals[i], eb) << "Mismatch at index " << i;
    }
}

template <typename Quantizer, typename T>
void runAllTest() {
    runQuantizeRecoverTest<Quantizer, T>();
    runFunctionalTest<Quantizer, T>();
}

TEST(SZ3_QuantizerTest, LinearQuantizer) { runAllTest<SZ3::LinearQuantizer<float>, float>(); }

// Smoke-level pass through the shared harness (default curve). The thorough per-curve coverage,
// including the error-bound sweep, lives in test_quantizer_level.cpp.
TEST(SZ3_QuantizerTest, LevelQuantizer) { runAllTest<SZ3::LevelQuantizer<float>, float>(); }

// FixedPointQuantizer does not match the generic eb-only ctor harness:
// it is constructed with `num_bits`, then must be `calibrate(max_abs)`-d
// before use. Custom test below.
template <typename T>
static void runFixedPointTest() {
    constexpr int num_bits = 24;
    constexpr int N = 100;
    std::vector<T> originals(N);
    std::vector<int64_t> bins(N);

    SZ3::FixedPointQuantizer<T> q(num_bits);
    T max_abs = 0;
    for (int i = 0; i < N; i++) {
        T v = static_cast<T>(10.0 + 0.0005 * i);  // ascending positives
        originals[i] = v;
        if (std::fabs(v) > max_abs) max_abs = std::fabs(v);
    }
    q.calibrate(static_cast<double>(max_abs));
    const double tol = q.max_abs_error();

    for (int i = 0; i < N; i++) {
        T data = originals[i];
        bins[i] = q.quantize_and_overwrite(data, T(0));
    }

    std::vector<unsigned char> buf(64);
    unsigned char* sp = buf.data();
    q.save(sp);
    size_t saved = sp - buf.data();
    EXPECT_GT(saved, 0u);

    SZ3::FixedPointQuantizer<T> q2;
    const unsigned char* lp = buf.data();
    size_t remaining = saved;
    q2.load(lp, remaining);

    for (int i = 0; i < N; i++) {
        T recovered = q2.recover(T(0), bins[i]);
        EXPECT_NEAR(recovered, originals[i], tol) << "i=" << i;
    }
}

TEST(SZ3_QuantizerTest, FixedPointQuantizer) {
    runFixedPointTest<float>();
    runFixedPointTest<double>();
}

// ScalarQuantizer is constructed with (step, one_bin_reconstruct, tail_offset)
// rather than the eb-only ctor expected by `runAllTest`. Custom test below.
template <typename Ti>
static void runScalarQuantizerTest() {
    constexpr int N = 100;
    const double step = 0.5;
    SZ3::ScalarQuantizer<Ti, int64_t> q(step, /*one_bin_reconstruct=*/1.0, /*tail_offset=*/0.0);

    std::vector<Ti> originals(N);
    std::vector<int64_t> bins(N);
    for (int i = 0; i < N; i++) {
        // Use predictions of 0 and slowly varying data values; ScalarQuantizer's
        // rule means recovery error is bounded by step.
        originals[i] = static_cast<Ti>(0.123 * i - 5.0);
        Ti scratch = originals[i];
        bins[i] = q.quantize_and_overwrite(scratch, Ti(0));
        // After overwrite, scratch == reconstructed value.
        EXPECT_NEAR(scratch, originals[i], static_cast<Ti>(step) * 1.001) << "i=" << i;
    }

    std::vector<unsigned char> buf(64);
    unsigned char* sp = buf.data();
    q.save(sp);
    size_t saved = sp - buf.data();
    EXPECT_GT(saved, 0u);

    // Reload + recover -- reconstruct should match in-place quantize result.
    SZ3::ScalarQuantizer<Ti, int64_t> q2;
    const unsigned char* lp = buf.data();
    size_t remaining = saved;
    q2.load(lp, remaining);

    for (int i = 0; i < N; i++) {
        Ti recovered = q2.recover(Ti(0), bins[i]);
        EXPECT_NEAR(recovered, originals[i], static_cast<Ti>(step) * 1.001) << "i=" << i;
    }
}

TEST(SZ3_QuantizerTest, ScalarQuantizer) {
    runScalarQuantizerTest<float>();
    runScalarQuantizerTest<double>();
}

// BitTruncationQuantizer: AND'ing the float bit pattern with a top-bytes mask
// must (1) be idempotent on already-truncated values, (2) round-trip exactly
// through save/load, and (3) zero exactly the dropped low bytes.
template <typename T>
static void runBitTruncationTest(int keep_bytes) {
    constexpr int N = 64;
    SZ3::BitTruncationQuantizer<T> q(keep_bytes);

    std::vector<T> originals(N);
    std::vector<uint64_t> bins(N);
    for (int i = 0; i < N; i++) {
        originals[i] = static_cast<T>((i - 32) * 1.234567 + 0.0078125);  // mix of magnitudes & signs
    }

    for (int i = 0; i < N; i++) {
        T data = originals[i];
        bins[i] = q.quantize_and_overwrite(data, T(0));
        // Idempotence: re-truncating the reconstructed value yields the same bin.
        T data_again = data;
        int64_t bin_again = q.quantize_and_overwrite(data_again, T(0));
        EXPECT_EQ(bin_again, bins[i]) << "i=" << i << " not idempotent";
        EXPECT_EQ(data_again, data) << "i=" << i << " reconstruction drift";
    }

    std::vector<unsigned char> buf(64);
    unsigned char* sp = buf.data();
    q.save(sp);
    size_t saved = sp - buf.data();
    EXPECT_GT(saved, 0u);

    SZ3::BitTruncationQuantizer<T> q2(1);  // wrong initial keep_bytes; load() must overwrite
    const unsigned char* lp = buf.data();
    size_t remaining = saved;
    q2.load(lp, remaining);
    EXPECT_EQ(q2.keep_bytes(), keep_bytes);

    for (int i = 0; i < N; i++) {
        T recovered = q2.recover(T(0), bins[i]);
        // Must match what quantize_and_overwrite wrote back.
        T data = originals[i];
        T expected = data;
        q.quantize_and_overwrite(expected, T(0));
        EXPECT_EQ(recovered, expected) << "i=" << i;
    }
}

TEST(SZ3_QuantizerTest, BitTruncationQuantizer) {
    runBitTruncationTest<float>(2);   // keep top 2 of 4 bytes
    runBitTruncationTest<float>(3);   // keep top 3 of 4 bytes (low byte zeroed)
    runBitTruncationTest<double>(4);  // keep top 4 of 8 bytes
    runBitTruncationTest<double>(6);  // keep top 6 of 8 bytes
}

// Every quantizer writes its uid as the first byte of save(); load() rejects a mismatch, so two
// quantizers sharing a uid would silently accept each other's stream.
TEST(SZ3_QuantizerTest, UidsAreDistinct) {
    std::vector<SZ3::uchar> buffer(4096);
    std::map<SZ3::uchar, std::string> seen;
    auto record = [&](const std::string& name, auto&& quantizer) {
        SZ3::uchar* c = buffer.data();
        quantizer.save(c);
        SZ3::uchar uid = buffer[0];
        auto it = seen.find(uid);
        EXPECT_EQ(it, seen.end()) << name << " and " << (it == seen.end() ? "" : it->second) << " share uid 0b"
                                  << std::bitset<8>(uid);
        seen[uid] = name;
    };

    const double eb = 1e-3;
    record("LinearQuantizer", SZ3::LinearQuantizer<float>(eb));
    record("FixedPointQuantizer", SZ3::FixedPointQuantizer<float>(16));
    record("BitTruncationQuantizer", SZ3::BitTruncationQuantizer<float>(2));
    record("LevelQuantizer", SZ3::LevelQuantizer<float>(eb));
    record("ClusterQuantizer", SZ3::ClusterQuantizer<float>(0.0f, 1.0f, 8, eb));
    record("LogDomainQuantizer", SZ3::LogDomainQuantizer<float>(eb));
    record("GranularBitRoundQuantizer", SZ3::GranularBitRoundQuantizer<float>(3));
    record("OutlierQuantizer", SZ3::OutlierQuantizer<float>(eb));

    EXPECT_EQ(seen.size(), 8u);
}
