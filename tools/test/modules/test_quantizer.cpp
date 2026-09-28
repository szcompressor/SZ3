#include <cmath>
#include <cstdint>
#include <limits>

#include "SZ3/quantizer/LinearQuantizer.hpp"
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

TEST(QuantizerTest, LinearQuantizer) {
    runAllTest<SZ3::LinearQuantizer<float>, float>();
}

// Data at one end of an integer type predicted from the other end: data - pred and the reconstruction are
// outside the type. A kept value must be within the bound and be what recover() gives back.
template <typename T>
void runIntegerExtremesTest() {
    const T lo = std::numeric_limits<T>::lowest(), hi = std::numeric_limits<T>::max();
    const T pairs[][2] = {{hi, lo}, {lo, hi}, {hi, static_cast<T>(hi - 1)}, {lo, static_cast<T>(lo + 1)},
                          {0, lo},  {0, hi}};
    for (double eb : {0.5, 1.0, 3.0}) {
        for (const auto& p : pairs) {
            SZ3::LinearQuantizer<T> quantizer(eb);
            T data = p[0];
            int q = quantizer.quantize_and_overwrite(data, p[1]);
            if (q == 0) {
                EXPECT_EQ(data, p[0]);
                continue;
            }
            EXPECT_LE(std::fabs(static_cast<long double>(data) - p[0]), eb);
            EXPECT_EQ(quantizer.recover(p[1], q), data);
        }
    }
}

TEST(QuantizerTest, LinearQuantizerIntegerExtremes) {
    runIntegerExtremesTest<int8_t>();
    runIntegerExtremesTest<uint8_t>();
    runIntegerExtremesTest<int16_t>();
    runIntegerExtremesTest<uint16_t>();
    runIntegerExtremesTest<int32_t>();
    runIntegerExtremesTest<int64_t>();
}
