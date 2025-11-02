#include <cmath>
#include <cstdint>

#include "SZ3/encoder/ArithmeticEncoder.hpp"
#include "SZ3/encoder/BypassEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/RunlengthEncoder.hpp"
#include "gtest/gtest.h"

template <typename Encoder, typename T>
void runFunctionalTest() {
    int N = 1000;
    std::vector<SZ3::uchar> buffer_data(N * sizeof(T) * 2);
    std::vector<SZ3::uchar> buffer_conf(N * sizeof(T) * 2);
    std::vector<T> data(N);
    for (int i = 0; i < N; i++) {
        data[i] = i % 100;
    }

    size_t data_len = 0, conf_len = 0;
    {
        SZ3::uchar *buffer_data_pos = buffer_data.data();
        SZ3::uchar *buffer_conf_pos = buffer_conf.data();
        Encoder coder;
        coder.preprocess_encode(data, 100);
        coder.encode(data, buffer_data_pos);
        coder.save(buffer_conf_pos);
        data_len = buffer_data_pos - buffer_data.data();
        conf_len = buffer_conf_pos - buffer_conf.data();
    }
    {
        const SZ3::uchar *buffer_data_pos = buffer_data.data();
        const SZ3::uchar *buffer_conf_pos = buffer_conf.data();
        Encoder coder;
        coder.load(buffer_conf_pos, conf_len);
        auto dataDecoded = coder.decode(buffer_data_pos, N);
        for (int i = 0; i < N; i++) {
            EXPECT_EQ(data[i], dataDecoded[i]);
        }
    }
}

template <typename Encoder, typename T>
void runAllTest() {
    runFunctionalTest<Encoder, T>();
}

TEST(EncoderTest, HuffmanEncoder) { runAllTest<SZ3::HuffmanEncoder<int>, int>(); }

TEST(EncoderTest, RunlengthEncoder) { runAllTest<SZ3::RunlengthEncoder<int>, int>(); }

TEST(EncoderTest, ArithmeticEncoder) { runAllTest<SZ3::ArithmeticEncoder<int>, int>(); }

TEST(EncoderTest, BypassEncoder) { runAllTest<SZ3::BypassEncoder<int>, int>(); }
