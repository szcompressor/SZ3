#include "gtest/gtest.h"
#include "SZ3/decomposition/SVDDecomposition.hpp"
#include "SZ3/utils/Config.hpp"
#include <vector>
#include <cmath>

// Test fixture for SVDDecomposition
class SVDDecompositionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up a 4x4 data block for testing
        data = {1, 2, 3, 4,
                5, 6, 7, 8,
                9, 10, 11, 12,
                13, 14, 15, 16};
        conf = SZ3::Config(4, 4);
        conf.cmprAlgo = SZ3::ALGO_SVD;
        conf.absErrorBound = 0.1; // Set a reasonable error bound
        conf.svd_target_rank = 2; // Use a rank of 2
    }

    std::vector<float> data;
    SZ3::Config conf;
};

TEST_F(SVDDecompositionTest, CompressDecompress) {
    // Instantiate SVDDecomposition
    auto decomposer = SZ3::SVDDecomposition<float, 2>(conf);

    // Perform compression
    auto quantized_diff = decomposer.compress(conf, data.data());

    // Perform decompression
    std::vector<float> reconstructed_data(data.size());
    decomposer.decompress(conf, quantized_diff, reconstructed_data.data());

    // Check the error
    float max_error = 0;
    for (size_t i = 0; i < data.size(); ++i) {
        float error = std::abs(data[i] - reconstructed_data[i]);
        if (error > max_error) {
            max_error = error;
        }
    }

    // The error should be within the absolute error bound
    EXPECT_LE(max_error, conf.absErrorBound * 2); // Allow for some quantization error
}

TEST_F(SVDDecompositionTest, SaveAndLoad) {
    auto decomposer1 = SZ3::SVDDecomposition<float, 2>(conf);
    std::vector<unsigned char> buffer(4096);
    auto decomposer2 = SZ3::SVDDecomposition<float, 2>(conf);

    // Perform compression with decomposer1
    auto quantized_diff1 = decomposer1.compress(conf, data.data());

    // Serialize decomposer1
    unsigned char* c = buffer.data();
    decomposer1.save(c);

    // Deserialize into decomposer2
    const unsigned char* c2 = buffer.data();
    size_t remaining_length = buffer.size();
    decomposer2.load(c2, remaining_length);

    // Decompress using decomposer2 and compare with original data
    std::vector<float> reconstructed_data(data.size());
    decomposer2.decompress(conf, quantized_diff1, reconstructed_data.data());

    // Check the error (same logic as CompressDecompress test)
    float max_error = 0;
    for (size_t i = 0; i < data.size(); ++i) {
        float error = std::abs(data[i] - reconstructed_data[i]);
        if (error > max_error) {
            max_error = error;
        }
    }

    // The error should be within the absolute error bound
    EXPECT_LE(max_error, conf.absErrorBound * 2); // Allow for some quantization error
}
