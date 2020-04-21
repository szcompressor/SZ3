#include <quantizer/IntegerQuantizer.hpp>
#include <compressor/Compressor.hpp>
#include <quantizer/Quantizer.hpp>
#include <predictor/Predictor.hpp>
#include <utils/Iterator.hpp>
#include <utils/fileUtil.h>
#include <utils/Verification.hpp>
#include <def.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include "zstd.h"

unsigned long sz_lossless_compress(unsigned char *data, size_t dataLength) {
    unsigned long outSize = 0;
    size_t estimatedCompressedSize = 0;
    if (dataLength < 100)
        estimatedCompressedSize = 200;
    else
        estimatedCompressedSize = dataLength * 1.2;
    auto compressBytes = SZ::compat::make_unique<unsigned char[]>(estimatedCompressedSize);
    outSize = ZSTD_compress(compressBytes.get(), estimatedCompressedSize, data, dataLength,
                            3); //default setting of level is 3
    return outSize;
}


template<typename T, class Predictor, class... Args>
float
compress(std::unique_ptr<T[]> &data, size_t num, uint block_size, uint stride, Predictor predictor, bool lossless, T eb,
         Args... args
) {
    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + num);

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    auto sz = SZ::make_sz_general<T>(block_size, stride,
                                     predictor,
                                     SZ::LinearQuantizer<float>(eb),
                                     SZ::HuffmanEncoder<int>(),
                                     args...
    );

    size_t compressed_size = 0;
    auto compressed = sz.compress(data.get(), compressed_size);

    std::cout << "Compressed size before zstd = " << compressed_size << std::endl;
    if (lossless) {
        compressed_size = sz_lossless_compress(compressed, compressed_size * sizeof(float));
        std::cout << "Compressed size after zstd = " << compressed_size << std::endl;
    }

    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Compression time: "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << std::endl;

    auto ratio = num * sizeof(float) * 1.0 / compressed_size;
    std::cout << "********************Compression Ratio******************* = " << ratio << std::endl;

    clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<float[]> dec_data;
    dec_data.reset(sz.decompress(compressed, compressed_size));
    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Decompression time: "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << std::endl;

    SZ::verify<float>(data_.data(), data.get(), num);
    delete[] compressed;
    return ratio;
}

template<typename T, uint N, class... Args>
float
choose_compressor_and_compress(bool lorenzo_1, bool lorenzo_2, bool regression_1, bool regression_2, bool lossless,
                               std::unique_ptr<T[]> &data, size_t num,
                               uint block_size, uint stride, T eb, Args... args
) {
    std::vector<std::shared_ptr<SZ::VirtualPredictor<float, N>>> predictors;
    if (lorenzo_1) {
        auto P_l = std::make_shared<SZ::RealPredictor<float, N, SZ::LorenzoPredictor<float, N, 1>>>(
                std::make_shared<SZ::LorenzoPredictor<float, N, 1>>(eb));
        predictors.push_back(P_l);
    }
    if (lorenzo_2) {
        auto P_l2 = std::make_shared<SZ::RealPredictor<float, N, SZ::LorenzoPredictor<float, N, 2>>>(
                std::make_shared<SZ::LorenzoPredictor<float, N, 2>>(eb));
        predictors.push_back(P_l2);
    }
    if (regression_1) {
        auto P_r = std::make_shared<SZ::RealPredictor<float, N, SZ::RegressionPredictor<float, N>>>(
                std::make_shared<SZ::RegressionPredictor<float, N>>(block_size, eb));
        predictors.push_back(P_r);
    }
    if (regression_2) {
        auto P_r2 = std::make_shared<SZ::RealPredictor<float, N, SZ::PolyRegressionPredictor<float, N>>>(
                std::make_shared<SZ::PolyRegressionPredictor<float, N>>(block_size, eb));
        predictors.push_back(P_r2);
    }
    if (predictors.size() == 1) {
        return compress<T>(data, num, block_size, stride, predictors[0], lossless, eb, args...);
    } else {
        auto P_composed = std::make_shared<SZ::ComposedPredictor<T, N>>(predictors);
        return compress<T>(data, num, block_size, stride, P_composed, lossless, eb, args...);
    }
}

template<typename T>
float
sz(bool lorenzo_1, bool lorenzo_2, bool regression_1, bool regression_2, bool lossless, uint block_size, uint stride,
   uint pred_dim, T eb, std::unique_ptr<T[]> &data, size_t num, uint r1, uint r2, uint r3) {
    std::cout << "Options: "
              << "eb = " << eb
              << ", block_size = " << block_size
              << ", stride = " << stride
              << ", dim = " << pred_dim
              << ", lorenzo_1 = " << lorenzo_1
              << ", lorenzo_2 = " << lorenzo_2
              << ", regression_1 = " << regression_1
              << ", regression_2 = " << regression_2
              << ", lossless= " << lossless
              << std::endl;

    if (pred_dim == 3) {
        return choose_compressor_and_compress<float, 3>(lorenzo_1, lorenzo_2, regression_1, regression_2, lossless,
                                                        data, num, block_size, stride, eb, r1, r2, r3);
    } else if (pred_dim == 2) {
        return choose_compressor_and_compress<float, 2>(lorenzo_1, lorenzo_2, regression_1, regression_2, lossless,
                                                        data, num, block_size, stride, eb, r1 * r2, r3);
    } else {
        return choose_compressor_and_compress<float, 1>(lorenzo_1, lorenzo_2, regression_1, regression_2, lossless,
                                                        data, num, block_size, stride, eb, r1 * r2 * r3);
    }

}

int main(int argc, char **argv) {
    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
    std::cout << "Read " << num << " elements\n";

    int r1 = atoi(argv[2]);
    int r2 = atoi(argv[3]);
    int r3 = atoi(argv[4]);
    float reb = atof(argv[5]);
    int block_size = atoi(argv[6]);
    int pred_dim = atoi(argv[7]);
    int lorenzo_op = atoi(argv[8]);
    int regression_op = atoi(argv[9]);

    int stride = block_size;

    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    float eb = reb * (max - min);
    bool lossless = true;
    bool lorenzo_1 = lorenzo_op == 1 || lorenzo_op == 3;
    bool lorenzo_2 = lorenzo_op == 2 || lorenzo_op == 3;
    bool regression_1 = regression_op == 1 || regression_op == 3;
    bool regression_2 = regression_op == 2 || regression_op == 3;

    std::cout << "value range = " << max - min << std::endl;
    std::cout << "abs error bound = " << eb << std::endl;

    auto ratio = sz(lorenzo_1, lorenzo_2, regression_1, regression_2, lossless, block_size, stride, pred_dim, eb, data, num, r1,
                    r2, r3);
//    std::cerr << ratio;

    return 0;
}
