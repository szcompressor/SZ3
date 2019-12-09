#include <quantizer/IntegerQuantizer.hpp>
#include <compressor/Compressor.hpp>
#include <quantizer/Quantizer.hpp>
#include <utils/Iterator.hpp>
#include <utils/fileUtil.h>
#include <predictor/Predictor.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>

#include "zstd.h"

unsigned long sz_lossless_compress(unsigned char *data, unsigned long dataLength) {
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



template<typename Type, class Predictor>
void compress(std::unique_ptr<Type[]> &data, Predictor predictor, Type eb, uint r1, uint r2, uint r3, size_t num) {
    auto sz = SZ::make_sz_general<Type>(
            // SZ::LorenzoPredictor<float, 3>(),
            // SZ::RegressionPredictor<float, 3>(0.1*eb),
            predictor,
            SZ::LinearQuantizer<float>(eb),
            SZ::HuffmanEncoder<int>(),
            r1,
            r2,
            r3
    );

    //SZ::SZ_General_Compressor<float, 3> sz(lorenzo, linear_quantizer, 100, 500, 500);

    size_t compressed_size = 0;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<unsigned char[]> compressed;
    compressed.reset(sz.compress(data.get(), eb, compressed_size));
    err = clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Compression time: "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;

    auto lossless_size = sz_lossless_compress(compressed.get(), compressed_size * sizeof(float));
    std::cout << "***Compression Ratio*** = " << num * sizeof(float) * 1.0 / lossless_size << std::endl;
    SZ::writefile("test.dat", compressed.get(), compressed_size);


    err = clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<float[]> dec_data;
    dec_data.reset(sz.decompress(compressed.get(), compressed_size));
    err = clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Decompression time: "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << std::endl;
    float max_err = 0;
    for (int i = 0; i < num; i++) {
        max_err = std::max(max_err, std::abs(data[i] - dec_data[i]));
    }
    std::cout << "Max error = " << max_err << std::endl;

}


int main(int argc, char **argv) {
    //usage: sz_test2 ../data/288x115x69x69.f32.qmc 33120 69 69  1e-7 6 3 0 0
    size_t num = 0;
    // use Hurricane for testing
    auto data = SZ::readfile<float>(argv[1], num);
    std::cout << "Read " << num << " elements\n";

    int r1 = atoi(argv[2]);
    int r2 = atoi(argv[3]);
    int r3 = atoi(argv[4]);
    float aeb = atof(argv[5]);
    int block_size = atoi(argv[6]);
    int pred_dim = atoi(argv[7]);
    int all_lorenzo = atoi(argv[8]);
    int all_regression = atoi(argv[9]);
    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    std::cout << "Options: " << "block_size = " << block_size << ", pred_dim = " << pred_dim << ", all_lorenzo = "
              << all_lorenzo << ", all_regression= " << all_regression << std::endl;
//    std::cout << std::endl;
    std::cout << "value range = " << max - min << std::endl;
    std::cout << "precision = " << aeb * (max - min) << std::endl;

    float eb = aeb * (max - min);


    std::shared_ptr<SZ::VirtualPredictor<float, 3>> best_predictor;

    auto P_l = std::make_shared<SZ::RealPredictor<float, 3, SZ::LorenzoPredictor<float, 3, 1>>>(
            std::make_shared<SZ::LorenzoPredictor<float, 3, 1>>(eb));
    auto P_l2 = std::make_shared<SZ::RealPredictor<float, 3, SZ::LorenzoPredictor<float, 3, 2>>>(
            std::make_shared<SZ::LorenzoPredictor<float, 3, 2>>(eb));
    auto P_reg = std::make_shared<SZ::RealPredictor<float, 3, SZ::RegressionPredictor<float, 3>>>(
            std::make_shared<SZ::RegressionPredictor<float, 3>>(0.1 * eb));
    auto cp = std::make_shared<SZ::ComposedPredictor<float, 3>>(P_l, P_reg);

    if (all_regression == 1) {
        compress<float>(data, P_reg, eb, r1, r2, r3, num);
    } else {
        if (all_lorenzo == 0) {
            compress<float>(data, cp, eb, r1, r2, r3, num);
        } else if (all_lorenzo == 1) {
            compress<float>(data, P_l, eb, r1, r2, r3, num);
        } else if (all_lorenzo == 2) {
            compress<float>(data, P_l2, eb, r1, r2, r3, num);
        }
    }

    return 0;
}
