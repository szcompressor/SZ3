#include <quantizer/IntegerQuantizer.hpp>
#include <compressor/SZ_General_Compressor.hpp>
#include <predictor/ComposedPredictor.hpp>
#include <predictor/LorenzoPredictor.hpp>
#include <predictor/RegressionPredictor.hpp>
#include <lossless/Lossless_zstd.hpp>
#include <utils/Iterator.hpp>
#include <utils/FileUtil.h>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <memory>
#include <type_traits>

int main(int argc, char **argv) {
    size_t num = 0;
    // use Hurricane for testing
    auto data = SZ::readfile<float>(argv[1], num);
    std::cout << "Read " << num << " elements\n";
    float eb = 10;

    auto P_l = std::make_shared<SZ::LorenzoPredictor<float, 3, 1>>(eb);
    auto P_reg = std::make_shared<SZ::RegressionPredictor<float, 3>>(6, 0.1 * eb);
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<float, 3>>> predictors_;
    predictors_.push_back(P_l);
    predictors_.push_back(P_reg);
//    auto cp = std::make_shared<SZ::ComposedPredictor<float, 3>>(predictors_);
    SZ::Config<float, 3> conf(eb, std::array<size_t, 3>{100, 500, 500});
    auto sz = SZ::SZ_General_Compressor<float, 3, SZ::ComposedPredictor<float, 3>, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            conf,
            SZ::ComposedPredictor<float, 3>(predictors_),
            SZ::LinearQuantizer<float>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd()
    );
    size_t compressed_size = 0;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<unsigned char[]> compressed;
    compressed.reset(sz.compress(data.get(), compressed_size));
    err = clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Compression time: "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;
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
    return 0;
}
