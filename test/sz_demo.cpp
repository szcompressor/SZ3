#include "compressor/SZGeneralCompressor.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "encoder/HuffmanEncoder.hpp"
#include "lossless/Lossless_zstd.hpp"
#include "utils/FileUtil.h"
#include "utils/Config.hpp"
#include "utils/Verification.hpp"
#include "utils/Timer.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>

std::string src_file_name;
float relative_error_bound = 0;

template<typename T, class Predictor, uint N>
float SZ_Compress(std::unique_ptr<T[]> const &data,
                  const SZ::Config<T, N> &conf,
                  Predictor predictor) {

    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", block_size = " << conf.block_size
              << ", stride = " << conf.stride
              << ", quan_bin = " << conf.quant_bin
              << std::endl
              << "lorenzo = " << conf.enable_lorenzo
              << ", 2ndlorenzo = " << conf.enable_2ndlorenzo
              << ", regression = " << conf.enable_regression
              << ", lossless = " << conf.enable_lossless
              << std::endl;

    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

    SZ::Timer timer;
    timer.start();
    std::cout << "****************** Compression ******************" << std::endl;

    auto sz = SZ::make_sz_general_compressor(conf, predictor, SZ::LinearQuantizer<T>(conf.eb, conf.quant_bin),
                                             SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());

    size_t compressed_size = 0;
    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz.compress(data.get(), compressed_size));

    timer.stop("Compression");

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;

    std::stringstream ss;
    ss << src_file_name.substr(src_file_name.rfind('/') + 1) << "." << relative_error_bound << ".sz3";
    auto compressed_file_name = ss.str();
    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    std::cout << "Compressed file = " << compressed_file_name << std::endl;

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    std::unique_ptr<T[]> dec_data;
    dec_data.reset(sz.decompress(compressed.get(), compressed_size));
    timer.stop("Decompression");

    SZ::verify<T>(data_.data(), dec_data.get(), conf.num);

    auto decompressed_file_name = compressed_file_name + ".out";
    SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

    return ratio;
}

template<typename T, uint N>
float SZ_Compress(std::unique_ptr<T[]> const &data, const SZ::Config<T, N> &conf) {


    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
    int use_single_predictor =false;
//            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return SZ_Compress<T>(data, conf, SZ::LorenzoPredictor<T, N, 1>(conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return SZ_Compress<T>(data, conf, SZ::LorenzoPredictor<T, N, 2>(conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.eb));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return SZ_Compress<T>(data, conf, SZ::RegressionPredictor<T, N>(conf.block_size, conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
        }
    }

    return SZ_Compress<T>(data, conf, SZ::ComposedPredictor<T, N>(predictors));
}

template<typename T>
float SZ_Compress(std::unique_ptr<T[]> const &data, T eb, size_t r1, size_t r2, size_t r3, size_t r4) {
    return SZ_Compress(data, SZ::Config<T, 4>(eb, std::array<size_t, 4>{r1, r2, r3, r4}));
}

template<typename T>
float SZ_Compress(std::unique_ptr<T[]> const &data, T eb, size_t r1, size_t r2, size_t r3) {
    return SZ_Compress(data, SZ::Config<T, 3>(eb, std::array<size_t, 3>{r1, r2, r3}));
}

template<typename T>
float SZ_Compress(std::unique_ptr<T[]> const &data, T eb, size_t r1, size_t r2) {
    return SZ_Compress(data, SZ::Config<T, 2>(eb, std::array<size_t, 2>{r1, r2}));
}

template<typename T>
float SZ_Compress(std::unique_ptr<T[]> const &data, T eb, size_t r1) {
    return SZ_Compress(data, SZ::Config<T, 1>(eb, std::array<size_t, 1>{r1}));
}

template<class T, uint N>
float SZ_Compress_by_config(int argc, char **argv, int argp, std::unique_ptr<T[]> &data, float eb,
                            std::array<size_t, N> dims) {
    SZ::Config<float, N> conf(eb, dims);
    if (argp < argc) {
        int block_size = atoi(argv[argp++]);
        conf.block_size = block_size;
        conf.stride = block_size;
    }

    if (argp + 1 < argc) {
        int lorenzo_op = atoi(argv[argp++]);
        int regression_op = atoi(argv[argp++]);
        conf.enable_lorenzo = lorenzo_op == 1 || lorenzo_op == 3;
        conf.enable_2ndlorenzo = lorenzo_op == 2 || lorenzo_op == 3;
        conf.enable_regression = regression_op == 1 || regression_op == 3;
    }

    if (argp < argc) {
        conf.quant_bin = atoi(argv[argp++]);
    }

    return SZ_Compress(data, conf);
}

int main(int argc, char **argv) {
    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
    src_file_name = argv[1];
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

//    char *eb_op = argv[argp++] + 1;
    float eb = 0;
//    if (*eb_op == 'a') {
//        eb = atof(argv[argp++]);
//    } else {
    relative_error_bound = atof(argv[argp++]);
    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    eb = relative_error_bound * (max - min);
//    }

    if (argp == argc) {
        if (dim == 1) {
            SZ_Compress(data, eb, dims[0]);
        } else if (dim == 2) {
            SZ_Compress(data, eb, dims[0], dims[1]);
        } else if (dim == 3) {
            SZ_Compress(data, eb, dims[0], dims[1], dims[2]);
        } else if (dim == 4) {
            SZ_Compress(data, eb, dims[0], dims[1], dims[2], dims[3]);
        }
        return 0;
    }

    if (dim == 1) {
        SZ_Compress_by_config<float, 1>(argc, argv, argp, data, eb, std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        SZ_Compress_by_config<float, 2>(argc, argv, argp, data, eb, std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        SZ_Compress_by_config<float, 3>(argc, argv, argp, data, eb, std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        SZ_Compress_by_config<float, 4>(argc, argv, argp, data, eb,
                                        std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }

    return 0;
}