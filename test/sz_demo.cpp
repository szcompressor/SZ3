#include "compressor/SZGeneralCompressor.hpp"
#include "compressor/PQLCompressor.hpp"
#include "frontend/Frontend.hpp"
#include "frontend/SZ3Frontend.hpp"
#include "frontend/SZMetaFrontend.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "encoder/HuffmanEncoder.hpp"
#include "encoder/ArithmeticEncoder.hpp"
#include "lossless/Lossless_zstd.hpp"
#include "lossless/Lossless_bypass.hpp"
#include "utils/FileUtil.h"
#include "utils/Config.hpp"
#include "utils/Verification.hpp"
#include "utils/Timer.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>
#include <random>

std::string src_file_name;
float relative_error_bound = 0;
float compression_time = 0;

template<typename T, class Frontend, class Lossless, uint N>
float SZ_compress_final(std::unique_ptr<T[]> const &data,
                        const SZ::Config<T, N> &conf,
                        Frontend frontend, Lossless lossless) {

    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", block_size = " << conf.block_size
              << ", stride = " << conf.stride
              << ", quan_state_num = " << conf.quant_state_num
              << std::endl
              << "lorenzo = " << conf.enable_lorenzo
              << ", 2ndlorenzo = " << conf.enable_2ndlorenzo
              << ", regression = " << conf.enable_regression
              << ", encoder = " << conf.encoder_op
              << ", lossless = " << conf.lossless_op
              << std::endl;

    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

    auto sz_huffman = make_sz_general_compressor(conf, frontend, SZ::HuffmanEncoder<int>(), lossless);
    auto sz_arithmetic = make_sz_general_compressor(conf, frontend, SZ::ArithmeticEncoder<int>(true), lossless);
//    auto sz_noencoder = make_pql_compressor(conf, predictor, SZ::LinearQuantizer<T>(conf.eb, 64),
//                                            lossless);
    SZ::Timer timer;
    timer.start();
    std::cout << "****************** Compression ******************" << std::endl;

    SZ::concepts::CompressorInterface<T> *sz;
    if (conf.encoder_op == 0) {
        std::cout << "encoder_op == 0 not support yet." << std::endl;
        return 0;
//        sz = &sz_noencoder;
    } else if (conf.encoder_op == 1) {
        sz = &sz_huffman;
    } else {
        sz = &sz_arithmetic;
    }

    size_t compressed_size = 0;
    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz->compress(data.get(), compressed_size));

    compression_time = timer.stop("Compression");

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 10000);
    std::stringstream ss;
    ss << src_file_name.substr(src_file_name.rfind('/') + 1)
       << "." << relative_error_bound << "." << dis(gen) << ".sz3";
    auto compressed_file_name = ss.str();
    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    std::cout << "Compressed file = " << compressed_file_name << std::endl;

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    T *dec_data = sz->decompress(compressed.get(), compressed_size);
    timer.stop("Decompression");

    SZ::verify<T>(data_.data(), dec_data, conf.num);

//    auto decompressed_file_name = compressed_file_name + ".out";
//    SZ::writefile(decompressed_file_name.c_str(), dec_data, conf.num);
//    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;
    delete[] dec_data;
    return ratio;
}

template<typename T, class Frontend, uint N>
float SZ_compress_step3(std::unique_ptr<T[]> const &data,
                        const SZ::Config<T, N> &conf,
                        Frontend frontend) {
    if (conf.lossless_op == 0) {
        return SZ_compress_final<T>(data, conf, frontend, SZ::Lossless_bypass());
    } else {
        return SZ_compress_final<T>(data, conf, frontend, SZ::Lossless_zstd());
    }
}

template<typename T, uint N>
float SZ_compress_step2(std::unique_ptr<T[]> const &data, const SZ::Config<T, N> &conf) {


    auto quantizer = SZ::LinearQuantizer<T>(conf.eb, conf.quant_state_num / 2);
    return SZ_compress_step3<T>(data, conf, make_sz_meta_frontend(conf, quantizer));
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
    int use_single_predictor =
            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return SZ_compress_step3<T>(data, conf,
                                        make_sz3_frontend(conf, SZ::LorenzoPredictor<T, N, 1>(conf.eb), quantizer));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return SZ_compress_step3<T>(data, conf,
                                        make_sz3_frontend(conf, SZ::LorenzoPredictor<T, N, 2>(conf.eb), quantizer));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.eb));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return SZ_compress_step3<T>(data, conf,
                                        make_sz3_frontend(conf, SZ::RegressionPredictor<T, N>(conf.block_size, conf.eb),
                                                          quantizer));
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
        }
    }

    return SZ_compress_step3<T>(data, conf,
                                make_sz3_frontend(conf, SZ::ComposedPredictor<T, N>(predictors), quantizer));
}


template<class T, uint N>
float SZ_compress_step1(int argc, char **argv, int argp, std::unique_ptr<T[]> &data, float eb,
                        std::array<size_t, N> dims) {
    SZ::Config<float, N> conf(eb, dims);
    if (argp < argc) {
        int block_size = atoi(argv[argp++]);
        conf.block_size = block_size;
        conf.stride = block_size;
    }

    int lorenzo_op = 1;
    if (argp < argc) {
        lorenzo_op = atoi(argv[argp++]);
        conf.enable_lorenzo = lorenzo_op == 1 || lorenzo_op == 3;
        conf.enable_2ndlorenzo = lorenzo_op == 2 || lorenzo_op == 3;
    }

    int regression_op = 1;
    if (argp < argc) {
        regression_op = atoi(argv[argp++]);
        conf.enable_regression = regression_op == 1;
    }

    if (argp < argc) {
        conf.encoder_op = atoi(argv[argp++]);
        if (conf.encoder_op == 2) {
            conf.quant_state_num = 1024;
        }
    }

    if (argp < argc) {
        conf.lossless_op = atoi(argv[argp++]);
    }

    if (argp < argc) {
        conf.quant_state_num = atoi(argv[argp++]);
    }

    auto ratio = SZ_compress_step2(data, conf);
    printf("%s %.0E CR=%.2f Time=%.2f L=%d R=%d E=%d Lo=%d\n",
           argv[1], relative_error_bound, ratio, compression_time,
           lorenzo_op, regression_op, conf.encoder_op, conf.lossless_op);
    return ratio;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "usage: " << argv[0] <<
                  " data_file -num_dim dim0 .. dimn relative_eb [blocksize lorenzo_op regression_op encoder_op lossless_op quant_state_num]"
                  << std::endl;
        std::cout << "example: " << argv[0] <<
                  " qmcpack.dat -3 33120 69 69 1e-3 [6 1 1 1 1 32768]" << std::endl;
        return 0;
    }

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

    if (dim == 1) {
        SZ_compress_step1<float, 1>(argc, argv, argp, data, eb, std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        SZ_compress_step1<float, 2>(argc, argv, argp, data, eb, std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        SZ_compress_step1<float, 3>(argc, argv, argp, data, eb,
                                    std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        SZ_compress_step1<float, 4>(argc, argv, argp, data, eb,
                                    std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }


    return 0;
}