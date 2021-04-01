#include "compressor/SZGeneralCompressor.hpp"
#include "frontend/SZ3TimeBasedFrontend.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "lossless/Lossless_zstd.hpp"
#include "utils/FileUtil.h"
#include "utils/Config.hpp"
#include "utils/Verification.hpp"
#include "utils/Timer.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>

template<typename T, uint N, class Predictor>
SZ::concepts::CompressorInterface<T> *
make_sz_timebased2(const SZ::Config<T, N> &conf, Predictor predictor, T *data_ts0) {

    return new SZ::SZGeneralCompressor<T, N, SZ::SZ3TimeBasedFrontend<T, N, Predictor, SZ::LinearQuantizer<T>>,
            SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            conf,
            make_sz3_timebased_frontend(
                    conf, predictor, SZ::LinearQuantizer<T>(conf.eb, conf.quant_state_num / 2), data_ts0),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
}

template<typename T, uint N>
SZ::concepts::CompressorInterface<T> *
make_sz_timebased(const SZ::Config<T, N> &conf, T *data_ts0) {
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N - 1>>> predictors;

    int use_single_predictor =
            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return make_sz_timebased2(conf, SZ::LorenzoPredictor<T, N - 1, 1>(conf.eb), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N - 1, 1>>(conf.eb));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return make_sz_timebased2(conf, SZ::LorenzoPredictor<T, N - 1, 2>(conf.eb), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N - 1, 2>>(conf.eb));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return make_sz_timebased2(conf, SZ::RegressionPredictor<T, N - 1>(conf.block_size, conf.eb), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N - 1>>(conf.block_size, conf.eb));
        }
    }

    return make_sz_timebased2(conf, SZ::ComposedPredictor<T, N - 1>(predictors), data_ts0);
}

template<typename T, uint N>
float Compress(SZ::Config<T, N> conf) {
    assert(N == 2);
    if (conf.timestep_batch == 0) {
        conf.timestep_batch = conf.dims[0];
    }
    std::unique_ptr<T[]> data_ts0;
    if (conf.timestep_op == 1) {
        auto ts0_name_str = conf.src_file_name + ".ts0";
        auto ts0_name = ts0_name_str.data();
        if (SZ::file_exist(ts0_name)) {
            size_t num_ts0;
            data_ts0 = SZ::readfile<T>(ts0_name, num_ts0);
        } else {
            conf.timestep_op = 0;
        }
    }
    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << (conf.eb_mode == 0 ? conf.eb : conf.relative_eb)
              << ", timestep_batch = " << conf.timestep_batch
              << ", timestep_op = " << conf.timestep_op
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

    double total_compressed_size = 0;
    double total_compress_time = 0;
    double total_decompress_time = 0;
    auto dims = conf.dims;
    auto total_num = conf.num;
    std::vector<T> dec_data(total_num);

    for (size_t ts = 0; ts < dims[0]; ts += conf.timestep_batch) {
        conf.dims[0] = (ts + conf.timestep_batch > dims[0] ? dims[0] - ts : conf.timestep_batch);
        conf.num = conf.dims[0] * conf.dims[1];

        auto data = SZ::readfile<T>(conf.src_file_name.data(), ts * conf.dims[1], conf.num);
        T max = *std::max_element(data.get(), data.get() + conf.num);
        T min = *std::min_element(data.get(), data.get() + conf.num);
        if (conf.eb_mode == 0) {
            conf.relative_eb = conf.eb / (max - min);
        } else if (conf.eb_mode == 1) {
            conf.eb = conf.relative_eb * (max - min);
        }

        std::cout << "****************** Compression From " << ts << " to " << ts + conf.dims[0] - 1
                  << " ******************" << std::endl;
        SZ::Timer timer(true);
        auto sz = make_sz_timebased(conf, data_ts0.get());

        size_t compressed_size = 0;
        std::unique_ptr<SZ::uchar[]> compressed;
        compressed.reset(sz->compress(data.get(), compressed_size));
        total_compress_time += timer.stop("Compression");

        auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
        std::cout << "Compression Ratio = " << ratio << std::endl;
        std::cout << "Compressed size = " << compressed_size << std::endl;

        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> dis(0, 10000);
        std::stringstream ss;
        ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
           << "." << conf.relative_eb << "." << dis(gen) << ".sz3";
        auto compressed_file_name = ss.str();
        SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
        std::cout << "Compressed file = " << compressed_file_name << std::endl;

        std::cout << "****************** Decompression ****************" << std::endl;
        compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

        timer.start();
        auto ts_dec_data = sz->decompress(compressed.get(), compressed_size);
        total_decompress_time += timer.stop("Decompression");

        auto decompressed_file_name = compressed_file_name + ".out";
//    SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;
        remove(compressed_file_name.c_str());
        memcpy(&dec_data[ts * conf.dims[1]], ts_dec_data, conf.num * sizeof(T));
        total_compressed_size += compressed_size;
        delete sz;

    }

    std::cout << "****************** Final ****************" << std::endl;
    float ratio = total_num * sizeof(T) / total_compressed_size;
    std::cout << "Total Compression Ratio = " << ratio << std::endl;
    std::cout << "Total Compression Time = " << total_compress_time << std::endl;
    std::cout << "Total Decompression Time = " << total_decompress_time << std::endl;
    std::cout << "Total Compression Size = " << total_compressed_size << std::endl;
    auto data = SZ::readfile<T>(conf.src_file_name.data(), 0, total_num);
    SZ::verify<T>(data.get(), dec_data.data(), conf.num);
    return ratio;
}


template<class T, uint N>
float SZ_compress_parse_args(int argc, char **argv, int argp, std::array<size_t, N> dims) {
    SZ::Config<float, N> conf(dims);

    conf.src_file_name = argv[1];

    char *eb_op = argv[argp++] + 1;
    if (*eb_op == 'a') {
        conf.eb_mode = 0;
        conf.eb = atof(argv[argp++]);
    } else {
        conf.eb_mode = 1;
        conf.relative_eb = atof(argv[argp++]);
    }

    if (argp < argc) {
        conf.timestep_batch = atoi(argv[argp++]);
    }
    if (argp < argc) {
        conf.timestep_op = atoi(argv[argp++]);
    }
    conf.block_size = 128;
    conf.stride = 128;
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
        if (conf.encoder_op == 0) {
            conf.quant_state_num = 250;
        } else if (conf.encoder_op == 2) {
            conf.quant_state_num = 1024;
        }
    }

    if (argp < argc) {
        conf.lossless_op = atoi(argv[argp++]);
    }

    if (argp < argc) {
        conf.quant_state_num = atoi(argv[argp++]);
    }

    auto ratio = Compress(conf);
    return ratio;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "SZ" << SZ_versionString() << std::endl;
        std::cout << "usage: " << argv[0] <<
                  " data_file -num_dim dim0 .. dimn relative_eb [blocksize lorenzo_op regression_op encoder_op lossless_op quant_state_num]"
                  << std::endl;
        std::cout << "example: " << argv[0] <<
                  " qmcpack.dat -3 33120 69 69 1e-3 [6 1 1 1 1 32768]" << std::endl;
        return 0;
    }


    int dim = atoi(argv[2] + 1);
    assert(2 <= dim && dim <= 2);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

//    if (dim == 1) {
//        SZ_compress_parse_args<float, 1>(argc, argv, argp, data, eb, std::array<size_t, 1>{dims[0]});
//    } else if (dim == 2) {
    SZ_compress_parse_args<float, 2>(argc, argv, argp, std::array<size_t, 2>{dims[0], dims[1]});
//    } else if (dim == 3) {
//        SZ_compress_parse_args<float, 3>(argc, argv, argp, data, eb,
//                                         std::array<size_t, 3>{dims[0], dims[1], dims[2]});
//    } else if (dim == 4) {
//        SZ_compress_parse_args<float, 4>(argc, argv, argp, data, eb,
//                                         std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
//    }


    return 0;
}