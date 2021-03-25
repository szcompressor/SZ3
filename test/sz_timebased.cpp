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


template<typename T, class Predictor, uint N>
float SZ_Compress(std::unique_ptr<T[]> const &data, SZ::Config<T, N> conf, Predictor predictor) {
    assert(N == 2);
    if (conf.timestep_batch == 0) {
        conf.timestep_batch = conf.dims[0];
    }

    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", timestep_batch = " << conf.timestep_batch
              << ", quan_state_num = " << conf.quant_state_num
              << ", encoder = " << conf.encoder_op
              << ", lossless = " << conf.lossless_op
              << std::endl;

    std::vector<T> data_(data.get(), data.get() + conf.num);
    std::vector<T> dec_data(conf.num);

    double total_compressed_size = 0;
    auto dims = conf.dims;
    auto num = conf.num;

    for (size_t ts = 0; ts < dims[0]; ts += conf.timestep_batch) {
        conf.dims[0] = (ts + conf.timestep_batch - 1 > dims[0] ? dims[0] - ts + 1 : conf.timestep_batch);
        conf.num = conf.dims[0] * conf.dims[1];

        std::cout << "****************** Compression From " << ts << " to " << ts + conf.dims[0] - 1
                  << " ******************" << std::endl;
        SZ::Timer timer(true);
        auto sz = SZ::make_sz_general_compressor(
                conf,
                make_sz3_timebased_frontend(conf, predictor, SZ::LinearQuantizer<T>(conf.eb, conf.quant_state_num / 2)),
                SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());

        size_t compressed_size = 0;
        std::unique_ptr<SZ::uchar[]> compressed;
        auto ts_data = data.get() + ts * dims[1];
        compressed.reset(sz.compress(ts_data, compressed_size));
        timer.stop("Compression");

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
        auto ts_dec_data = sz.decompress(compressed.get(), compressed_size);
        timer.stop("Decompression");
        std::copy(ts_dec_data, ts_dec_data + conf.num, dec_data.begin() + ts * dims[1]);

        auto decompressed_file_name = compressed_file_name + ".out";
//    SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;
        remove(compressed_file_name.c_str());
        total_compressed_size += compressed_size;
    }
    SZ::verify<T>(data_.data(), dec_data.data(), num);

    std::cout << "****************** Final ****************" << std::endl;
    float ratio = num * sizeof(T) / total_compressed_size;
    std::cout << "Total Compression Size = " << total_compressed_size << std::endl;
    std::cout << "Total Compression Ratio = " << ratio << std::endl;
    return ratio;
}

template<typename T, uint N>
float SZ_compress_build_frontend(std::unique_ptr<T[]> const &data, const SZ::Config<T, N> &conf) {
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N - 1>>> predictors;

    int use_single_predictor =
            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return SZ_Compress<T>(data, conf, SZ::LorenzoPredictor<T, N - 1, 1>(conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N - 1, 1>>(conf.eb));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return SZ_Compress<T>(data, conf, SZ::LorenzoPredictor<T, N - 1, 2>(conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N - 1, 2>>(conf.eb));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return SZ_Compress<T>(data, conf, SZ::RegressionPredictor<T, N - 1>(conf.block_size, conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N - 1>>(conf.block_size, conf.eb));
        }
    }

    return SZ_Compress<T>(data, conf, SZ::ComposedPredictor<T, N - 1>(predictors));
}


template<class T, uint N>
float SZ_compress_parse_args(int argc, char **argv, int argp, std::array<size_t, N> dims) {
    SZ::Config<float, N> conf(dims);

    conf.src_file_name = argv[1];
    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
    std::cout << "Read " << num << " elements\n";
    assert(conf.num == num);

    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    char *eb_op = argv[argp++] + 1;
    if (*eb_op == 'a') {
        conf.eb = atof(argv[argp++]);
        conf.relative_eb = conf.eb / (max - min);
    } else {
        conf.relative_eb = atof(argv[argp++]);
        conf.eb = conf.relative_eb * (max - min);
    }

    if (argp < argc) {
        conf.timestep_batch = atoi(argv[argp++]);
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

    auto ratio = SZ_compress_build_frontend(data, conf);
    return ratio;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "SZ v" << SZ_versionString() << std::endl;
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