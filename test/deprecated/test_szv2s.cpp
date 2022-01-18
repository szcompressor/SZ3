#include "SZ3/frontend/SZMetaFrontend.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include "test_szv2.hpp"

template<typename T, uint N>
float SZ_compress_build_frontend(std::unique_ptr<T[]> const &data, const SZ::Config &conf) {
    auto quantizer = SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quant_state_num / 2);
    return SZ_compress_build_backend<T, N>(data, conf, SZ::make_sz_meta_frontend<T, N>(conf, quantizer));
}

template<class T, uint N>
float SZ_compress_parse_args(int argc, char **argv, int argp, std::unique_ptr<T[]> &data, float eb,
                             std::array<size_t, N> dims) {
    SZ::Config conf;
    conf.update_dims(dims.begin(), dims.end());
    conf.absErrorBound = eb;
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

    auto ratio = SZ_compress_build_frontend<T, N>(data, conf);
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
        SZ_compress_parse_args<float, 1>(argc, argv, argp, data, eb, std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        SZ_compress_parse_args<float, 2>(argc, argv, argp, data, eb, std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        SZ_compress_parse_args<float, 3>(argc, argv, argp, data, eb,
                                         std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        SZ_compress_parse_args<float, 4>(argc, argv, argp, data, eb,
                                         std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }


    return 0;
}
