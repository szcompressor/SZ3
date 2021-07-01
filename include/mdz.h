//
// Created by Kai Zhao on 7/1/21.
//

#ifndef SZ3_MDZ_H
#define SZ3_MDZ_H

#include <vector>
//#include <sstream>
//#include "utils/FileUtil.h"
#include "utils/Timer.hpp"
#include "def.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "lossless/Lossless_zstd.hpp"
//#include "lossless/Lossless_bypass.hpp"
#include "encoder/HuffmanEncoder.hpp"
//#include "encoder/BypassEncoder.hpp"
#include "compressor/SZExaaltCompressor.hpp"
#include "compressor/SZGeneralCompressor.hpp"
#include "frontend/SZ3TimeBasedFrontend.hpp"
#include "frontend/SZ3Frontend.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "utils/KmeansUtil.h"
//#include "utils/QuantOptimizatioin.hpp"
//#include <cstdio>



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

template<typename T, uint N, class Predictor>
SZ::concepts::CompressorInterface<T> *
make_sz2(const SZ::Config<T, N> &conf, Predictor predictor) {

    return new SZ::SZGeneralCompressor<T, N, SZ::SZ3Frontend<T, N, Predictor, SZ::LinearQuantizer<T>>,
            SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            conf, make_sz3_frontend(conf, predictor, SZ::LinearQuantizer<T>(conf.eb, conf.quant_state_num / 2)),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
}

template<typename T, uint N>
SZ::concepts::CompressorInterface<T> *
make_sz(const SZ::Config<T, N> &conf) {
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;

    int use_single_predictor =
            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return make_sz2(conf, SZ::LorenzoPredictor<T, N, 1>(conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return make_sz2(conf, SZ::LorenzoPredictor<T, N, 2>(conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.eb));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return make_sz2(conf, SZ::RegressionPredictor<T, N>(conf.block_size, conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N >>(conf.block_size, conf.eb));
        }
    }

    return make_sz2(conf, SZ::ComposedPredictor<T, N>(predictors));
}

template<typename T, uint N>
SZ::uchar *compress(SZ::Config<T, N> conf, T *data, int method, size_t &compressed_size,
                    float level_start, float level_offset, int level_num, T *ts0) {
    if ((method == 0 || method == 1) && level_num == 0) {
        printf("VQ/VQT not availble on current dataset, please use ADP or MT\n");
        exit(0);
    }
    SZ::uchar *compressed_data;
    if (method == 0 || method == 1) {
        auto sz = SZ::SZ_Exaalt_Compressor<T, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                conf, SZ::LinearQuantizer<float>(conf.eb, conf.quant_state_num / 2),
                SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd(), method);
        sz.set_level(level_start, level_offset, level_num);
        compressed_data = sz.compress(data, compressed_size);
    } else if (method == 2 || method == 4) {
        auto sz = make_sz_timebased(conf, ts0);
        compressed_data = sz->compress(data, compressed_size);
    } else {
        auto sz = make_sz(conf);
        compressed_data = sz->compress(data, compressed_size);
    }
//    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
//    std::cout << "Compression Ratio = " << ratio << std::endl;
//    std::cout << "Compressed size = " << compressed_size << std::endl;
    return compressed_data;
}

template<typename T, uint N>
int select_compressor(SZ::Config<T, N> conf, T *data, bool firsttime,
                      float level_start, float level_offset, int level_num, T *data_ts0) {
    std::vector<size_t> compressed_size(10, std::numeric_limits<size_t>::max());
    if (conf.dims[0] > 10) {
        conf.dims[0] = 10;
        conf.num = conf.dims[0] * conf.dims[1];
    }
    if (firsttime) {
        conf.dims[0] /= 2;
        conf.num = conf.dims[0] * conf.dims[1];
        data += conf.num;
    }
    std::vector<T> data1;
    SZ::uchar *cmpr;
    if (level_num > 0) {

        data1 = std::vector<T>(data, data + conf.num);
        cmpr = compress(conf, data1.data(), 0, compressed_size[0], level_start, level_offset, level_num, data_ts0);
        delete[] cmpr;

        data1 = std::vector<T>(data, data + conf.num);
        cmpr = compress(conf, data1.data(), 1, compressed_size[1], level_start, level_offset, level_num, data_ts0);
        delete[] cmpr;
    } else {
        data1 = std::vector<T>(data, data + conf.num);
        cmpr = compress(conf, data1.data(), 3, compressed_size[3], level_start, level_offset, level_num, data_ts0);
        delete[] cmpr;
    }

    data1 = std::vector<T>(data, data + conf.num);
    cmpr = compress(conf, data1.data(), 2, compressed_size[2], level_start, level_offset, level_num, data_ts0);
    delete[] cmpr;

//    data1 = std::vector(&data[t * conf.dims[1]], &data[t * conf.dims[1]] + conf.num);
//    MT(conf, t, data1.data(), compressed_size[4], false, (T *) nullptr);

    int method = std::distance(compressed_size.begin(),
                               std::min_element(compressed_size.begin(), compressed_size.end()));
    const char *compressor_names[] = {"VQ", "VQT", "MT", "LR", "TS"};
    printf("Select %s as Compressor, method=%d, %lu %lu %lu %lu %lu\n",
           compressor_names[method], method,
           compressed_size[0], compressed_size[1], compressed_size[2], compressed_size[3], compressed_size[4]);
    std::cout << "****************** END Selection ****************" << std::endl;
    return method;
}

#endif //SZ3_MDZ_H
