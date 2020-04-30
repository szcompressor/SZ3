//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_SZ_IMPL_HPP
#define SZ_SZ_IMPL_HPP


#include "quantizer/IntegerQuantizer.hpp"
#include "compressor/SZ_General_Compressor.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/PolyRegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "utils/fileUtil.h"
#include "utils/Config.hpp"
#include "utils/Verification.hpp"
#include <def.hpp>
#include <cstdio>
#include <iostream>
#include <memory>

namespace SZ {

    template<typename T, class Predictor, uint N>
    float SZ_Compress_step2(std::unique_ptr<T[]> const &data,
                            const Config<T, N> &conf,
                            Predictor predictor) {
        std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);
        std::cout << "****************** Compression ******************" << std::endl;

        auto sz = SZ::make_sz_general_compressor(conf, predictor, SZ::LinearQuantizer<T>(conf.eb), SZ::HuffmanEncoder<int>());

        size_t compressed_size = 0;
        std::unique_ptr<SZ::uchar[]> compressed;
        compressed.reset(sz.compress(data.get(), compressed_size));

        clock_gettime(CLOCK_REALTIME, &end);
        std::cout << "Compression time = "
                  << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                  << "s"
                  << std::endl;

        auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
        std::cout << "Compression Ratio = " << ratio << std::endl;
        SZ::writefile("compressed.out", compressed.get(), compressed_size);

        std::cout << "****************** Decompression ******************" << std::endl;
        compressed = SZ::readfile<SZ::uchar>("compressed.out", compressed_size);

        clock_gettime(CLOCK_REALTIME, &start);
        std::unique_ptr<T[]> dec_data;
        dec_data.reset(sz.decompress(compressed.get(), compressed_size));
        clock_gettime(CLOCK_REALTIME, &end);
        std::cout << "Decompression time: "
                  << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                  << "s"
                  << std::endl;

        SZ::verify<T>(data_.data(), data.get(), conf.num);
        return ratio;
    }

    template<typename T, uint N>
    float SZ_Compress_step1(std::unique_ptr<T[]> const &data, const Config<T, N> &conf) {
        std::vector<std::shared_ptr<SZ::concepts::VirtualPredictor<T, N>>> predictors;
        int use_single_predictor =
                (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression + conf.enable_2ndregression) == 1;
        if (conf.enable_lorenzo) {
            if (use_single_predictor) {
                return SZ_Compress_step2<T>(data, conf, SZ::LorenzoPredictor<T, N, 1>(conf.eb));
            } else {
                predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
            }
        }
        if (conf.enable_2ndlorenzo) {
            if (use_single_predictor) {
                return SZ_Compress_step2<T>(data, conf, SZ::LorenzoPredictor<T, N, 2>(conf.eb));
            } else {
                predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.eb));
            }
        }
        if (conf.enable_regression) {
            if (use_single_predictor) {
                return SZ_Compress_step2<T>(data, conf, SZ::RegressionPredictor<T, N>(conf.block_size, conf.eb));
            } else {
                predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
            }
        }
        if (conf.enable_2ndregression) {
            if (use_single_predictor) {
                return SZ_Compress_step2<T>(data, conf, SZ::PolyRegressionPredictor<T, N>(conf.block_size, conf.eb));
            } else {
                predictors.push_back(std::make_shared<SZ::PolyRegressionPredictor<T, N>>(conf.block_size, conf.eb));
            }
        }
        return SZ_Compress_step2<T>(data, conf, SZ::ComposedPredictor<T, N>(predictors));
    }
}

#endif //SZ_SZ_IMPL_HPP
