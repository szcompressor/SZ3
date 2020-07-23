//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_SZ_HPP
#define SZ_SZ_HPP

#include "compressor/SZ_General_Compressor.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/PolyRegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "encoder/HuffmanEncoder.hpp"
#include "lossless/Lossless_zstd.hpp"
#include "lossless/Lossless_bypass.hpp"
#include "lossless/Lossless_lzfse.hpp"
#include "lossless/Lossless_lizard.hpp"
#include "lossless/Lossless_tornado.hpp"
#include "utils/fileUtil.h"
#include "utils/Config.hpp"
#include "utils/Verification.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>

namespace SZ {

    template<typename T, class Predictor, uint N>
    float SZ_Compress(std::unique_ptr<T[]> const &data,
                      const Config<T, N> &conf,
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
                  << ", 2ndregression = " << conf.enable_2ndregression
                  << ", lossless = " << conf.enable_lossless
                  << std::endl;

        std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);
        std::cout << "****************** Compression ******************" << std::endl;

        auto sz = SZ::make_sz_general_compressor(conf, predictor, SZ::LinearQuantizer<T>(conf.eb, conf.quant_bin),
                                                 SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());

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
        SZ::writefile("compressed.dat", compressed.get(), compressed_size);

        std::cout << "****************** Decompression ****************" << std::endl;
        compressed = SZ::readfile<SZ::uchar>("compressed.dat", compressed_size);

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
    float SZ_Compress(std::unique_ptr<T[]> const &data, const SZ::Config<T, N> &conf) {


        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
        int use_single_predictor =
                (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression + conf.enable_2ndregression) == 1;
        if (conf.enable_lorenzo) {
            if (use_single_predictor) {
                return SZ::SZ_Compress<T>(data, conf, SZ::LorenzoPredictor<T, N, 1>(conf.eb));
            } else {
                predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
            }
        }
        if (conf.enable_2ndlorenzo) {
            if (use_single_predictor) {
                return SZ::SZ_Compress<T>(data, conf, SZ::LorenzoPredictor<T, N, 2>(conf.eb));
            } else {
                predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.eb));
            }
        }
        if (conf.enable_regression) {
            if (use_single_predictor) {
                return SZ::SZ_Compress<T>(data, conf, SZ::RegressionPredictor<T, N>(conf.block_size, conf.eb));
            } else {
                predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
            }
        }
        if (conf.enable_2ndregression) {
            if (use_single_predictor) {
                return SZ::SZ_Compress<T>(data, conf, SZ::PolyRegressionPredictor<T, N>(conf.block_size, conf.eb));
            } else {
                predictors.push_back(std::make_shared<SZ::PolyRegressionPredictor<T, N>>(conf.block_size, conf.eb));
            }
        }
        return SZ::SZ_Compress<T>(data, conf, SZ::ComposedPredictor<T, N>(predictors));
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


}
#endif //SZ_SZ_HPP
