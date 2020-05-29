//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_SZ_HPP
#define SZ_SZ_HPP


#include "SZ_Impl.hpp"
//#include "def.hpp"
#include "utils/Config.hpp"
#include <cstdio>
#include <iostream>
#include <memory>

template<typename T, uint N>
float SZ_Compress(std::unique_ptr<T[]> const &data, const SZ::Config<T, N> &conf) {
    std::cout << "Options: "
              << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", block_size = " << conf.block_size
              << ", stride = " << conf.stride
              << "  use_lorenzo = " << conf.enable_lorenzo
              << ", use_2ndlorenzo = " << conf.enable_2ndlorenzo
              << ", use_regression = " << conf.enable_regression
              << ", use_2ndregression = " << conf.enable_2ndregression
              << std::endl
              << "  use_lossless= " << conf.enable_lossless
              << ", quan_bin= " << conf.quant_bin
              << ", zone = " << conf.zone
              << ", decompress_zone_idx = " << conf.decompress_zone_idx
              << std::endl;

    std::vector<std::shared_ptr<SZ::concepts::VirtualPredictor<T, N>>> predictors;
    int use_single_predictor =
            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression + conf.enable_2ndregression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return SZ::SZ_Compress_Impl<T>(data, conf, SZ::LorenzoPredictor<T, N, 1>(conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return SZ::SZ_Compress_Impl<T>(data, conf, SZ::LorenzoPredictor<T, N, 2>(conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.eb));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return SZ::SZ_Compress_Impl<T>(data, conf, SZ::RegressionPredictor<T, N>(conf.block_size, conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
        }
    }
    if (conf.enable_2ndregression) {
        if (use_single_predictor) {
            return SZ::SZ_Compress_Impl<T>(data, conf, SZ::PolyRegressionPredictor<T, N>(conf.block_size, conf.eb));
        } else {
            predictors.push_back(std::make_shared<SZ::PolyRegressionPredictor<T, N>>(conf.block_size, conf.eb));
        }
    }
    return SZ::SZ_Compress_Impl<T>(data, conf, SZ::ComposedPredictor<T, N>(predictors));
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

#endif //SZ_SZ_HPP
