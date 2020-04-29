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
              << "eb = " << conf.eb
              << ", block_size = " << conf.block_size
              << ", stride = " << conf.stride
              //                  << ", dim = " << conf.pred_dim
              << ", lorenzo_1 = " << conf.enable_lorenzo
              << ", lorenzo_2 = " << conf.enable_2ndlorenzo
              << ", regression_1 = " << conf.enable_regression
              << ", regression_2 = " << conf.enable_2ndregression
              << ", lossless= " << conf.enable_lossless
              << std::endl;

    return SZ::SZ_Compress_step1<T, N>(data, conf);
}

template<typename T>
float SZ_Compress(std::unique_ptr<T[]> const &data, T eb, uint r1, uint r2, uint r3, uint r4) {
    return SZ_Compress(data, SZ::Config<T, 3>(eb, std::array<size_t, 3>{r1 * r2, r3, r4}));
}

template<typename T>
float SZ_Compress(std::unique_ptr<T[]> const &data, T eb, uint r1, uint r2, uint r3) {
    return SZ_Compress(data, SZ::Config<T, 3>(eb, std::array<size_t, 3>{r1, r2, r3}));
}

template<typename T>
float SZ_Compress(std::unique_ptr<T[]> const &data, T eb, uint r1, uint r2) {
    return SZ_Compress(data, SZ::Config<T, 2>(eb, std::array<size_t, 2>{r1, r2}));
}

template<typename T>
float SZ_Compress(std::unique_ptr<T[]> const &data, T eb, uint r1) {
    return SZ_Compress(data, SZ::Config<T, 1>(eb, std::array<size_t, 1>{r1}));
}

#endif //SZ_SZ_HPP
