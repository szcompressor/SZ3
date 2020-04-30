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
    float SZ_Compress_Impl(std::unique_ptr<T[]> const &data,
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
}

#endif //SZ_SZ_IMPL_HPP
