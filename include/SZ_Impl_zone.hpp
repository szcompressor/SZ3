//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_SZ_IMPL_ZONE_HPP
#define SZ_SZ_IMPL_ZONE_HPP


#include "quantizer/IntegerQuantizer.hpp"
#include "compressor/SZ_Zone_Compressor.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/PolyRegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "utils/fileUtil.h"
#include "utils/Config.hpp"
#include "utils/Verification.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>

namespace SZ {

    template<typename T, class Predictor, uint N>
    float SZ_Compress_Impl_zone(std::unique_ptr<T[]> const &data,
                                const Config<T, N> &conf,
                                Predictor predictor) {
        std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);
        std::cout << "****************** Compression ******************" << std::endl;

        auto sz = SZ::make_sz_zone_compressor(conf, predictor, SZ::LinearQuantizer<T>(conf.eb), SZ::HuffmanEncoder<int>());

        std::vector<std::unique_ptr<SZ::uchar[]>> compressed(conf.zone + 1);
        std::vector<size_t> compressed_size(conf.zone + 1, 0);
        auto szresult = sz.compress(data.get(), compressed_size, conf.zone);
        for (int i = 0; i < conf.zone + 1; i++) {
            compressed[i].reset(szresult[i]);
        }

        clock_gettime(CLOCK_REALTIME, &end);
        std::cout << "Compression time = "
                  << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                  << "s"
                  << std::endl;

        auto compressed_size_total = std::accumulate(compressed_size.begin(), compressed_size.end(), 0);
        auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size_total;
        std::cout << "Compression Ratio = " << ratio << std::endl;
        for (int i = 0; i < conf.zone; i++) {
            std::string filename = "compressed.data." + std::to_string(i) + ".out";
            SZ::writefile(filename.data(), compressed[i].get(), compressed_size[i]);
        }
        SZ::writefile("compressed.encoder.out", compressed[conf.zone].get(), compressed_size[conf.zone]);

        std::cout << "****************** Decompression ******************" << std::endl;
        size_t compressed_encoder_size, compressed_data_size;
        auto compressed_encoder = SZ::readfile<SZ::uchar>("compressed.encoder.out", compressed_encoder_size);
        std::string filename = "compressed.data." + std::to_string(conf.decompress_zone_idx) + ".out";
        auto compressed_data = SZ::readfile<SZ::uchar>(filename.data(), compressed_data_size);

        clock_gettime(CLOCK_REALTIME, &start);
        sz.decompress_encoder(compressed_encoder.get(), compressed_encoder_size);
        clock_gettime(CLOCK_REALTIME, &end);
        std::cout << "Decompression encoder time = "
                  << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                  << "s"
                  << std::endl;

        clock_gettime(CLOCK_REALTIME, &start);
//        std::unique_ptr<T[]> dec_data;
        std::vector<T> dec_data;
//        dec_data.reset(sz.decompress_zone(compressed_data.get(), compressed_data_size));
        sz.decompress_zone(compressed_data.get(), compressed_data_size, dec_data);
        clock_gettime(CLOCK_REALTIME, &end);
        std::cout << "Decompression time = "
                  << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                  << "s"
                  << std::endl;

//        SZ::verify<T>(data_.data(), data.get(), conf.num);
        return ratio;
    }
}

#endif //SZ_SZ_IMPL_ZONE_HPP
