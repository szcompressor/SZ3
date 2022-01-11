//
// Created by Kai Zhao on 3/2/21.
//

#ifndef SZ3_SZ_V2_HPP
#define SZ3_SZ_V2_HPP

#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/ArithmeticEncoder.hpp"
#include "SZ3/encoder/BypassEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/utils/Verification.hpp"
#include "SZ3/utils/Timer.hpp"
#include <sstream>
#include <random>

std::string src_file_name;
float relative_error_bound = 0;

template<typename T, class Frontend, class Encoder, class Lossless, uint N>
float SZ_compress(std::unique_ptr<T[]> const &data,
                  const SZ::Config<T, N> &conf,
                  Frontend frontend, Encoder encoder, Lossless lossless) {

    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", block_size = " << conf.block_size
              << ", stride = " << conf.stride
              << ", quan_state_num = " << conf.quant_state_num
              << std::endl
              << "lorenzo = " << conf.enable_lorenzo
              << ", 2ndlorenzo = " << conf.enable_2ndlorenzo
              << ", regression = " << conf.enable_regression
              << ", 2ndregression = " << conf.enable_2ndregression
              << ", encoder = " << conf.encoder_op
              << ", lossless = " << conf.lossless_op
              << std::endl;

    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

    auto sz = make_sz_general_compressor(conf, frontend, encoder, lossless);

    SZ::Timer timer(true);
    std::cout << "****************** Compression ******************" << std::endl;

    size_t compressed_size = 0;
    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz.compress(data.get(), compressed_size));

    timer.stop("Compression");

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    auto compressed_file_name = SZ::compressed_path(src_file_name, false);
    std::cout << "Compression ratio = " << ratio << std::endl;
    std::cout << "Compression size = " << compressed_size << std::endl;
    std::cout << "Compression file = " << compressed_file_name << std::endl;
    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    T *dec_data = sz.decompress(compressed.get(), compressed_size);
    timer.stop("Decompression");

//    remove(compressed_file_name.c_str());
    auto decompressed_file_name = SZ::decompressed_path(src_file_name, false);
    SZ::writefile(decompressed_file_name.c_str(), dec_data, conf.num);
    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

    SZ::verify<T>(data_.data(), dec_data, conf.num);

    delete[] dec_data;
    return ratio;
}

template<typename T, class Frontend, uint N>
float SZ_compress_build_backend(std::unique_ptr<T[]> const &data,
                                const SZ::Config<T, N> &conf,
                                Frontend frontend) {
    if (conf.lossless_op == 1) {
        if (conf.encoder_op == 1) {
            return SZ_compress<T>(data, conf, frontend, SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());
        } else if (conf.encoder_op == 2) {
            return SZ_compress<T>(data, conf, frontend, SZ::ArithmeticEncoder<int>(true), SZ::Lossless_zstd());
        } else {
            return SZ_compress<T>(data, conf, frontend, SZ::BypassEncoder<int>(), SZ::Lossless_zstd());
        }
    } else {
        if (conf.encoder_op == 1) {
            return SZ_compress<T>(data, conf, frontend, SZ::HuffmanEncoder<int>(), SZ::Lossless_bypass());
        } else if (conf.encoder_op == 2) {
            return SZ_compress<T>(data, conf, frontend, SZ::ArithmeticEncoder<int>(true), SZ::Lossless_bypass());
        } else {
            return SZ_compress<T>(data, conf, frontend, SZ::BypassEncoder<int>(), SZ::Lossless_bypass());
        }
    }
}

#endif //SZ3_SZ_V2_HPP

