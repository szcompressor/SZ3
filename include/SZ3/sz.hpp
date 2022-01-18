#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/SZBlockInterpolationCompressor.hpp"
#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZMetaFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/SimplePredictor.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Verification.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/ConfigNew.hpp"
#include <cmath>
#include <memory>


template<class T, uint N>
char *SZ_compress_interp_N(SZ::ConfigNew conf, T *data, size_t &outSize) {

//    std::cout << "****************** Interp Compression ****************" << std::endl;
//    std::cout << "Interp Op          = " << interp_op << std::endl
//              << "Direction          = " << direction_op << std::endl
//              << "SZ block size      = " << block_size << std::endl
//              << "Interp block size  = " << interp_block_size << std::endl;

    assert(N == conf.N);
    if (conf.errorBoundMode != ABS) {
        if (conf.errorBoundMode == REL) {
            conf.absErrorBound = conf.relErrorBound * SZ::data_range(data, conf.num);
        } else {
            printf("Error, error bound mode not supported\n");
            exit(0);
        }
    }
    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(conf.absErrorBound),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    return (char *) sz.compress(conf, data, outSize);
}


template<class T>
char *SZ_compress_interp(SZ::ConfigNew conf, T *data, size_t &outSize) {
    if (conf.N == 1) {
        return SZ_compress_interp_N<T, 1>(conf, data, outSize);
    } else if (conf.N == 2) {
        return SZ_compress_interp_N<T, 2>(conf, data, outSize);
    } else if (conf.N == 3) {
        return SZ_compress_interp_N<T, 3>(conf, data, outSize);
    } else if (conf.N == 4) {
        return SZ_compress_interp_N<T, 4>(conf, data, outSize);
    } else {
        for (int i = 4; i < conf.N; i++) {
            conf.dims[3] *= conf.dims[i];
        }
        conf.dims.resize(4);
        conf.N = 4;
        return SZ_compress_interp_N<T, 4>(conf, data, outSize);
    }
}


template<class T, uint N>
T *SZ_decompress_interp_N(const SZ::ConfigNew &conf, char *cmpData, size_t cmpSize) {
    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    return sz.decompress(conf, (SZ::uchar *) cmpData, cmpSize);
}

template<class T>
T *SZ_decompress_interp(const SZ::ConfigNew &conf, char *cmpData, size_t cmpSize) {
    if (conf.N == 1) {
        return SZ_decompress_interp_N<T, 1>(conf, cmpData, cmpSize);
    } else if (conf.N == 2) {
        return SZ_decompress_interp_N<T, 2>(conf, cmpData, cmpSize);
    } else if (conf.N == 3) {
        return SZ_decompress_interp_N<T, 3>(conf, cmpData, cmpSize);
    } else if (conf.N == 4) {
        return SZ_decompress_interp_N<T, 4>(conf, cmpData, cmpSize);
    } else {
        return SZ_decompress_interp_N<T, 4>(conf, cmpData, cmpSize);
    }
}