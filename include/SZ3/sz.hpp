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
#include "SZ3/utils/Config.hpp"
#include <cmath>
#include <memory>


template<class T, uint N>
char *SZ_compress_N(SZ::Config conf, T *data, size_t &outSize) {

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

    char *cmpData;
    if (conf.cmprMethod == METHOD_LORENZO) {
        auto quantizer = SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quant_state_num / 2);
        auto sz = SZ::make_sz_general_compressor<T, N>(conf, SZ::make_sz_meta_frontend<T, N>(conf, quantizer), SZ::HuffmanEncoder<int>(),
                                                       SZ::Lossless_zstd());
        cmpData = (char *) sz.compress(conf, data, outSize);

    } else if (conf.cmprMethod == METHOD_INTERP) {
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(conf.absErrorBound),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd());
        cmpData = (char *) sz.compress(conf, data, outSize);
    }
    cmpData[outSize] = conf.cmprMethod;
    outSize += 1;
    return cmpData;
}


template<class T, uint N>
T *SZ_decompress_N(const SZ::Config &conf, char *cmpData, size_t cmpSize) {
    char cmprMethod = cmpData[cmpSize - 1];
    SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData;

    if (cmprMethod == METHOD_LORENZO) {
        SZ::LinearQuantizer<T> quantizer;
        auto sz = SZ::make_sz_general_compressor<T, N>(conf, SZ::make_sz_meta_frontend<T, N>(conf, quantizer), SZ::HuffmanEncoder<int>(),
                                                       SZ::Lossless_zstd());
        return sz.decompress(cmpDataPos, cmpSize, conf.num);

    } else {
//        if (cmprMethod == METHOD_INTERP) {
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd());
        return sz.decompress(cmpDataPos, cmpSize, conf.num);

    }

}


template<class T>
char *SZ_compress(SZ::Config conf, T *data, size_t &outSize) {
    if (conf.N == 1) {
        return SZ_compress_N<T, 1>(conf, data, outSize);
    } else if (conf.N == 2) {
        return SZ_compress_N<T, 2>(conf, data, outSize);
    } else if (conf.N == 3) {
        return SZ_compress_N<T, 3>(conf, data, outSize);
    } else if (conf.N == 4) {
        return SZ_compress_N<T, 4>(conf, data, outSize);
    } else {
        for (int i = 4; i < conf.N; i++) {
            conf.dims[3] *= conf.dims[i];
        }
        conf.dims.resize(4);
        conf.N = 4;
        return SZ_compress_N<T, 4>(conf, data, outSize);
    }
}

template<class T>
T *SZ_decompress(const SZ::Config &conf, char *cmpData, size_t cmpSize) {
    if (conf.N == 1) {
        return SZ_decompress_N<T, 1>(conf, cmpData, cmpSize);
    } else if (conf.N == 2) {
        return SZ_decompress_N<T, 2>(conf, cmpData, cmpSize);
    } else if (conf.N == 3) {
        return SZ_decompress_N<T, 3>(conf, cmpData, cmpSize);
    } else if (conf.N == 4) {
        return SZ_decompress_N<T, 4>(conf, cmpData, cmpSize);
    } else {
        return SZ_decompress_N<T, 4>(conf, cmpData, cmpSize);
    }
}