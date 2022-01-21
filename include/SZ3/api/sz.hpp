#ifndef SZ3_SZ_HPP
#define SZ3_SZ_HPP

#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZFastFrontend.hpp"
#include "SZ3/frontend/SZGeneralFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/predictor/PolyRegressionPredictor.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/api/impl/szn.hpp"
#include <cmath>
#include <memory>

/**
 * API for compression
 * @tparam T source data type
 * @param conf compression configuration. Please update the config with 1). data dimension and shape and 2). desired settings.
 * @param data source data
 * @param outSize compressed data size in bytes
 * @return compressed data
The compression algorithms are:
METHOD_INTERP_LORENZO:
 The default algorithm in SZ3. It is the implementation of our ICDE'21 paper.
 The whole dataset will be compressed by interpolation or lorenzo predictor with *optimized* settings.
METHOD_INTERP:
 The whole dataset will be compressed by interpolation predictor with default settings.
METHOD_LORENZO_REG:
 The whole dataset will be compressed by lorenzo/regression based predictors block by block with default settings.
 The four predictors ( 1st-order lorenzo, 2nd-order lorenzo, 1st-order regression, 2nd-order regression)
 can be enabled or disabled independently by conf settings (enable_lorenzo, enable_2ndlorenzo, enable_regression, enable_2ndregression).

example:
SZ::Config conf(100, 200, 300); // 300 is the fastest dimension
conf.cmprMethod = METHOD_INTERP_LORENZO; // use interp+lorenzo as the algorithm
conf.errorBoundMode = ABS; // refer to def.hpp for all supported error bound mode
conf.absErrorBound = 1E-3; // absolute error bound 1e-3
char *compressedData = SZ_compress(conf, data, outSize);
 */
template<class T>
char *SZ_compress(SZ::Config &conf, T *data, size_t &outSize) {
    char *cmpData;
    if (conf.N == 1) {
        cmpData = SZ_compress_N<T, 1>(conf, data, outSize);
    } else if (conf.N == 2) {
        cmpData = SZ_compress_N<T, 2>(conf, data, outSize);
    } else if (conf.N == 3) {
        cmpData = SZ_compress_N<T, 3>(conf, data, outSize);
    } else if (conf.N == 4) {
        cmpData = SZ_compress_N<T, 4>(conf, data, outSize);
    } else {
        for (int i = 4; i < conf.N; i++) {
            conf.dims[3] *= conf.dims[i];
        }
        conf.dims.resize(4);
        conf.N = 4;
        cmpData = SZ_compress_N<T, 4>(conf, data, outSize);
    }
    {
        //save config
        SZ::uchar *cmpDataPos = (SZ::uchar *) cmpData + outSize;
        conf.save(cmpDataPos);
        size_t newSize = (char *) cmpDataPos - cmpData;
        SZ::write(int(newSize - outSize), cmpDataPos);
        outSize = (char *) cmpDataPos - cmpData;
    }
    return cmpData;
}

/**
 * API for decompression
 * @tparam T decompressed data type
 * @param conf compression configuration. Setting the correct config is NOT needed for decompression.
 * The correct config will be loaded from compressed data and returned.
 * @param cmpData compressed data
 * @param cmpSize compressed data size in bytes
 * @return decompressed data
 example:
 SZ::Config conf;
 float decompressedData = SZ_decompress(conf, char *cmpData, size_t cmpSize)
 */
template<class T>
T *SZ_decompress(SZ::Config &conf, char *cmpData, size_t cmpSize) {
    {
        //load config
        int confSize;
        memcpy(&confSize, cmpData + (cmpSize - sizeof(int)), sizeof(int));
        SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData + (cmpSize - sizeof(int) - confSize);
        conf.load(cmpDataPos);
    }

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


#endif