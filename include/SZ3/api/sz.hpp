#ifndef SZ3_SZ_HPP
#define SZ3_SZ_HPP

#include "SZ3/api/impl/SZImpl.hpp"
#include "SZ3/version.hpp"
#include <memory>

/**
 * API for compression
 * @tparam T source data type
 * @param config compression configuration. Please update the config with 1). data dimension and shape and 2). desired settings.
 * @param data source data
 * @param cmpSize compressed data size in bytes
 * @return compressed data, remember to 'delete []' when the data is no longer needed.

The compression algorithms are:
ALGO_INTERP_LORENZO:
 The default algorithm in SZ3. It is the implementation of our ICDE'21 paper.
 The whole dataset will be compressed by interpolation or lorenzo predictor with auto-optimized settings.
ALGO_INTERP:
 The whole dataset will be compressed by interpolation predictor with default settings.
ALGO_LORENZO_REG:
 The whole dataset will be compressed by lorenzo and/or regression based predictors block by block with default settings.
 The four predictors ( 1st-order lorenzo, 2nd-order lorenzo, 1st-order regression, 2nd-order regression)
 can be enabled or disabled independently by conf settings (lorenzo, lorenzo2, regression, regression2).

Interpolation+lorenzo example:
SZ3::Config conf(100, 200, 300); // 300 is the fastest dimension
conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
conf.errorBoundMode = SZ3::EB_ABS; // refer to def.hpp for all supported error bound mode
conf.absErrorBound = 1E-3; // absolute error bound 1e-3
char *compressedData = SZ_compress(conf, data, outSize);

Interpolation example:
SZ3::Config conf(100, 200, 300); // 300 is the fastest dimension
conf.cmprAlgo = SZ3::ALGO_INTERP;
conf.errorBoundMode = SZ3::EB_REL; // refer to def.hpp for all supported error bound mode
conf.relErrorBound = 1E-3; // value-rang-based error bound 1e-3
char *compressedData = SZ_compress(conf, data, outSize);

Lorenzo/regression example :
SZ3::Config conf(100, 200, 300); // 300 is the fastest dimension
conf.cmprAlgo = SZ3::ALGO_LORENZO_REG;
conf.lorenzo = true; // only use 1st order lorenzo
conf.lorenzo2 = false;
conf.regression = false;
conf.regression2 = false;
conf.errorBoundMode = SZ3::EB_ABS; // refer to def.hpp for all supported error bound mode
conf.absErrorBound = 1E-3; // absolute error bound 1e-3
char *compressedData = SZ_compress(conf, data, outSize);
 */
template<class T>
size_t SZ_compress(const SZ3::Config &conf_, const T *data, char *cmpData, size_t cmpCap) {
    using namespace SZ3;
    Config conf(conf_);
    
    if (cmpCap < conf.num * sizeof(T)) {
        throw std::invalid_argument(
            "cmpCap too small, remember to initialize the cmpCap with at least the same size of the original data");
    }
    
    auto dst = (uchar*)cmpData + conf.size_est();
    auto dstCap = cmpCap - conf.size_est();
    
    size_t dstLen = 0;
    if (conf.N == 1) {
        dstLen = SZ_compress_impl<T, 1>(conf, data, dst, dstCap);
    } else if (conf.N == 2) {
        dstLen = SZ_compress_impl<T, 2>(conf, data, dst, dstCap);
    } else if (conf.N == 3) {
        dstLen = SZ_compress_impl<T, 3>(conf, data, dst, dstCap);
    } else if (conf.N == 4) {
        dstLen = SZ_compress_impl<T, 4>(conf, data, dst, dstCap);
    } else {
        printf("Data dimension higher than 4 is not supported.\n");
        exit(0);
    }
    
    
    auto confPos = (uchar *) cmpData;
    conf.save(confPos);
    return conf.size_est() + dstLen;
}

template<class T>
char *SZ_compress(const SZ3::Config &conf, const T *data, size_t &cmpSize) {
    using namespace SZ3;
    
    size_t bufferLen = conf.num * sizeof(T) * 1.2;
    auto buffer = new char[bufferLen];
    cmpSize = SZ_compress(conf, data, buffer, bufferLen);
    
    return buffer;
}


/**
 * API for decompression
 * Similar with SZ_decompress(SZ3::Config &conf, char *cmpData, size_t cmpSize)
 * The only difference is this one needs pre-allocated decData as input
 * @tparam T decompressed data type
 * @param conf configuration placeholder. It will be overwritten by the compression configuration
 * @param cmpData compressed data
 * @param cmpSize compressed data size in bytes
 * @param decData pre-allocated memory space for decompressed data

 example:
 auto decData = new float[100*200*300];
 SZ3::Config conf;
 SZ_decompress(conf, cmpData, cmpSize, decData);

 */
template<class T>
void SZ_decompress(SZ3::Config &conf, char *cmpData, size_t cmpSize, T *&decData) {
    using namespace SZ3;
    auto confPos = (const uchar *) cmpData;
    auto cmpDataPos = confPos + conf.size_est();
    conf.load(confPos);
    
    if (decData == nullptr) {
        decData = new T[conf.num];
    }
    if (conf.N == 1) {
        SZ_decompress_impl<T, 1>(conf, cmpDataPos, cmpSize, decData);
    } else if (conf.N == 2) {
        SZ_decompress_impl<T, 2>(conf, cmpDataPos, cmpSize, decData);
    } else if (conf.N == 3) {
        SZ_decompress_impl<T, 3>(conf, cmpDataPos, cmpSize, decData);
    } else if (conf.N == 4) {
        SZ_decompress_impl<T, 4>(conf, cmpDataPos, cmpSize, decData);
    } else {
        printf("Data dimension higher than 4 is not supported.\n");
        exit(0);
    }
}

/**
 * API for decompression
 * @tparam T decompressed data type
 * @param conf configuration placeholder. It will be overwritten by the compression configuration
 * @param cmpData compressed data
 * @param cmpSize compressed data size in bytes
 * @return decompressed data, remember to 'delete []' when the data is no longer needed.

 example:
 SZ3::Config conf;
 float decompressedData = SZ_decompress(conf, cmpData, cmpSize)
 */
template<class T>
T *SZ_decompress(SZ3::Config &conf, char *cmpData, size_t cmpSize) {
    using namespace SZ3;
    T *decData = nullptr;
    SZ_decompress<T>(conf, cmpData, cmpSize, decData);
    return decData;
}

#endif