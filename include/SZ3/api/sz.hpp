#ifndef SZ3_SZ_HPP
#define SZ3_SZ_HPP

#include "SZ3/api/impl/SZImpl.hpp"
#include "SZ3/version.hpp"

/**
 * API for compression
 * @tparam T source data type
 * @param config compression configuration. Please update the config with 1). data dimension and shape and 2). desired
settings.
 * @param data source data
 * @param cmpData pre-allocated buffer for compressed data
 * @param cmpCap pre-allocated buffer size (in bytes) for compressed data
 * @return compressed data size (in bytes)

The compression algorithms are:
ALGO_INTERP_LORENZO:
 The default algorithm in SZ3. It is the implementation of our ICDE'21 paper.
 The whole dataset will be compressed by interpolation or lorenzo predictor with auto-optimized settings.
ALGO_INTERP:
 The whole dataset will be compressed by interpolation predictor with default settings.
ALGO_LORENZO_REG:
 The whole dataset will be compressed by lorenzo and/or regression based predictors block by block with default
settings. The four predictors ( 1st-order lorenzo, 2nd-order lorenzo, 1st-order regression, 2nd-order regression) can be
enabled or disabled independently by conf settings (lorenzo, lorenzo2, regression, regression2).

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
template <class T>
size_t SZ_compress(const SZ3::Config &config, const T *data, char *cmpData, size_t cmpCap) {
    using namespace SZ3;
    Config conf(config);
    
    if (cmpCap < SZ_compress_size_bound<T>(conf)) {
        fprintf(stderr, "%s\n", SZ_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        throw std::invalid_argument(SZ_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
    }
    
    auto cmpDataPos = reinterpret_cast<uchar *>(cmpData) + conf.size_est();
    auto cmpDataCap = cmpCap - conf.size_est();

    size_t cmpDataLen = 0;
    if (conf.N == 1) {
        cmpDataLen = SZ_compress_impl<T, 1>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 2) {
        cmpDataLen = SZ_compress_impl<T, 2>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 3) {
        cmpDataLen = SZ_compress_impl<T, 3>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 4) {
        cmpDataLen = SZ_compress_impl<T, 4>(conf, data, cmpDataPos, cmpDataCap);
    } else {
        fprintf(stderr, "Data dimension higher than 4 is not supported.\n");
        throw std::invalid_argument("Data dimension higher than 4 is not supported.");
    }

    auto cmpConfPos = reinterpret_cast<uchar *>(cmpData);
    conf.save(cmpConfPos);
    return conf.size_est() + cmpDataLen;
}

/**
 * API for compression
 * @tparam T  source data type
 * @param config config compression configuration
 * @param data source data
 * @param cmpSize compressed data size (in bytes)
 * @return compressed data, remember to 'delete []' when the data is no longer needed.
 *
 * Similar with SZ_compress(SZ3::Config &conf, const T *data, char *cmpData, size_t cmpCap)
 * The only difference is this one doesn't need the pre-allocated buffer (thus remember to do 'delete []' yourself)
 */
template <class T>
char *SZ_compress(const SZ3::Config &config, const T *data, size_t &cmpSize) {
    using namespace SZ3;

    size_t bufferLen = SZ_compress_size_bound<T>(config);
    auto buffer = new char[bufferLen];
    cmpSize = SZ_compress(config, data, buffer, bufferLen);

    return buffer;
}

/**
 * API for decompression
 * @tparam T decompressed data type
 * @param config configuration placeholder. It will be overwritten by the compression configuration
 * @param cmpData compressed data
 * @param cmpSize compressed data size in bytes
 * @param decData pre-allocated buffer for decompressed data

 example:
 auto decData = new float[100*200*300];
 SZ3::Config conf;
 SZ_decompress(conf, cmpData, cmpSize, decData);

 */
template <class T>
void SZ_decompress(SZ3::Config &config, const char *cmpData, size_t cmpSize, T *&decData) {
    using namespace SZ3;
    auto cmpConfPos = reinterpret_cast<const uchar *>(cmpData);
    config.load(cmpConfPos);

    auto cmpDataPos = reinterpret_cast<const uchar *>(cmpData) + config.size_est();
    auto cmpDataSize = cmpSize - config.size_est();

    if (decData == nullptr) {
        decData = new T[config.num];
    }
    if (config.N == 1) {
        SZ_decompress_impl<T, 1>(config, cmpDataPos, cmpDataSize, decData);
    } else if (config.N == 2) {
        SZ_decompress_impl<T, 2>(config, cmpDataPos, cmpDataSize, decData);
    } else if (config.N == 3) {
        SZ_decompress_impl<T, 3>(config, cmpDataPos, cmpDataSize, decData);
    } else if (config.N == 4) {
        SZ_decompress_impl<T, 4>(config, cmpDataPos, cmpDataSize, decData);
    } else {
        fprintf(stderr, "Data dimension higher than 4 is not supported.\n");
        throw std::invalid_argument("Data dimension higher than 4 is not supported.");
    }
}

/**
 * API for decompression
 * Similar with SZ_decompress(SZ3::Config &config, char *cmpData, size_t cmpSize, T *&decData)
 * The only difference is this one doesn't need pre-allocated buffer for decompressed data
 *
 * @tparam T decompressed data type
 * @param config configuration placeholder. It will be overwritten by the compression configuration
 * @param cmpData compressed data
 * @param cmpSize compressed data size in bytes
 * @return decompressed data, remember to 'delete []' when the data is no longer needed.

 example:
 SZ3::Config conf;
 float decompressedData = SZ_decompress(conf, cmpData, cmpSize)
 */
template <class T>
T *SZ_decompress(SZ3::Config &config, const char *cmpData, size_t cmpSize) {
    using namespace SZ3;
    T *decData = nullptr;
    SZ_decompress<T>(config, cmpData, cmpSize, decData);
    return decData;
}

#endif
