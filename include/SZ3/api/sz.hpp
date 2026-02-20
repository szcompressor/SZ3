/**
 * @file sz.hpp
 * @ingroup API
 * @brief SZ3 compression and decompression API.
 *
 * This header provides the main API functions for compressing and decompressing data using SZ3.
 *
 * Compressed Data Format of SZ3:
 * The compressed data is always stored in little-endian order.
 * The compressed data consists of three main sections:
 * 1. Header (16 bytes): Contains metadata about the compressed data.
 *    - Magic Number (4 bytes): Identifies the data as SZ3-compressed.
 *    - Version (4 bytes): Indicates the version of the SZ3 format.
 *    - Compressed Size (8 bytes): Specifies the size of the compressed payload.
 * 2. Compressed Payload: The actual compressed data.
 * 3. Configuration: Stores the compression configuration used.
 *
 * The layout can be visualized as follows:
 * [ Header (16 bytes) | Compressed Payload | Configuration ]
 */

#ifndef SZ3_SZ_HPP
#define SZ3_SZ_HPP

#include "SZ3/api/impl/SZImpl.hpp"
#include "SZ3/version.hpp"


/**
 * @brief Compress data into a pre-allocated buffer.
 *
 * Compresses the input data using the provided configuration and stores the result in a pre-allocated buffer.
 * Use `SZ_compress_size_bound` to determine the minimum required buffer size.
 *
 * @tparam T The data type of the source data.
 * @param config The compression configuration (e.g., dims, error bound mode, algorithm).
 * @param data Pointer to the source data array.
 * @param cmpData Pointer to the pre-allocated buffer for compressed data.
 * @param cmpCap The capacity of the pre-allocated buffer in bytes.
 * @return The size of the compressed data in bytes.
 * @code
 * SZ3::Config conf(100, 200, 300); // 300 is the fastest dimension
 * conf.errorBoundMode = SZ3::EB_ABS;
 * conf.absErrorBound = 1E-3;
 * size_t outSize = SZ_compress(conf, data, outBuff, outBuffCap);
 * @endcode
 */
template <class T>
size_t SZ_compress(const SZ3::Config& config, const T* data, char* cmpData, size_t cmpCap) {
    using namespace SZ3;
    Config conf(config);

    if (cmpCap < SZ_compress_size_bound<T>(conf)) {
        throw std::invalid_argument(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
    }

    auto cmpDataPos = reinterpret_cast<uchar*>(cmpData);

    // save 16 bytes header
    write(config.sz3MagicNumber, cmpDataPos); // magic number (4 bytes)
    write(config.sz3DataVer, cmpDataPos);     // data version (4 bytes)
    auto sizeHeaderPos = cmpDataPos;
    cmpDataPos += 8; // reserve space for cmp data size (8 bytes)

    // begin compression
    auto cmpDataCap = cmpCap - 16 - conf.size_est() * 2;
    uint64_t cmpDataSize = 0;
    if (conf.N == 1) {
        cmpDataSize = SZ_compress_impl<T, 1>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 2) {
        cmpDataSize = SZ_compress_impl<T, 2>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 3) {
        cmpDataSize = SZ_compress_impl<T, 3>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 4) {
        cmpDataSize = SZ_compress_impl<T, 4>(conf, data, cmpDataPos, cmpDataCap);
    } else {
        throw std::invalid_argument("Data dimension higher than 4 is not supported.");
    }

    // save compressed data size back in header
    write(cmpDataSize, sizeHeaderPos);

    // save config
    cmpDataPos += cmpDataSize;
    auto confSize = conf.save(cmpDataPos);

    return 16 + cmpDataSize + confSize;
}

/**
 * @brief Compress data into a newly allocated buffer.
 *
 * Compresses the input data and returns a newly heap-allocated buffer.
 * The caller is responsible for freeing this buffer using `delete[]`.
 *
 * @tparam T The data type of the source data.
 * @param config The compression configuration.
 * @param data Pointer to the source data array.
 * @param cmpSize Output parameter set to the size of the compressed data in bytes.
 * @return Pointer to the newly allocated buffer. Must be freed by the caller using `delete[]`.
 */
template <class T>
char* SZ_compress(const SZ3::Config& config, const T* data, size_t& cmpSize) {
    using namespace SZ3;

    size_t bufferLen = SZ_compress_size_bound<T>(config);
    auto buffer = new char[bufferLen];
    cmpSize = SZ_compress(config, data, buffer, bufferLen);

    return buffer;
}

/**
 * @brief Decompress data into a pre-allocated (or auto-allocated) buffer.
 *
 * Reads the configuration from the compressed data stream, then decompresses into `decData`.
 * If `decData` is null, a new buffer will be allocated automatically.
 *
 * @tparam T The data type of the decompressed data.
 * @param config Output: overwritten with the configuration embedded in the compressed data.
 * @param cmpData Pointer to the compressed data.
 * @param cmpSize Size of the compressed data in bytes.
 * @param decData Pre-allocated output buffer, or null to auto-allocate.
 * @code
 * auto decData = new float[100 * 200 * 300];
 * SZ3::Config conf;
 * SZ_decompress(conf, cmpData, cmpSize, decData);
 * @endcode
 */
template <class T>
void SZ_decompress(SZ3::Config& config, const char* cmpData, size_t cmpSize, T*& decData) {
    using namespace SZ3;

    auto cmpDataPos = reinterpret_cast<const uchar*>(cmpData);

    read(config.sz3MagicNumber, cmpDataPos);
    if (config.sz3MagicNumber != SZ3_MAGIC_NUMBER) {
        throw std::invalid_argument("magic number mismatch, the input data is not compressed by SZ3");
    }

    read(config.sz3DataVer, cmpDataPos);
    if (versionStr(config.sz3DataVer) != SZ3_DATA_VER) {
        std::stringstream ss;
        printf("program v%s , program-data %s , input data v%s\n", SZ3_VER, SZ3_DATA_VER,
               versionStr(config.sz3DataVer).data());
        ss << "Please use SZ3 v" << versionStr(config.sz3DataVer) << " to decompress the data" << std::endl;
        std::cerr << ss.str();
        throw std::invalid_argument(ss.str());
    }

    uint64_t cmpDataSize = 0;
    read(cmpDataSize, cmpDataPos);

    auto cmpConfPos = cmpDataPos + cmpDataSize;
    config.load(cmpConfPos);

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
        throw std::invalid_argument("Data dimension higher than 4 is not supported.");
    }
}

/**
 * @brief Decompress data into a newly allocated buffer.
 *
 * Convenience overload that auto-allocates the output buffer and returns it.
 * The caller is responsible for freeing the returned buffer with `delete[]`.
 *
 * @tparam T The data type of the decompressed data.
 * @param config Output: overwritten with the configuration embedded in the compressed data.
 * @param cmpData Pointer to the compressed data.
 * @param cmpSize Size of the compressed data in bytes.
 * @return Pointer to the newly allocated decompressed data. Must be freed with `delete[]`.
 */
template <class T>
T* SZ_decompress(SZ3::Config& config, const char* cmpData, size_t cmpSize) {
    using namespace SZ3;
    T* decData = nullptr;
    SZ_decompress<T>(config, cmpData, cmpSize, decData);
    return decData;
}

#endif
