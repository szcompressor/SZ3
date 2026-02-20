/**
 * @file Lossless.hpp
 */

#ifndef SZ3_LOSSLESS_HPP
#define SZ3_LOSSLESS_HPP

#include "SZ3/def.hpp"

namespace SZ3::concepts {

/**
 * @brief Interface for Lossless compressors
 * 
 * Lossless compression is applied after the lossy stage to further reduce data size
 * without losing any additional information. Typically wraps existing libraries (e.g., Zstd, Gzip).
 */
class LosslessInterface {
   public:
    virtual ~LosslessInterface() = default;

    /**
     * @brief Compress data using a lossless algorithm
     * 
     * @param src Input data to be compressed
     * @param srcLen Length of input data in bytes
     * @param dst Buffer to store compressed data
     * @param dstCap Capacity of the destination buffer
     * @return size_t Size of compressed data in bytes
     */
    virtual size_t compress(const uchar *src, size_t srcLen, uchar *dst, size_t dstCap) = 0;

    /**
     * @brief Decompress data using a lossless algorithm
     * 
     * @param src Compressed data buffer
     * @param srcLen Length of compressed data
     * @param dst Buffer to store decompressed data (reference to pointer)
     * @param dstLen Length of decompressed data (reference, updated by function)
     * @return size_t Size of decompressed data in bytes
     */
    virtual size_t decompress(const uchar *src, const size_t srcLen, uchar *&dst, size_t &dstLen) = 0;
};
}  // namespace SZ3::concepts

#endif  // SZ_LOSSLESS_HPP
