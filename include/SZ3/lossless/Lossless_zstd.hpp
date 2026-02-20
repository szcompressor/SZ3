/**
 * @file Lossless_zstd.hpp
 * @ingroup Lossless
 */

#ifndef SZ3_LOSSLESS_ZSTD_HPP
#define SZ3_LOSSLESS_ZSTD_HPP

#include <stdexcept>

#include "SZ3/def.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "zstd.h"

namespace SZ3 {

/**
 * @brief Zstd Lossless Compressor wrapper
 * 
 * Uses the Zstd library for lossless compression.
 */
class Lossless_zstd : public concepts::LosslessInterface {
   public:
    Lossless_zstd() = default;

    /**
     * @brief Construct a new Zstd Lossless object
     * 
     * @param comp_level Zstd compression level (default is 3)
     */
    Lossless_zstd(int comp_level) : compression_level(comp_level) {}

    /**
     * @brief Compress data using Zstd
     * 
     * Note: Checks if the destination buffer is large enough using ZSTD_compressBound.
     * Throws an error if insufficient.
     * Zstd itself will not throw error and instead will write a portion of data when the destination buffer is not large enough.
     * 
     * @param src Input data
     * @param srcLen Input length
     * @param dst Output buffer
     * @param dstCap Output capacity
     * @return size_t Compressed size
     */
    size_t compress(const uchar *src, size_t srcLen, uchar *dst, size_t dstCap) override {
        write(srcLen, dst);
        dstCap -= sizeof(size_t);  // reserve space for srcLen
        if (dstCap < ZSTD_compressBound(srcLen)) {
            throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        }
        size_t dstLen = ZSTD_compress(dst, dstCap, src, srcLen, compression_level);
        return dstLen + sizeof(size_t);
    }

    size_t decompress(const uchar *src, const size_t srcLen, uchar *&dst, size_t &dstLen) override {
        read(dstLen, src);
        if (dst == nullptr) {
            dst = static_cast<uchar *>(malloc(dstLen));
        }
        return ZSTD_decompress(dst, dstLen, src, srcLen - sizeof(dstLen));
    }

   private:
    int compression_level = 3;  // default setting of level is 3
};
}  // namespace SZ3
#endif  // SZ_LOSSLESS_ZSTD_HPP
