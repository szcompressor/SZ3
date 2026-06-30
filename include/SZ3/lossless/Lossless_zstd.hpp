//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ3_LOSSLESS_ZSTD_HPP
#define SZ3_LOSSLESS_ZSTD_HPP

#include <stdexcept>

#include "SZ3/def.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "zstd.h"

namespace SZ3 {
class Lossless_zstd : public concepts::LosslessInterface {
   public:
    Lossless_zstd() = default;

    Lossless_zstd(int comp_level) : compression_level(comp_level) {}

    /**
     * Attention
     * When dstCap is smaller than the space needed, ZSTD will not throw any errors.
     * Instead, it will write a portion of the compressed data to dst and stops.
     * This behavior is not desirable in SZ, as we need the whole compressed data for decompression.
     * Therefore, we need to check if the dst buffer (dstCap) is large enough for zstd
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
        /// The compressed buffer starts with the decompressed size, followed by the zstd stream.
        /// Validate the inputs so corrupted data can not read out of bounds, allocate an untrusted
        /// amount, or pass a zstd error code back as a size.
        if (srcLen < sizeof(dstLen)) {
            throw std::out_of_range("SZ3 lossless: compressed data is smaller than the size header");
        }
        /// When the caller provides the destination buffer, the incoming dstLen is its capacity. The
        /// decompressed size is read from the (untrusted) compressed data, so it must not exceed that
        /// capacity, otherwise zstd would write past the end of the buffer.
        const size_t dst_capacity = dstLen;
        read(dstLen, src);
        if (dst == nullptr) {
            /// dst == nullptr means the caller asks us to allocate the output buffer, and its size is read
            /// from the (untrusted) compressed payload. When the caller supplies a non-zero capacity it is an
            /// upper bound on how large that allocation may be; reject a payload that declares a larger size
            /// before allocating it, so corrupted data can not force an arbitrary (and untracked) allocation.
            /// A zero capacity means the caller did not supply a bound (legacy behavior).
            if (dst_capacity != 0 && dstLen > dst_capacity) {
                throw std::out_of_range("SZ3 lossless: declared decompressed size exceeds the allowed capacity");
            }
            dst = static_cast<uchar *>(malloc(dstLen));
            if (dst == nullptr) {
                throw std::runtime_error("SZ3 lossless: can not allocate the decompression buffer");
            }
        } else if (dstLen > dst_capacity) {
            throw std::out_of_range("SZ3 lossless: decompressed size exceeds the destination buffer");
        }
        size_t res = ZSTD_decompress(dst, dstLen, src, srcLen - sizeof(dstLen));
        if (ZSTD_isError(res)) {
            throw std::runtime_error("SZ3 lossless: zstd decompression failed");
        }
        /// The declared size is read from the (untrusted) payload; require zstd to actually produce that many
        /// bytes, so a frame that expands to fewer bytes can not leave the tail of the output buffer
        /// uninitialized (which a caller would otherwise copy out as if it were decompressed data).
        if (res != dstLen) {
            throw std::out_of_range("SZ3 lossless: decompressed size does not match the declared size");
        }
        return res;
    }

   private:
    int compression_level = 3;  // default setting of level is 3
};
}  // namespace SZ3
#endif  // SZ_LOSSLESS_ZSTD_HPP
