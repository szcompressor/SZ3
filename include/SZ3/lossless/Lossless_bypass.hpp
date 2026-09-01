//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ3_LOSSLESS_BYPASS_HPP
#define SZ3_LOSSLESS_BYPASS_HPP

#include <cstring>
#include <stdexcept>

#include "SZ3/def.hpp"
#include "SZ3/lossless/Lossless.hpp"

namespace SZ3 {
class Lossless_bypass : public concepts::LosslessInterface {
public:
    size_t compress(const uchar *src, size_t srcLen, uchar *dst, size_t dstCap) override {
        // Nothing here shrinks the payload, so a caller that sized its destination from a bound assuming
        // compression can be smaller than srcLen. Lossless_zstd already checks and throws; do the same
        // instead of memcpy'ing past the end of the destination.
        if (dstCap < srcLen) {
            throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        }
        std::memcpy(dst, src, srcLen);
        return srcLen;
    }

    size_t decompress(const uchar *src, const size_t srcLen, uchar *&dst, size_t &dstLen) override {
        // Mirror Lossless_zstd: when the caller asks us to allocate, a non-zero incoming dstLen is an
        // upper bound on the allocation; zero means no bound. dstLen stays a pure output parameter for a
        // caller-provided buffer.
        const size_t dst_capacity = dstLen;
        dstLen = srcLen;
        if (dst == nullptr) {
            if (dst_capacity != 0 && dstLen > dst_capacity) {
                throw std::out_of_range("SZ3 bypass lossless: payload exceeds the allowed capacity");
            }
            dst = static_cast<uchar *>(malloc(dstLen));
            if (dst == nullptr) {
                throw std::runtime_error("SZ3 bypass lossless: can not allocate the decompression buffer");
            }
        }
        std::memcpy(dst, src, dstLen);
        return dstLen;
    }
};
}  // namespace SZ3
#endif  // SZ_LOSSLESS_BYPASS_HPP
