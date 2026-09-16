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
        if (dstCap < srcLen) {
            throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        }
        std::memcpy(dst, src, srcLen);
        return srcLen;
    }

    size_t decompress(const uchar *src, const size_t srcLen, uchar *&dst, size_t &dstLen) override {
        // dstLen caps a self-allocation only, as in Lossless_zstd.
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
