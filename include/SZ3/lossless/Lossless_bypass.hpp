/**
 * @file Lossless_bypass.hpp
 * @ingroup Lossless
 */

#ifndef SZ3_LOSSLESS_BYPASS_HPP
#define SZ3_LOSSLESS_BYPASS_HPP

#include <cstdlib>
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
        dstLen = srcLen;
        if (dst == nullptr) {
            dst = static_cast<uchar *>(malloc(dstLen));
        }
        std::memcpy(dst, src, dstLen);
        return dstLen;
    }
};
}  // namespace SZ3
#endif  // SZ_LOSSLESS_BYPASS_HPP
