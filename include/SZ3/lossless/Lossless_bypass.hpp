//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_BYPASS_HPP
#define SZ_LOSSLESS_BYPASS_HPP

#include "SZ3/def.hpp"
#include "SZ3/lossless/Lossless.hpp"

namespace SZ3 {
class Lossless_bypass : public concepts::LosslessInterface {
   public:
    size_t compress(const uchar *src, size_t srcLen, uchar *dst, size_t dstCap) override {
        std::memcpy(dst, src, srcLen);
        // dst = src;
        return srcLen;
    }

    size_t decompress(const uchar *src, const size_t srcLen, uchar *&dst, size_t &dstLen) override {
        std::memcpy(dst, src, srcLen);
        // dst = (uchar *)src;
        return srcLen;
    }
};
}  // namespace SZ3
#endif  // SZ_LOSSLESS_BYPASS_HPP
