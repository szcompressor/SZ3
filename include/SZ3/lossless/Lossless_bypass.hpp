//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_BYPASS_HPP
#define SZ_LOSSLESS_BYPASS_HPP

#include "zstd.h"
#include "SZ3/def.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/lossless/Lossless.hpp"

namespace SZ3 {
    class Lossless_bypass : public concepts::LosslessInterface {
     
     public:
        
        size_t compress(uchar *src, size_t srcLen, uchar *dst, size_t dstCap) {
            dst = src;
            return srcLen;
        }
        
        size_t decompress(const uchar *src, const size_t srcLen, uchar *dst, size_t dstCap) {
            dst = (uchar *) src;
            return srcLen;
        }
    };
}
#endif //SZ_LOSSLESS_BYPASS_HPP
