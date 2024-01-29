//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_ZSTD_HPP
#define SZ_LOSSLESS_ZSTD_HPP

#include "zstd.h"
#include "SZ3/def.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/lossless/Lossless.hpp"

namespace SZ3 {
    class Lossless_zstd : public concepts::LosslessInterface {

    public:
        Lossless_zstd() = default;

        Lossless_zstd(int comp_level) : compression_level(comp_level) {};

        void compress(uchar *src, size_t srcLen, uchar *dst, size_t &dstCap) {
//            size_t estimatedCompressedSize = std::max(size_t(srcLen * 1.2), size_t(400));
//            uchar *compressBytes = new uchar[estimatedCompressedSize];
//            uchar *dstPos = dst;
//            write(srcLen, dstPos);
            if (dstCap < srcLen) {
                throw std::invalid_argument(
                        "dstCap/dstLen too small, remember to initialize the dstCap/dstLen with the array size when you call compression func");
            }
            dstCap = ZSTD_compress(dst, dstCap, src, srcLen, compression_level);
//            dstLen += sizeof(size_t);
//            return compressBytes;
        }


        void decompress(const uchar *src, const size_t srcLen, uchar *dst, size_t &dstCap) {
//            const uchar *dataPos = data;
//            size_t dataLength = 0;
//            read(dataLength, dataPos, compressedSize);

//            uchar *oriData = new uchar[dataLength];
            dstCap = ZSTD_decompress(dst, dstCap, src, srcLen);
//            compressedSize = dataLength;
//            return oriData;
        }

    private:
        int compression_level = 3;  //default setting of level is 3
    };
}
#endif //SZ_LOSSLESS_ZSTD_HPP
