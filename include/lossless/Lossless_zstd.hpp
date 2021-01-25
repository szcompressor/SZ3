//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_ZSTD_HPP
#define SZ_LOSSLESS_ZSTD_HPP

#include "zstd.h"
#include "def.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/FileUtil.h"
#include "lossless/Lossless.hpp"

namespace SZ {
    class Lossless_zstd : public concepts::LosslessInterface {

    public:
        uchar *compress(uchar *data, size_t dataLength, size_t &outSize) {
            size_t estimatedCompressedSize = dataLength < 100 ? 200 : dataLength * 1.2;
            uchar *compressBytes = new uchar[estimatedCompressedSize];
            uchar *compressBytesPos = compressBytes;
            write(dataLength, compressBytesPos);

            outSize = ZSTD_compress(compressBytesPos, estimatedCompressedSize, data, dataLength,
                                    3); //default setting of level is 3
            outSize += sizeof(size_t);
            return compressBytes;
        }

        void postcompress_data(uchar *data) {
            delete[] data;
        }

        uchar *decompress(const uchar *data, size_t &compressedSize) {
            const uchar *dataPos = data;
            size_t dataLength = 0;
            read(dataLength, dataPos, compressedSize);

            uchar *oriData = new uchar[dataLength];
            ZSTD_decompress(oriData, dataLength, dataPos, compressedSize);
            compressedSize = dataLength;
            return oriData;
        }

        void postdecompress_data(uchar *data) {
            delete[] data;
        }

    };
}
#endif //SZ_LOSSLESS_ZSTD_HPP