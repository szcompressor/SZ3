//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_HPP
#define SZ_LOSSLESS_HPP

#include "zstd.h"
#include "MemoryOps.hpp"

namespace SZ {
    uchar *lossless_compress(uchar *data, size_t dataLength, size_t &outSize) {
        size_t estimatedCompressedSize = dataLength < 100 ? 200 : dataLength * 1.2;
        uchar *compressBytes = new uchar[estimatedCompressedSize];
        uchar *compressBytesPos = compressBytes;
        write(dataLength, compressBytesPos);

        outSize = ZSTD_compress(compressBytesPos, estimatedCompressedSize, data, dataLength,
                                3); //default setting of level is 3
        outSize += sizeof(size_t);
        return compressBytes;
    }

    uchar *lossless_decompress(const uchar *data, size_t compressedSize) {
        const uchar *dataPos = data;
        size_t dataLength = 0, tmp = compressedSize;
        read(dataLength, dataPos, tmp);

        uchar *oriData = new uchar[dataLength];
        ZSTD_decompress(oriData, dataLength, dataPos, compressedSize);
        return oriData;
    }
}
#endif //SZ_LOSSLESS_HPP
