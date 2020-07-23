//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_LZ4_HPP
#define SZ_LOSSLESS_LZ4_HPP

#include "lz4/lz4.h"
#include "def.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/fileUtil.h"
#include "lossless/Lossless.hpp"

namespace SZ {
    class Lossless_LZ4 : public concepts::LosslessInterface {

    public:
        uchar *compress(uchar *data, size_t dataLength, size_t &outSize) {
            size_t estimatedCompressedSize = dataLength < 100 ? 200 : dataLength * 1.2;
            uchar *compressBytes = new uchar[estimatedCompressedSize];
            uchar *compressBytesPos = compressBytes;
            write(dataLength, compressBytesPos);

            LZ4_compress_default(data, compressBytesPos, dataLength, estimatedCompressedSize);
            outSize += sizeof(size_t);
            return compressBytes;
        }

        void postcompress_data(uchar *data) {
            delete[] data;
        }

        uchar *decompress(const uchar *data, size_t compressedSize) {
            const uchar *dataPos = data;
            size_t dataLength = 0, tmp = compressedSize;
            read(dataLength, dataPos, tmp);

            uchar *oriData = new uchar[dataLength];
            LZ4_decompress_safe(dataPos, oriData, compressedSize, dataLength);
            return oriData;
        }

        void postdecompress_data(uchar *data) {
            delete[] data;
        }

    };
}
#endif //SZ_LOSSLESS_LZ4_HPP