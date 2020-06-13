//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_TORNADO_HPP
#define SZ_LOSSLESS_TORNADO_HPP

#ifdef ENABLE_LZBENCH

#include "tornado/tor_test.h"

#include "def.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/fileUtil.h"
#include "Lossless.hpp"

namespace SZ {
    class Lossless_TORNADO : public concepts::LosslessInterface {

    public:
        uchar *compress(uchar *data, size_t dataLength, size_t &outSize) {
            size_t estimatedCompressedSize = dataLength < 100 ? 200 : dataLength * 1.2;
            uchar *compressBytes = new uchar[estimatedCompressedSize];
            uchar *compressBytesPos = compressBytes;
            write(dataLength, compressBytesPos);

            outSize = tor_compress(3, data, dataLength, compressBytesPos, estimatedCompressedSize);

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
            tor_decompress((uchar *) dataPos, compressedSize, oriData, dataLength);
            return oriData;
        }

        void postdecompress_data(uchar *data) {
            delete[] data;
        }

    };
}
#endif //ENABLE_LZBENCH
#endif //SZ_LOSSLESS_TORNADO_HPP