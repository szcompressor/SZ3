//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_LZFSE_HPP
#define SZ_LOSSLESS_LZFSE_HPP
#ifdef ENABLE_LZBENCH

#include "lzfse/lzfse.h"

#include "def.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/fileUtil.h"
#include "lossless/Lossless.hpp"

namespace SZ {
    class Lossless_LZFSE : public concepts::LosslessInterface {

    public:
        uchar *compress(uchar *data, size_t dataLength, size_t &outSize) {
            size_t estimatedCompressedSize = dataLength < 100 ? 200 : dataLength * 1.2;
            uchar *compressBytes = new uchar[estimatedCompressedSize];
            uchar *compressBytesPos = compressBytes;
            write(dataLength, compressBytesPos);

            outSize = lzfse_encode_buffer(compressBytesPos, estimatedCompressedSize, data, dataLength, NULL);

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
            lzfse_decode_buffer(oriData, dataLength, dataPos, compressedSize, NULL);
            return oriData;
        }

        void postdecompress_data(uchar *data) {
            delete[] data;
        }

    };
}
#endif //ENABLE_LZBENCH
#endif //SZ_LOSSLESS_LZFSE_HPP