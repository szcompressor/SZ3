//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_LIZARD_HPP
#define SZ_LOSSLESS_LIZARD_HPP

#ifdef LZBENCH
#include "lizard/lizard_compress.h"
#include "lizard/lizard_decompress.h"
#include "def.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/fileUtil.h"
#include "Lossless.hpp"

namespace SZ {
    class Lossless_LIZARD : public concepts::LosslessInterface {

    public:
        uchar *compress(uchar *data, size_t dataLength, size_t &outSize) {
            size_t estimatedCompressedSize = dataLength < 100 ? 200 : dataLength * 1.2;
            uchar *compressBytes = new uchar[estimatedCompressedSize];
            uchar *compressBytesPos = compressBytes;
            write(dataLength, compressBytesPos);

            outSize = Lizard_compress(reinterpret_cast<const char *>(data), reinterpret_cast<char *>(compressBytesPos),
                                      dataLength, estimatedCompressedSize, 12); //default setting of level is 3
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
            Lizard_decompress_safe(reinterpret_cast<const char *>(dataPos), reinterpret_cast<char *>(oriData), compressedSize,
                                   dataLength);
            return oriData;
        }

        void postdecompress_data(uchar *data) {
            delete[] data;
        }

    };
}
#endif //LZBENCH
#endif //SZ_LOSSLESS_LIZARD_HPP