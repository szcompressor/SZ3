//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_ZSTD_HPP
#define SZ_LOSSLESS_ZSTD_HPP

#include "zstd.h"
#include "def.hpp"
#include "utils/MemoryUtil.hpp"
#include "utils/FileUtil.hpp"
#include "lossless/Lossless.hpp"

namespace SZ {
    class Lossless_zstd : public concepts::LosslessInterface {

    public:
        Lossless_zstd() = default;

        Lossless_zstd(int comp_level) : compression_level(comp_level) {};

        uchar *compress(uchar *dataIn, size_t inSize, size_t &outSize) {
            size_t estimatedOutSize = inSize < 100 ? 200 : inSize * 1.2;
            uchar *buffer = new uchar[estimatedOutSize];
            outSize = compress(dataIn, inSize, buffer);
            return buffer;
        }

        size_t compress(uchar *dataIn, size_t inSize, uchar *dataOut) {
            size_t estimatedCompressedSize = inSize < 100 ? 200 : inSize * 1.2;
            uchar *dataOutPos = dataOut;
            write(inSize, dataOutPos);

            size_t outSize = ZSTD_compress(dataOutPos, estimatedCompressedSize,
                                           dataIn, inSize, compression_level);
            outSize += sizeof(size_t);
//            printf("[ZSTD] ratio = %.2f inSize = %lu outSize = %lu\n", inSize * 1.0 / outSize, inSize, outSize);
            return outSize;
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

    private:
        int compression_level = 3;  //default setting of level is 3
    };
}
#endif //SZ_LOSSLESS_ZSTD_HPP