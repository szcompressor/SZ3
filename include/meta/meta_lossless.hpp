#ifndef _meta_lossless_hpp
#define _meta_lossless_hpp

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "zstd.h"
#include "meta_utils.hpp"

namespace META {

#define ZSTD_COMPRESSOR 1

    unsigned long meta_lossless_compress(int losslessCompressor, int level, unsigned char *data, unsigned long dataLength,
                                       unsigned char **compressBytes) {
        unsigned long outSize = 0;
        size_t estimatedCompressedSize = 0;
        switch (losslessCompressor) {
            case ZSTD_COMPRESSOR:
                if (dataLength < 100)
                    estimatedCompressedSize = 200;
                else
                    estimatedCompressedSize = dataLength * 1.2;
                *compressBytes = (unsigned char *) malloc(estimatedCompressedSize);
                outSize = ZSTD_compress(*compressBytes, estimatedCompressedSize, data, dataLength,
                                        level); //default setting of level is 3
                break;
            default:
                printf("Error: Unrecognized lossless compressor in meta_lossless_compress()\n");
        }
        return outSize;
    }

    unsigned long
    meta_lossless_decompress(int losslessCompressor, unsigned char *compressBytes, unsigned long cmpSize, unsigned char **oriData,
                           unsigned long targetOriSize) {
        unsigned long outSize = 0;
        switch (losslessCompressor) {
            case ZSTD_COMPRESSOR:
                *oriData = (unsigned char *) malloc(targetOriSize);
                ZSTD_decompress(*oriData, targetOriSize, compressBytes, cmpSize);
                outSize = targetOriSize;
                break;
            default:
                printf("Error: Unrecognized lossless compressor in meta_lossless_decompress()\n");
        }
        return outSize;
    }
}
#endif