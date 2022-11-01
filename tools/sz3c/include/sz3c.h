//
// Created by Kai Zhao on 10/27/22.
//

#ifndef SZ3_SZ3C_H
#define SZ3_SZ3C_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif
unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize,
                                int errBoundMode, double absErrBound, double relBoundRatio, double pwrBoundRatio,
                                size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength,
                    size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

#ifdef __cplusplus
}
#endif

#endif //SZ3_SZ3C_H
