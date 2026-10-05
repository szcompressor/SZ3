//
// Created by Kai Zhao on 10/27/22.
//

#include "sz3c.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

#include "SZ3/api/sz.hpp"

using namespace SZ3;

// Errors return NULL with a message on stderr: a library must not end its caller's process, and an exception must not
// cross into C.
unsigned char *SZ_compress_args(int dataType, void *data, size_t *outSize, int errBoundMode, double absErrBound,
                                double relBoundRatio, double pwrBoundRatio, size_t r5, size_t r4, size_t r3, size_t r2,
                                size_t r1) try {
    SZ3::Config conf;
    if (r2 == 0) {
        conf = SZ3::Config(r1);
    } else if (r3 == 0) {
        conf = SZ3::Config(r2, r1);
    } else if (r4 == 0) {
        conf = SZ3::Config(r3, r2, r1);
    } else if (r5 == 0) {
        conf = SZ3::Config(r4, r3, r2, r1);
    } else {
        conf = SZ3::Config(r5 * r4, r3, r2, r1);
    }
    //    conf.loadcfg(conPath);
    conf.absErrorBound = absErrBound;
    conf.relErrorBound = relBoundRatio;
    //    conf.pwrErrorBound = pwrBoundRatio;
    if (errBoundMode == ABS) {
        conf.errorBoundMode = EB_ABS;
    } else if (errBoundMode == REL) {
        conf.errorBoundMode = EB_REL;
    } else if (errBoundMode == ABS_AND_REL) {
        conf.errorBoundMode = EB_ABS_AND_REL;
    } else if (errBoundMode == ABS_OR_REL) {
        conf.errorBoundMode = EB_ABS_OR_REL;
    } else {
        fprintf(stderr, "SZ3: errBoundMode %d not supported\n", errBoundMode);
        return nullptr;
    }

    unsigned char *cmpr_data = nullptr;
    if (dataType == SZ_FLOAT) {
        cmpr_data = reinterpret_cast<unsigned char *>(SZ_compress<float>(conf, static_cast<float *>(data), *outSize));
#if (!SZ3_DEBUG_TIMINGS)
    } else if (dataType == SZ_DOUBLE) {
        cmpr_data = reinterpret_cast<unsigned char *>(SZ_compress<double>(conf, static_cast<double *>(data), *outSize));
#endif
    } else {
        fprintf(stderr, "SZ3: dataType %d not supported\n", dataType);
        return nullptr;
    }

    // convert c++ memory (by 'new' operator) to c memory (by malloc)
    auto *cmpr = static_cast<unsigned char *>(malloc(*outSize));
    if (cmpr) memcpy(cmpr, cmpr_data, *outSize);
    delete[] cmpr_data;

    return cmpr;
} catch (const std::exception &e) {
    fprintf(stderr, "SZ3: %s\n", e.what());
    return nullptr;
}

// Decompresses into SZ3's own buffer first: the caller's dimensions are only checked against the stream once it is
// read.
template <class T>
static void *decompress(unsigned char *bytes, size_t byteLength, size_t n) {
    SZ3::Config conf;
    std::unique_ptr<T[]> dec(SZ_decompress<T>(conf, reinterpret_cast<char *>(bytes), byteLength));
    if (conf.num != n) {
        fprintf(stderr, "SZ3: the data holds %zu values, not the %zu the dimensions give\n", conf.num, n);
        return nullptr;
    }
    auto *out = static_cast<T *>(malloc(n * sizeof(T)));
    if (out) memcpy(out, dec.get(), n * sizeof(T));
    return out;
}

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2,
                    size_t r1) try {
    size_t n = 0;
    if (r2 == 0) {
        n = r1;
    } else if (r3 == 0) {
        n = r1 * r2;
    } else if (r4 == 0) {
        n = r1 * r2 * r3;
    } else if (r5 == 0) {
        n = r1 * r2 * r3 * r4;
    } else {
        n = r1 * r2 * r3 * r4 * r5;
    }

    if (dataType == SZ_FLOAT) {
        return decompress<float>(bytes, byteLength, n);
#if (!SZ3_DEBUG_TIMINGS)
    } else if (dataType == SZ_DOUBLE) {
        return decompress<double>(bytes, byteLength, n);
#endif
    } else {
        fprintf(stderr, "SZ3: dataType %d not supported\n", dataType);
        return nullptr;
    }
} catch (const std::exception &e) {
    fprintf(stderr, "SZ3: %s\n", e.what());
    return nullptr;
}

void free_buf(void *p) { free(p); }
