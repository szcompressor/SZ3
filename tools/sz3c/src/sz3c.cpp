//
// Created by Kai Zhao on 10/27/22.
//

#include "sz3c.h"
#include "SZ3/api/sz.hpp"
/** Begin errorbound mode in SZ2 (defines.h) **/
#define ABS 0
#define REL 1
#define VR_REL 1  //alternative name to REL
#define ABS_AND_REL 2
#define ABS_OR_REL 3
#define PSNR 4
#define NORM 5

#define PW_REL 10
#define ABS_AND_PW_REL 11
#define ABS_OR_PW_REL 12
#define REL_AND_PW_REL 13
#define REL_OR_PW_REL 14
/** End errorbound mode in SZ2 (defines.h) **/

/** Begin dataType in SZ2 (defines.h) **/
#define SZ_FLOAT 0
#define SZ_DOUBLE 1
#define SZ_UINT8 2
#define SZ_INT8 3
#define SZ_UINT16 4
#define SZ_INT16 5
#define SZ_UINT32 6
#define SZ_INT32 7
#define SZ_UINT64 8
#define SZ_INT64 9
/** End dataType in SZ2 (defines.h) **/

using namespace SZ;

unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize,
                                int errBoundMode, double absErrBound, double relBoundRatio, double pwrBoundRatio,
                                size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {

    SZ::Config conf;
    if (r2 == 0) {
        conf = SZ::Config(r1);
    } else if (r3 == 0) {
        conf = SZ::Config(r2, r1);
    } else if (r4 == 0) {
        conf = SZ::Config(r3, r2, r1);
    } else if (r5 == 0) {
        conf = SZ::Config(r4, r3, r2, r1);
    } else {
        conf = SZ::Config(r5 * r4, r3, r2, r1);
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
        printf("errBoundMode %d not support\n ", errBoundMode);
        exit(0);
    }

    unsigned char *cmpr_data = NULL;
    if (dataType == SZ_FLOAT) {
        cmpr_data = (unsigned char *) SZ_compress<float>(conf, (float *) data, *outSize);
    } else if (dataType == SZ_DOUBLE) {
        cmpr_data = (unsigned char *) SZ_compress<double>(conf, (double *) data, *outSize);
    } else {
        printf("dataType %d not support\n", dataType);
        exit(0);
    }

    auto *cmpr = (unsigned char *) malloc(*outSize);
    memcpy(cmpr, cmpr_data, *outSize);
    delete[]cmpr_data;

    return cmpr;

}

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength,
                    size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
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

    SZ::Config conf;
    if (dataType == SZ_FLOAT) {
        auto dec_data = (float *) malloc(n * sizeof(float));
        SZ_decompress<float>(conf, (char *) bytes, byteLength, dec_data);
        return dec_data;
    } else if (dataType == SZ_DOUBLE) {
        auto dec_data = (double *) malloc(n * sizeof(double));
        return SZ_decompress<double>(conf, (char *) bytes, byteLength);
        return dec_data;
    } else {
        printf("dataType %d not support\n", dataType);
        exit(0);
    }
}