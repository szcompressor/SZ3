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

unsigned char *SZ_compress_args(int dataType, void *data, size_t *outSize,
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

    if (dataType == SZ_FLOAT) {
        return (unsigned char *) SZ_compress<float>(conf, (float *) data, *outSize);
    } else if (dataType == SZ_DOUBLE) {
        return (unsigned char *) SZ_compress<double>(conf, (double *) data, *outSize);
    } else {
        printf("dataType %d not support\n", dataType);
        exit(0);
    }

}

//template<class T>
//void compress(char *inPath, char *cmpPath, SZ::Config conf) {
//    T *data = new T[conf.num];
//    SZ::readfile<T>(inPath, conf.num, data);
//
//    size_t outSize;
//    SZ::Timer timer(true);
//    char *bytes = SZ_compress<T>(conf, data, outSize);
//    double compress_time = timer.stop();
//
//    char outputFilePath[1024];
//    if (cmpPath == nullptr) {
//        sprintf(outputFilePath, "%s.sz", inPath);
//    } else {
//        strcpy(outputFilePath, cmpPath);
//    }
//    SZ::writefile(outputFilePath, bytes, outSize);
//
//    printf("compression ratio = %.2f \n", conf.num * 1.0 * sizeof(T) / outSize);
//    printf("compression time = %f\n", compress_time);
//    printf("compressed data file = %s\n", outputFilePath);
//
//    delete[]data;
//    delete[]bytes;
//}
//
//template<class T>
//void decompress(char *inPath, char *cmpPath, char *decPath,
//                SZ::Config conf,
//                int binaryOutput, int printCmpResults) {
//
//    size_t cmpSize;
//    auto cmpData = SZ::readfile<char>(cmpPath, cmpSize);
//
//    SZ::Timer timer(true);
//    T *decData = SZ_decompress<T>(conf, cmpData.get(), cmpSize);
//    double compress_time = timer.stop();
//
//    char outputFilePath[1024];
//    if (decPath == nullptr) {
//        sprintf(outputFilePath, "%s.out", cmpPath);
//    } else {
//        strcpy(outputFilePath, decPath);
//    }
//    if (binaryOutput == 1) {
//        SZ::writefile<T>(outputFilePath, decData, conf.num);
//    } else {
//        SZ::writeTextFile<T>(outputFilePath, decData, conf.num);
//    }
//    if (printCmpResults) {
//        //compute the distortion / compression errors...
//        size_t totalNbEle;
//        auto ori_data = SZ::readfile<T>(inPath, totalNbEle);
//        assert(totalNbEle == conf.num);
//        SZ::verify<T>(ori_data.get(), decData, conf.num);
//    }
//    delete[]decData;
//
//    printf("compression ratio = %f\n", conf.num * sizeof(T) * 1.0 / cmpSize);
//    printf("decompression time = %f seconds.\n", compress_time);
//    printf("decompressed file = %s\n", outputFilePath);
//}
