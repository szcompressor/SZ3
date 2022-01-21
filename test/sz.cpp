#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "SZ3/api/sz.hpp"


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

void usage() {
    printf("Usage: sz <options>\n");
    printf("Options:\n");
    printf("* operation type:\n");
    printf("	-z <compressed file>: the compression operation with an optionally specified output file.\n");
    printf("                          (the compressed file will be named as <input_file>.sz if not specified)\n");
    printf("	-x <decompressed file>: the decompression operation with an optionally specified output file\n");
    printf("                      (the decompressed file will be named as <cmpred_file>.out if not specified)\n");
    printf("	-p: print meta data (configuration info)\n");
    printf("	-h: print the help information\n");
    printf("	-v: print the version number\n");
    printf("* data type:\n");
    printf("	-f: single precision (float type)\n");
    printf("	-d: double precision (double type)\n");
//    printf("* configuration file: \n");
//    printf("	-c <configuration file> : configuration file sz.config\n");
    printf("* error control: (the error control parameters here will overwrite the setting in sz.config)\n");
    printf("	-M <error bound mode> : 10 options as follows. \n");
    printf("		ABS (absolute error bound)\n");
    printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
//    printf("		ABS_AND_REL (using min{ABS, REL})\n");
//    printf("		ABS_OR_REL (using max{ABS, REL})\n");
//    printf("		PSNR (peak signal-to-noise ratio)\n");
//    printf("		NORM (norm2 error : sqrt(sum(xi-xi')^2)\n");
//    printf("		PW_REL (point-wise relative error bound)\n");
    printf("	-A <absolute error bound>: specifying absolute error bound\n");
    printf("	-R <value_range based relative error bound>: specifying relative error bound\n");
//    printf("	-P <point-wise relative error bound>: specifying point-wise relative error bound\n");
//    printf("	-S <PSNR>: specifying PSNR\n");
//    printf("	-N <normErr>: specifying normErr\n");
    printf("* input data file:\n");
    printf("	-i <original data file> : original data file\n");
    printf("	-s <compressed data file> : compressed data file in decompression\n");
    printf("* output type of decompressed file: \n");
    printf("	-b (by default) : decompressed file stored in binary format\n");
//    printf("	-t : decompreadded file stored in text format\n");
//    printf("	-T : pre-processing with Tucker Tensor Decomposition\n");
    printf("* dimensions: \n");
    printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
    printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
    printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
    printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
    printf("* print compression results: \n");
    printf("	-a : print compression results such as distortions\n");
    printf("* examples: \n");
    printf("	sz -z -f -c sz.config -i testdata/x86/testfloat_8_8_128.dat -3 8 8 128\n");
    printf("	sz -z -f -c sz.config -M ABS -A 1E-3 -i testdata/x86/testfloat_8_8_128.dat -3 8 8 128\n");
    printf("	sz -x -f -s testdata/x86/testfloat_8_8_128.dat.sz -3 8 8 128\n");
    printf("	sz -x -f -s testdata/x86/testfloat_8_8_128.dat.sz -i testdata/x86/testfloat_8_8_128.dat -3 8 8 128 -a\n");
    printf("	sz -z -d -c sz.config -i testdata/x86/testdouble_8_8_128.dat -3 8 8 128\n");
    printf("	sz -x -d -s testdata/x86/testdouble_8_8_128.dat.sz -3 8 8 128\n");
    printf("	sz -p -s testdata/x86/testdouble_8_8_128.dat.sz\n");
    exit(0);
}

template<class T>
void compress(char *inPath, char *cmpPath, SZ::Config conf) {
    T *data = new T[conf.num];
    SZ::readfile<T>(inPath, conf.num, data);

    size_t outSize;
    SZ::Timer timer(true);
    char *bytes = SZ_compress<T>(conf, data, outSize);
    double compress_time = timer.stop();

    char outputFilePath[1024];
    if (cmpPath == NULL) {
        sprintf(outputFilePath, "%s.sz", inPath);
    } else {
        strcpy(outputFilePath, cmpPath);
    }
    SZ::writefile(outputFilePath, bytes, outSize);

    printf("compression ratio = %.2f \n", conf.num * 1.0 * sizeof(T) / outSize);
    printf("compression time = %f\n", compress_time);
    printf("compressed data file = %s\n", outputFilePath);

    delete[]data;
    delete[]bytes;
}

template<class T>
void decompress(char *inPath, char *cmpPath, char *decPath,
                SZ::Config conf,
                int binaryOutput, int printCmpResults) {

    size_t cmpSize;
    auto cmpData = SZ::readfile<char>(cmpPath, cmpSize);

    SZ::Timer timer(true);
    T *decData = SZ_decompress<T>(conf, cmpData.get(), cmpSize);
    double compress_time = timer.stop();

    char outputFilePath[1024];
    if (decPath == NULL) {
        sprintf(outputFilePath, "%s.out", cmpPath);
    } else {
        strcpy(outputFilePath, decPath);
    }
    if (binaryOutput == 1) {
        SZ::writefile<T>(outputFilePath, decData, conf.num);
    } else { //txt output
//        writeFloatData(decData, nbEle, outputFilePath, &status);
    }
    if (printCmpResults) {
        //compute the distortion / compression errors...
        size_t totalNbEle;
        auto ori_data = SZ::readfile<T>(inPath, totalNbEle);
        assert(totalNbEle == conf.num);
        SZ::verify<T>(ori_data.get(), decData, conf.num);
    }
    delete[]decData;

    printf("compression ratio = %f\n", conf.num * sizeof(T) * 1.0 / cmpSize);
    printf("decompression time = %f seconds.\n", compress_time);
    printf("decompressed file = %s\n", outputFilePath);
}

int main(int argc, char *argv[]) {
    int binaryOutput = 1;
    int printCmpResults = 0;
    int compressionMode = 0; // &1 : compression ; &2: decompression
    int dataType = SZ_FLOAT;
    char *inPath = NULL;
    char *cmpPath = NULL;
    char *conPath = NULL;
    char *decPath = NULL;

    char *errBoundMode = NULL;
    char *absErrorBound = NULL;
    char *relErrorBound = NULL;
    char *pwrErrorBound = NULL;
    char *psnr_ = NULL;
    char *normError = NULL;

    size_t r5 = 0;
    size_t r4 = 0;
    size_t r3 = 0;
    size_t r2 = 0;
    size_t r1 = 0;

    size_t i = 0;
    int status;
    size_t nbEle;
    if (argc == 1)
        usage();

    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-' || argv[i][2])
            usage();
        switch (argv[i][1]) {
            case 'h':
                usage();
                exit(0);
            case 'v':
                printf("version: %d.%d.%d\n", SZ_VERSION_MAJOR, SZ_VERSION_MINOR, SZ_VERSION_RELEASE);
                exit(0);
            case 'b':
                binaryOutput = 1;
                break;
//            case 't':
//                binaryOutput = 0;
//                break;
            case 'a':
                printCmpResults = 1;
                break;
            case 'z':
                compressionMode |= 1;
                if (i + 1 < argc) {
                    cmpPath = argv[i + 1];
                    if (cmpPath[0] != '-')
                        i++;
                    else
                        cmpPath = NULL;
                }
                break;
            case 'x':
                compressionMode |= 2;
                if (i + 1 < argc) {
                    decPath = argv[i + 1];
                    if (decPath[0] != '-')
                        i++;
                    else
                        decPath = NULL;
                }
                break;
            case 'f':
                dataType = SZ_FLOAT;
                break;
            case 'd':
                dataType = SZ_DOUBLE;
                break;
            case 'i':
                if (++i == argc)
                    usage();
                inPath = argv[i];
                break;
            case 's':
                if (++i == argc)
                    usage();
                cmpPath = argv[i];
                break;
            case 'c':
                if (++i == argc)
                    usage();
                conPath = argv[i];
                break;
            case '1':
                if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1)
                    usage();
                break;
            case '2':
                if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r2) != 1)
                    usage();
                break;
            case '3':
                if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r3) != 1)
                    usage();
                break;
            case '4':
                if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r3) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r4) != 1)
                    usage();
                break;
            case 'M':
                if (++i == argc)
                    usage();
                errBoundMode = argv[i];
                break;
            case 'A':
                if (++i == argc)
                    usage();
                absErrorBound = argv[i];
                break;
            case 'R':
                if (++i == argc)
                    usage();
                relErrorBound = argv[i];
                break;
//            case 'P':
//                if (++i == argc)
//                    usage();
//                pwrErrorBound = argv[i];
//                break;
//            case 'N':
//                if (++i == argc)
//                    usage();
//                normError = argv[i];
//                break;
//            case 'S':
//                if (++i == argc)
//                    usage();
//                psnr_ = argv[i];
//                break;
            default:
                usage();
                break;
        }
    }

    if ((inPath == NULL) & (cmpPath == NULL)) {
        printf("Error: you need to specify either a raw binary data file or a compressed data file as input\n");
        usage();
        exit(0);
    }


    SZ::Config conf;
    if (r2 == 0) {
        nbEle = r1;
        conf = SZ::Config(r1);
    } else if (r3 == 0) {
        nbEle = r1 * r2;
        conf = SZ::Config(r2, r1);
    } else if (r4 == 0) {
        nbEle = r1 * r2 * r3;
        conf = SZ::Config(r3, r2, r1);
    } else if (r5 == 0) {
        nbEle = r1 * r2 * r3 * r4;
        conf = SZ::Config(r4, r3, r2, r1);
    } else {
        nbEle = r1 * r2 * r3 * r4 * r5;
        conf = SZ::Config(r5, r4, r3, r2, r1);
    }

    if (errBoundMode != NULL) {
        if (strcmp(errBoundMode, "ABS") == 0)
            conf.errorBoundMode = ABS;
        else if (strcmp(errBoundMode, "REL") == 0 || strcmp(errBoundMode, "VR_REL") == 0)
            conf.errorBoundMode = REL;
//        else if (strcmp(errBoundMode, "ABS_AND_REL") == 0)
//            conf.errorBoundMode = ABS_AND_REL;
//        else if (strcmp(errBoundMode, "ABS_OR_REL") == 0)
//            conf.errorBoundMode = ABS_OR_REL;
//        else if (strcmp(errBoundMode, "PSNR") == 0)
//            conf.errorBoundMode = PSNR;
//        else if (strcmp(errBoundMode, "PW_REL") == 0)
//            conf.errorBoundMode = PW_REL;
//        else if (strcmp(errBoundMode, "NORM") == 0)
//            conf.errorBoundMode = NORM;
        else {
            printf("Error: wrong error bound mode setting by using the option '-M'\n");
            usage();
            exit(0);
        }
    }

    if (compressionMode & 1) {
        if (absErrorBound != NULL) {
            conf.absErrorBound = atof(absErrorBound);
        }
        if (relErrorBound != NULL) {
            conf.relErrorBound = atof(relErrorBound);
        }

        if (dataType == SZ_FLOAT) {
            compress<float>(inPath, cmpPath, conf);
        } else if (dataType == SZ_DOUBLE) {
            compress<double>(inPath, cmpPath, conf);
        } else {
            printf("Error: data type not supported \n");
            usage();
            exit(0);
        }
    }
    if (compressionMode & 2) { //decompression
        if (printCmpResults) {
            if (inPath == NULL) {
                printf("Error: Since you add -a option (analysis), please specify the original data path by -i <path>.\n");
                exit(0);
            }
        }

        if (dataType == SZ_FLOAT) {
            decompress<float>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);
        } else if (dataType == SZ_DOUBLE) {
            decompress<double>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);
        } else {
            printf("Error: data type not supported \n");
            usage();
            exit(0);
        }
    }
    return 0;
}
