#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

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

inline void usage() {
    printf("Note: SZ3 command line arguments are backward compatible with SZ2, \n");
    printf("      use -h2 to show the supported SZ2 command line arguments. \n");
    printf("Usage: sz3 <options>\n");
    printf("Options:\n");
    printf("* general options:\n");
    printf("	-h: print the help information\n");
    printf("	-h2: print the help information for SZ2 style command line\n");
    printf("	-v: print the version number\n");
    printf("	-a : print compression results such as distortions\n");
    printf("* input and output:\n");
    printf("	-i <path> : original input file in binary format\n");
    printf("	-o <path> : decompressed file in binary format\n");
    printf("	-z <path> : compressed file\n");
    printf("	-t : store decompressed file in text format\n");
    //    printf("	-p: print meta data (configuration info)\n");
    printf("* data type:\n");
    printf("	-f: single precision (float type)\n");
    printf("	-d: double precision (double type)\n");
    printf("	-I <width>: integer type (width = 32 or 64)\n");
    printf("* configuration file: \n");
    printf("	-c <configuration file> : configuration file sz.config\n");
    printf("* error control: (the error control parameters here will overwrite the setting in sz.config)\n");
    printf("	-M <error control mode> <error bound (optional)> \n");
    printf("	error control mode as follows: \n");
    printf("		ABS (absolute error bound)\n");
    printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
    printf("		PSNR (peak signal-to-noise ratio)\n");
    printf("		NORM (norm2 error : sqrt(sum(xi-xi')^2)\n");
    printf("		ABS_AND_REL (using min{ABS, REL})\n");
    printf("		ABS_OR_REL (using max{ABS, REL})\n");
    printf(
        "	error bound can be set directly after the error control mode, or separately with the following "
        "options:\n");
    printf("		-A <absolute error bound>: specifying absolute error bound\n");
    printf("		-R <value_range based relative error bound>: specifying relative error bound\n");
    //    printf("		-P <point-wise relative error bound>: specifying point-wise relative error bound\n");
    printf("		-S <PSNR>: specifying PSNR\n");
    printf("		-N <normErr>: specifying normErr\n");
    printf("* dimensions: \n");
    printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
    printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
    printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
    printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
    printf("* examples: \n");
    printf("	sz -f -i test.dat    -z test.dat.sz     -3 8 8 128 -M ABS 1e-3 \n");
    printf("	sz -f -z test.dat.sz -o test.dat.sz.out -3 8 8 128 -M REL 1e-3 -a \n");
    printf("	sz -f -i test.dat    -o test.dat.sz.out -3 8 8 128 -M ABS_AND_REL -A 1 -R 1e-3 -a \n");
    printf("	sz -f -i test.dat    -o test.dat.sz.out -3 8 8 128 -c sz.config \n");
    printf("	sz -f -i test.dat    -o test.dat.sz.out -3 8 8 128 -c sz.config -M ABS 1e-3 -a\n");
    exit(0);
}

inline void usage_sz2() {
    printf("Note: below are the supported command line arguments in SZ2 style\n");
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
    printf("* configuration file: \n");
    printf("	-c <configuration file> : configuration file sz.config\n");
    printf("* error control: (the error control parameters here will overwrite the setting in sz.config)\n");
    printf("	-M <error bound mode> : 10 options as follows. \n");
    printf("		ABS (absolute error bound)\n");
    printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
    printf("		ABS_AND_REL (using min{ABS, REL})\n");
    printf("		ABS_OR_REL (using max{ABS, REL})\n");
    printf("		PSNR (peak signal-to-noise ratio)\n");
    printf("		NORM (norm2 error : sqrt(sum(xi-xi')^2)\n");
    //    printf("		PW_REL (point-wise relative error bound)\n");
    printf("	-A <absolute error bound>: specifying absolute error bound\n");
    printf("	-R <value_range based relative error bound>: specifying relative error bound\n");
    //    printf("	-P <point-wise relative error bound>: specifying point-wise relative error bound\n");
    printf("	-S <PSNR>: specifying PSNR\n");
    printf("	-N <normErr>: specifying normErr\n");
    printf("* input data file:\n");
    printf("	-i <original data file> : original data file\n");
    printf("	-s <compressed data file> : compressed data file in decompression\n");
    printf("* output type of decompressed file: \n");
    printf("	-b (by default) : decompressed file stored in binary format\n");
    printf("	-t : decompreadded file stored in text format\n");
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
    printf(
        "	sz -x -f -s testdata/x86/testfloat_8_8_128.dat.sz -i testdata/x86/testfloat_8_8_128.dat -3 8 8 128 "
        "-a\n");
    printf("	sz -z -d -c sz.config -i testdata/x86/testdouble_8_8_128.dat -3 8 8 128\n");
    printf("	sz -x -d -s testdata/x86/testdouble_8_8_128.dat.sz -3 8 8 128\n");
    printf("	sz -p -s testdata/x86/testdouble_8_8_128.dat.sz\n");
    exit(0);
}

// A malformed argument is a caller error, not a help request: report it and let main() exit non-zero.
[[noreturn]] static void argError(const std::string &msg) { throw std::invalid_argument(msg); }

// A dimension is a positive element count. Reject signs, non-digits, trailing junk, zero, and overflow
// here, so the SZ3 sizes derived from it downstream never come from a corrupt or absurd value.
static size_t parse_dim(const char *s) {
    if (s[0] < '0' || s[0] > '9') argError(std::string("invalid dimension '") + s + "'");
    errno = 0;
    char *end = nullptr;
    unsigned long long v = strtoull(s, &end, 10);
    bool overflows_size_t = false;
    if constexpr (sizeof(size_t) < sizeof(unsigned long long)) {
        overflows_size_t = v > static_cast<unsigned long long>(SIZE_MAX);
    }
    if (*end != '\0' || errno == ERANGE || v == 0 || overflows_size_t)
        argError(std::string("invalid dimension '") + s + "'");
    return static_cast<size_t>(v);
}

template <class T>
void compress(char *inPath, char *cmpPath, SZ3::Config &conf) {
    T *data = new T[conf.num];
    SZ3::readfile<T>(inPath, conf.num, data);
    // SZ_compress refuses a capacity below its own bound, which a small input falls under.
    size_t bytesCap = SZ3::SZ_compress_size_bound<T>(conf);
    auto bytes = new char[bytesCap];

    SZ3::Timer timer(true);
    size_t outSize = SZ_compress<T>(conf, data, bytes, bytesCap);
    double compress_time = timer.stop();

    char outputFilePath[1024];
    if (cmpPath == nullptr) {
        snprintf(outputFilePath, 1024, "%s.sz", inPath);
    } else {
        strcpy(outputFilePath, cmpPath);
    }
    SZ3::writefile(outputFilePath, bytes, outSize);

    printf("compression ratio = %.2f \n", conf.num * 1.0 * sizeof(T) / outSize);
    printf("compression time = %f\n", compress_time);
    printf("compressed data file = %s\n", outputFilePath);

    delete[] data;
    delete[] bytes;
}

template <class T>
void decompress(char *inPath, char *cmpPath, char *decPath, SZ3::Config &conf, int binaryOutput, int printCmpResults) {
    size_t cmpSize;
    auto cmpData = SZ3::readfile<char>(cmpPath, cmpSize);

    SZ3::Timer timer(true);
    auto decData = SZ_decompress<T>(conf, cmpData.get(), cmpSize);
    double compress_time = timer.stop();

    char outputFilePath[1024];
    if (decPath == nullptr) {
        snprintf(outputFilePath, 1024, "%s.out", cmpPath);
    } else {
        strcpy(outputFilePath, decPath);
    }
    if (binaryOutput == 1) {
        SZ3::writefile<T>(outputFilePath, decData, conf.num);
    } else {
        SZ3::writeTextFile<T>(outputFilePath, decData, conf.num);
    }
    if (printCmpResults) {
        // compute the distortion / compression errors...
        size_t totalNbEle;
        auto ori_data = SZ3::readfile<T>(inPath, totalNbEle);
        assert(totalNbEle == conf.num);
        SZ3::verify<T>(ori_data.get(), decData, conf.num);
    }
    delete[] decData;

    printf("compression ratio = %f\n", conf.num * sizeof(T) * 1.0 / cmpSize);
    printf("decompression time = %f seconds.\n", compress_time);
    printf("decompressed file = %s\n", outputFilePath);
}

static int run(int argc, char *argv[]) {
    bool binaryOutput = true;
    int printCmpResults = 0;
    int printMeta = 0;
    bool compression = false;
    bool decompression = false;
    int dataType = SZ_FLOAT;
    char *inPath = nullptr;
    char *cmpPath = nullptr;
    char *conPath = nullptr;
    char *decPath = nullptr;
    bool delCmpPath = false;

    char *errBoundMode = nullptr;
    char *errBound = nullptr;
    char *absErrorBound = nullptr;
    char *relErrorBound = nullptr;
    // char *pwrErrorBound = nullptr;
    char *psnrErrorBound = nullptr;
    char *normErrorBound = nullptr;

    bool sz2mode = false;

    size_t r4 = 0;
    size_t r3 = 0;
    size_t r2 = 0;
    size_t r1 = 0;

    int i = 0;
    // int status;
    if (argc == 1) usage();
    int width = -1;

    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-' || argv[i][2]) {
            if (argv[i][1] == 'h' && argv[i][2] == '2') {
                usage_sz2();
            } else if (strcmp(argv[i], "--help") == 0) {
                usage();
            } else {
                argError(std::string("unrecognized argument: ") + argv[i]);
            }
        }
        switch (argv[i][1]) {
            case 'h':
                usage();
                exit(0);
            case 'v':
                printf("SZ3 Version: %s\n", SZ3_VER);
                printf("SZ3 Data Format Version: %s\n", SZ3_DATA_VER);
                printf("\nThird-party libraries copyright notices:\n");
                printf("----------------------------------------\n");
                printf("ska_hash:\n");
                printf("  Copyright (c) 2017 Malte Skarupke\n");
                printf("  Licensed under the Boost Software License - Version 1.0\n");
                exit(0);
            case 'b':
                binaryOutput = true;
                break;
            case 't':
                binaryOutput = false;
                break;
            case 'a':
                printCmpResults = 1;
                break;
            case 'p':
                printMeta = 1;
                break;
            case 'z':
                compression = true;
                if (i + 1 < argc) {
                    cmpPath = argv[i + 1];
                    if (cmpPath[0] != '-')
                        i++;
                    else
                        cmpPath = nullptr;
                }
                break;
            case 'x':
                sz2mode = true;
                decompression = true;
                if (i + 1 < argc) {
                    decPath = argv[i + 1];
                    if (decPath[0] != '-')
                        i++;
                    else
                        decPath = nullptr;
                }
                break;
            case 'f':
                dataType = SZ_FLOAT;
                break;
            case 'd':
                dataType = SZ_DOUBLE;
                break;
            case 'I':
                if (++i == argc || sscanf(argv[i], "%d", &width) != 1) {
                    argError("-I requires an integer width (32 or 64)");
                }
                if (width == 32) {
                    dataType = SZ_INT32;
                } else if (width == 64) {
                    dataType = SZ_INT64;
                } else {
                    argError("-I width must be 32 or 64");
                }
                break;
            case 'i':
                if (++i == argc) argError("-i requires a file path");
                inPath = argv[i];
                break;
            case 'o':
                if (++i == argc) argError("-o requires a file path");
                decPath = argv[i];
                break;
            case 's':
                sz2mode = true;
                if (++i == argc) argError("-s requires a file path");
                cmpPath = argv[i];
                break;
            case 'c':
                if (++i == argc) argError("-c requires a config file path");
                conPath = argv[i];
                break;
            case '1':
                if (++i == argc) argError("-1 requires 1 dimension");
                r1 = parse_dim(argv[i]);
                break;
            case '2':
                if (++i == argc) argError("-2 requires 2 dimensions");
                r1 = parse_dim(argv[i]);
                if (++i == argc) argError("-2 requires 2 dimensions");
                r2 = parse_dim(argv[i]);
                break;
            case '3':
                if (++i == argc) argError("-3 requires 3 dimensions");
                r1 = parse_dim(argv[i]);
                if (++i == argc) argError("-3 requires 3 dimensions");
                r2 = parse_dim(argv[i]);
                if (++i == argc) argError("-3 requires 3 dimensions");
                r3 = parse_dim(argv[i]);
                break;
            case '4':
                if (++i == argc) argError("-4 requires 4 dimensions");
                r1 = parse_dim(argv[i]);
                if (++i == argc) argError("-4 requires 4 dimensions");
                r2 = parse_dim(argv[i]);
                if (++i == argc) argError("-4 requires 4 dimensions");
                r3 = parse_dim(argv[i]);
                if (++i == argc) argError("-4 requires 4 dimensions");
                r4 = parse_dim(argv[i]);
                break;
            case 'M': {
                if (++i == argc) argError("-M requires an error bound mode");
                errBoundMode = argv[i];
                // match_enum silently keeps the default when the name is unknown, which would accept a typo as
                // ABS. Reject a name that is neither in the table nor the VR_REL alias handled during setup.
                bool knownMode = strcmp(errBoundMode, "VR_REL") == 0;
                for (const auto &kv : SZ3::EB_MAP) {
                    if (SZ3::to_lower(kv.first) == SZ3::to_lower(errBoundMode)) knownMode = true;
                }
                if (!knownMode) argError(std::string("unknown error bound mode: ") + errBoundMode);
                if (i + 1 < argc && argv[i + 1][0] != '-') {
                    errBound = argv[++i];
                }
                break;
            }
            case 'A':
                if (++i == argc) argError("-A requires an absolute error bound");
                absErrorBound = argv[i];
                break;
            case 'R':
                if (++i == argc) argError("-R requires a relative error bound");
                relErrorBound = argv[i];
                break;
                //            case 'P':
                //                if (++i == argc)
                //                    usage();
                //                pwrErrorBound = argv[i];
                //                break;
            case 'N':
                if (++i == argc) argError("-N requires a norm error bound");
                normErrorBound = argv[i];
                break;
            case 'S':
                if (++i == argc) argError("-S requires a PSNR value");
                psnrErrorBound = argv[i];
                break;
            default:
                argError(std::string("unknown option: ") + argv[i]);
                break;
        }
    }

    if ((inPath == nullptr) && (cmpPath == nullptr)) {
        argError("you need to specify either a raw binary data file (-i) or a compressed data file (-s/-z) as input");
    }

    if (!sz2mode && inPath != nullptr && cmpPath != nullptr) {
        compression = true;
    }
    if (cmpPath != nullptr && decPath != nullptr) {
        decompression = true;
    }
    char cmpPathTmp[1024];
    if (inPath != nullptr && cmpPath == nullptr && decPath != nullptr) {
        compression = true;
        decompression = true;
        snprintf(cmpPathTmp, 1024, "%s.sz.tmp", inPath);
        cmpPath = cmpPathTmp;
        delCmpPath = true;
    }
    if (inPath == nullptr || (errBoundMode == nullptr && conPath == nullptr)) {
        compression = false;
    }
    if (!compression && !decompression) {
        argError("nothing to do: specify compression (-i with -z/-o and an error bound) or decompression (-s/-x)");
    }

    SZ3::Config conf;
    if (r2 == 0) {
        conf = SZ3::Config(r1);
    } else if (r3 == 0) {
        conf = SZ3::Config(r2, r1);
    } else if (r4 == 0) {
        conf = SZ3::Config(r3, r2, r1);
    } else {
        conf = SZ3::Config(r4, r3, r2, r1);
    }
    if (compression && conPath != nullptr) {
        conf.loadcfg(conPath);
    }

    if (errBoundMode != nullptr) {
        {
            // backward compatible with SZ2
            if (relErrorBound != nullptr) {
                conf.relErrorBound = atof(relErrorBound);
            }
            if (absErrorBound != nullptr) {
                conf.absErrorBound = atof(absErrorBound);
            }
            if (psnrErrorBound != nullptr) {
                conf.psnrErrorBound = atof(psnrErrorBound);
            }
            if (normErrorBound != nullptr) {
                conf.l2normErrorBound = atof(normErrorBound);
            }
        }

        SZ3::match_enum(errBoundMode, SZ3::EB_MAP, conf.errorBoundMode);
        if (strcmp(errBoundMode, "VR_REL") == 0) {
            conf.errorBoundMode = SZ3::EB_REL;
        }
        if (conf.errorBoundMode == SZ3::EB_ABS) {
            if (errBound != nullptr) {
                conf.absErrorBound = atof(errBound);
            }
        } else if (conf.errorBoundMode == SZ3::EB_REL) {
            if (errBound != nullptr) {
                conf.relErrorBound = atof(errBound);
            }
        } else if (conf.errorBoundMode == SZ3::EB_PSNR) {
            if (errBound != nullptr) {
                conf.psnrErrorBound = atof(errBound);
            }
        } else if (conf.errorBoundMode == SZ3::EB_L2NORM) {
            if (errBound != nullptr) {
                conf.l2normErrorBound = atof(errBound);
            }
        } else if (conf.errorBoundMode == SZ3::EB_ABS_AND_REL) {
        } else if (conf.errorBoundMode == SZ3::EB_ABS_OR_REL) {
        } else {
            argError("wrong error bound mode setting by using the option '-M'");
        }
    }

    if (compression) {
        if (dataType == SZ_FLOAT) {
            compress<float>(inPath, cmpPath, conf);
#if (!SZ3_DEBUG_TIMINGS)
        } else if (dataType == SZ_DOUBLE) {
            compress<double>(inPath, cmpPath, conf);
        } else if (dataType == SZ_INT32) {
            compress<int32_t>(inPath, cmpPath, conf);
        } else if (dataType == SZ_INT64) {
            compress<int64_t>(inPath, cmpPath, conf);
#endif
        } else {
            argError("data type not supported");
        }
    }
    if (decompression) {
        if (printCmpResults && inPath == nullptr) {
            argError("the -a option (analysis) needs the original data path, specify it with -i <path>");
        }

        if (dataType == SZ_FLOAT) {
            decompress<float>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);
#if (!SZ3_DEBUG_TIMINGS)
        } else if (dataType == SZ_DOUBLE) {
            decompress<double>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);
        } else if (dataType == SZ_INT32) {
            decompress<int32_t>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);
        } else if (dataType == SZ_INT64) {
            decompress<int64_t>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);
#endif
        } else {
            argError("data type not supported");
        }
    }
    if (printMeta) {
        conf.print();
    }
    if (delCmpPath) {
        remove(cmpPath);
    }
    return 0;
}

// A configuration the algorithms refuse arrives here as an exception; a caller needs a message
// and an exit code, not a signal.
int main(int argc, char *argv[]) {
    try {
        return run(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << "sz3: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "sz3: unknown error" << std::endl;
    }
    return 1;
}
