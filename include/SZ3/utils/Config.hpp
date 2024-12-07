//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_Config_HPP
#define SZ_Config_HPP

#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/inih/INIReader.h"
#include "SZ3/version.hpp"

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

namespace SZ3 {

enum EB { EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL };
constexpr const char *EB_STR[] = {"ABS", "REL", "PSNR", "NORM", "ABS_AND_REL", "ABS_OR_REL"};
constexpr EB EB_OPTIONS[] = {EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL};

enum ALGO { ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP, ALGO_NOPRED, ALGO_LOSSLESS };
constexpr const char *ALGO_STR[] = {"ALGO_LORENZO_REG", "ALGO_INTERP_LORENZO", "ALGO_INTERP",
                                    "ALGO_NOPRED",      "ALGO_LOSSLESS"       };

constexpr const ALGO ALGO_OPTIONS[] = {ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP,
                                       ALGO_NOPRED,      ALGO_LOSSLESS};

enum INTERP_ALGO { INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC };
constexpr const char *INTERP_ALGO_STR[] = {"INTERP_ALGO_LINEAR", "INTERP_ALGO_CUBIC"};
constexpr INTERP_ALGO INTERP_ALGO_OPTIONS[] = {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC};

template <class T>
const char *enum2Str(T e) {
    if (std::is_same<T, ALGO>::value) {
        return ALGO_STR[e];
    } else if (std::is_same<T, INTERP_ALGO>::value) {
        return INTERP_ALGO_STR[e];
    } else if (std::is_same<T, EB>::value) {
        return EB_STR[e];
    } else {
        fprintf(stderr, "invalid enum type for enum2Str()\n ");
        throw std::invalid_argument("invalid enum type for enum2Str()");
    }
}

class Config {
   public:
    template <class... Dims>
    Config(Dims... args) {
        dims = std::vector<size_t>{static_cast<size_t>(std::forward<Dims>(args))...};
        setDims(dims.begin(), dims.end());
    }

    template <class Iter>
    size_t setDims(Iter begin, Iter end) {
        auto dims_ = std::vector<size_t>(begin, end);
        dims.clear();
        for (auto dim : dims_) {
            if (dim > 1) {
                dims.push_back(dim);
            }
        }
        if (dims.empty()) {
            dims = {1};
        }
        N = dims.size();
        num = std::accumulate(dims.begin(), dims.end(), static_cast<size_t>(1), std::multiplies<size_t>());
        pred_dim = N;
        blockSize = (N == 1 ? 128 : (N == 2 ? 16 : 6));
        stride = blockSize;
        return num;
    }

    void loadcfg(const std::string &cfgpath) {
        INIReader cfg(cfgpath);

        if (cfg.ParseError() != 0) {
            std::cerr << "Can't load cfg file " << cfgpath << std::endl;
            throw std::invalid_argument("Can't load cfg file");
        }

        auto cmprAlgoStr = cfg.Get("GlobalSettings", "CmprAlgo", "");
        if (!cmprAlgoStr.empty()) {
            if (cmprAlgoStr == ALGO_STR[ALGO_LORENZO_REG]) {
                cmprAlgo = ALGO_LORENZO_REG;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP_LORENZO]) {
                cmprAlgo = ALGO_INTERP_LORENZO;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_INTERP]) {
                cmprAlgo = ALGO_INTERP;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_NOPRED]) {
                cmprAlgo = ALGO_NOPRED;
            } else if (cmprAlgoStr == ALGO_STR[ALGO_LOSSLESS]) {
                cmprAlgo = ALGO_LOSSLESS;
            } else {
                fprintf(stderr, "Unknown compression algorithm %s\n", cmprAlgoStr.data());
                throw std::invalid_argument("Unknown compression algorithm");
            }
        }
        auto ebModeStr = cfg.Get("GlobalSettings", "ErrorBoundMode", "");
        if (!ebModeStr.empty()) {
            if (ebModeStr == EB_STR[EB_ABS]) {
                errorBoundMode = EB_ABS;
            } else if (ebModeStr == EB_STR[EB_REL]) {
                errorBoundMode = EB_REL;
            } else if (ebModeStr == EB_STR[EB_PSNR]) {
                errorBoundMode = EB_PSNR;
            } else if (ebModeStr == EB_STR[EB_L2NORM]) {
                errorBoundMode = EB_L2NORM;
            } else if (ebModeStr == EB_STR[EB_ABS_AND_REL]) {
                errorBoundMode = EB_ABS_AND_REL;
            } else if (ebModeStr == EB_STR[EB_ABS_OR_REL]) {
                errorBoundMode = EB_ABS_OR_REL;
            } else {
                fprintf(stderr, "Unknown error bound mode %s\n", ebModeStr.data());
                throw std::invalid_argument("Unknown error bound mode");
            }
        }
        absErrorBound = cfg.GetReal("GlobalSettings", "AbsErrorBound", absErrorBound);
        relErrorBound = cfg.GetReal("GlobalSettings", "RelErrorBound", relErrorBound);
        psnrErrorBound = cfg.GetReal("GlobalSettings", "PSNRErrorBound", psnrErrorBound);
        l2normErrorBound = cfg.GetReal("GlobalSettings", "L2NormErrorBound", l2normErrorBound);

        openmp = cfg.GetBoolean("GlobalSettings", "OpenMP", openmp);
        lorenzo = cfg.GetBoolean("AlgoSettings", "Lorenzo", lorenzo);
        lorenzo2 = cfg.GetBoolean("AlgoSettings", "Lorenzo2ndOrder", lorenzo2);
        regression = cfg.GetBoolean("AlgoSettings", "Regression", regression);
        regression2 = cfg.GetBoolean("AlgoSettings", "Regression2ndOrder", regression2);

        auto interpAlgoStr = cfg.Get("AlgoSettings", "InterpolationAlgo", "");
        if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_LINEAR]) {
            interpAlgo = INTERP_ALGO_LINEAR;
        } else if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_CUBIC]) {
            interpAlgo = INTERP_ALGO_CUBIC;
        }
        interpDirection = cfg.GetInteger("AlgoSettings", "InterpolationDirection", interpDirection);
        blockSize = cfg.GetInteger("AlgoSettings", "BlockSize", blockSize);
        quantbinCnt = cfg.GetInteger("AlgoSettings", "QuantizationBinTotal", quantbinCnt);
    }

    size_t save(unsigned char *&c) {
        auto c0 = c;
        write(sz3MagicNumber, c);
        write(sz3DataVer, c);
        write(N, c);
        // write(dims.data(), dims.size(), c);
        auto bitWidth = vector_bit_width(dims);
        write(bitWidth, c);
        vector2bytes(dims, bitWidth, c);

        write(num, c);
        write(cmprAlgo, c);

        write(errorBoundMode, c);
        if (errorBoundMode == EB_ABS) {
            write(absErrorBound, c);
        } else if (errorBoundMode == EB_REL) {
            write(relErrorBound, c);
        } else if (errorBoundMode == EB_PSNR) {
            write(psnrErrorBound, c);
        } else if (errorBoundMode == EB_L2NORM) {
            write(l2normErrorBound, c);
        } else if (errorBoundMode == EB_ABS_OR_REL) {
            write(absErrorBound, c);
            write(relErrorBound, c);
        } else if (errorBoundMode == EB_ABS_AND_REL) {
            write(absErrorBound, c);
            write(relErrorBound, c);
        }

        uint8_t boolvals = (lorenzo & 1) << 7 | (lorenzo2 & 1) << 6 | (regression & 1) << 5 | (regression2 & 1) << 4 |
                           (openmp & 1) << 3;
        write(boolvals, c);

        write(dataType, c);
        write(lossless, c);
        write(encoder, c);
        write(interpAlgo, c);
        write(interpDirection, c);

        write(quantbinCnt, c);
        write(blockSize, c);
        write(stride, c);
        write(pred_dim, c);

        // printf("%lu\n", c - c0);
        return c - c0;
    }

    void load(const unsigned char *&c) {
        // auto c0 = c;
        read(sz3MagicNumber, c);
        if (sz3MagicNumber != SZ3_MAGIC_NUMBER) {
            fprintf(stderr, "magic number mismatch, the input data is not compressed by SZ3\n");
            throw std::invalid_argument("magic number mismatch, the input data is not compressed by SZ3");
        }
        read(sz3DataVer, c);
        if (versionStr(sz3DataVer) != SZ3_DATA_VER) {
            std::stringstream ss;
            printf("program v%s , program-data %s , input data v%s\n", SZ3_VER, SZ3_DATA_VER,
                   versionStr(sz3DataVer).data());
            ss << "Please use SZ3 v" << versionStr(sz3DataVer) << " to decompress the data" << std::endl;
            std::cerr<< ss.str();
            throw std::invalid_argument(ss.str());
        }

        read(N, c);

        uint8_t bitWidth;
        read(bitWidth, c);
        dims = bytes2vector<size_t>(c, bitWidth, N);
        // dims.resize(N);
        // read(dims.data(), N, c);
        read(num, c);
        read(cmprAlgo, c);

        read(errorBoundMode, c);
        if (errorBoundMode == EB_ABS) {
            read(absErrorBound, c);
        } else if (errorBoundMode == EB_REL) {
            read(relErrorBound, c);
        } else if (errorBoundMode == EB_PSNR) {
            read(psnrErrorBound, c);
        } else if (errorBoundMode == EB_L2NORM) {
            read(l2normErrorBound, c);
        } else if (errorBoundMode == EB_ABS_OR_REL) {
            read(absErrorBound, c);
            read(relErrorBound, c);
        } else if (errorBoundMode == EB_ABS_AND_REL) {
            read(absErrorBound, c);
            read(relErrorBound, c);
        }

        uint8_t boolvals;
        read(boolvals, c);
        lorenzo = (boolvals >> 7) & 1;
        lorenzo2 = (boolvals >> 6) & 1;
        regression = (boolvals >> 5) & 1;
        regression2 = (boolvals >> 4) & 1;
        openmp = (boolvals >> 3) & 1;

        read(dataType, c);
        read(lossless, c);
        read(encoder, c);
        read(interpAlgo, c);
        read(interpDirection, c);

        read(quantbinCnt, c);
        read(blockSize, c);
        read(stride, c);
        read(pred_dim, c);

        // print();
        // printf("%d\n", c - c0);
        // if (c - c0 > size_est()) {
        //     throw std::invalid_argument("Config loaded size is larger than estimated size");
        // } else {
        //     // for padding
        //     c = c0 + size_est();
        // }
    }

    void print() {
        printf("===================== Begin SZ3 Configuration =====================\n");
        // printf("sz3MagicNumber = %u\n", sz3MagicNumber);
        printf("sz3DataVer = %s\n", versionStr(sz3DataVer).data());
        printf("N = %d\n", N);
        printf("dims = ");
        for (auto dim : dims) {
            printf("%zu ", dim);
        }
        printf("\nnum = %zu\n", num);
        printf("CmprAlgo = %s\n", enum2Str(static_cast<ALGO>(cmprAlgo)));
        printf("ErrorBoundMode = %s\n", enum2Str(static_cast<EB>(errorBoundMode)));
        printf("AbsErrorBound = %f\n", absErrorBound);
        printf("RelErrorBound = %f\n", relErrorBound);
        printf("PSNRErrorBound = %f\n", psnrErrorBound);
        printf("L2NormErrorBound = %f\n", l2normErrorBound);
        printf("Lorenzo = %d\n", lorenzo);
        printf("Lorenzo2ndOrder = %d\n", lorenzo2);
        printf("Regression = %d\n", regression);
        printf("Regression2ndOrder = %d\n", regression2);
        printf("OpenMP = %d\n", openmp);
        printf("DataType = %d\n", dataType);
        printf("Lossless = %d\n", lossless);
        printf("Encoder = %d\n", encoder);
        printf("InterpolationAlgo = %s\n", enum2Str(static_cast<INTERP_ALGO>(interpAlgo)));
        printf("InterpolationDirection = %d\n", interpDirection);
        printf("QuantizationBinTotal = %d\n", quantbinCnt);
        printf("BlockSize = %d\n", blockSize);
        printf("Stride = %d\n", stride);
        printf("PredDim = %d\n", pred_dim);
        printf("===================== End SZ3 Configuration =====================\n");
    }

    static size_t size_est() {
        return sizeof(Config) + sizeof(size_t) * 5;  // sizeof(size_t) * 5 is for dims vector
    }

    uint32_t sz3MagicNumber = SZ3_MAGIC_NUMBER;
    uint32_t sz3DataVer = versionInt(SZ3_DATA_VER);
    char N = 0;
    std::vector<size_t> dims;
    size_t num = 0;
    uint8_t cmprAlgo = ALGO_INTERP_LORENZO;
    uint8_t errorBoundMode = EB_ABS;
    double absErrorBound = 1e-3;
    double relErrorBound = 0.0;
    double psnrErrorBound = 0.0;
    double l2normErrorBound = 0.0;
    bool lorenzo = true;
    bool lorenzo2 = false;
    bool regression = true;
    bool regression2 = false;
    bool openmp = false;
    uint8_t dataType = SZ_FLOAT;  // dataType is only used in HDF5 filter
    uint8_t lossless = 1;         // 0-> skip lossless(use lossless_bypass); 1-> zstd
    uint8_t encoder = 1;          // 0-> skip encoder; 1->HuffmanEncoder; 2->ArithmeticEncoder
    uint8_t interpAlgo = INTERP_ALGO_CUBIC;
    uint8_t interpDirection = 0;
    int quantbinCnt = 65536;
    int blockSize = 0;
    int stride = 0;        // not used now
    uint8_t pred_dim = 0;  // not used now
};

}  // namespace SZ3

#endif  // SZ_CONFIG_HPP
