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
//#include "SZ3/utils/Config_basic.hpp"
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



    enum EB {
        EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL
    };
    constexpr const char *EB_STR[] = {"ABS", "REL", "PSNR", "NORM", "ABS_AND_REL", "ABS_OR_REL"};
    constexpr EB EB_OPTIONS[] = {EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL};

    enum ALGO { ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP, ALGO_NOPRED, ALGO_LOSSLESS };
    constexpr const char *ALGO_STR[] = {"ALGO_LORENZO_REG", "ALGO_INTERP_LORENZO", "ALGO_INTERP",
                                        "ALGO_NOPRED",      "ALGO_LOSSLESS"       };

    constexpr const ALGO ALGO_OPTIONS[] = {ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP,
                                           ALGO_NOPRED,      ALGO_LOSSLESS};


    enum INTERP_ALGO {
        INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC, INTERP_ALGO_QUAD
    };
    constexpr const char *INTERP_ALGO_STR[] = {"INTERP_ALGO_LINEAR", "INTERP_ALGO_CUBIC","INTERP_ALGO_QUAD"};
    constexpr const INTERP_ALGO INTERP_ALGO_OPTIONS[] = { INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC, INTERP_ALGO_QUAD };

    enum TUNING_TARGET {
        TUNING_TARGET_RD, TUNING_TARGET_CR, TUNING_TARGET_SSIM, TUNING_TARGET_AC
    };
    constexpr const char *TUNING_TARGET_STR[] = {"TUNING_TARGET_RD", "TUNING_TARGET_CR", "TUNING_TARGET_SSIM", "TUNING_TARGET_AC"};
    constexpr const TUNING_TARGET TUNING_TARGET_OPTIONS[] = {
        TUNING_TARGET_RD, TUNING_TARGET_CR, TUNING_TARGET_SSIM, TUNING_TARGET_AC
    };


    struct Interp_Meta{
        uint8_t interpAlgo = INTERP_ALGO_CUBIC;
        uint8_t interpParadigm = 0;//1D, MD,HD
        uint8_t cubicSplineType = 0;//noknot,nat
        uint8_t interpDirection = 0;//0,N!-1
        uint8_t adjInterp=0;//0,1
        std::array<float,3> dimCoeffs={1.0/3.0,1.0/3.0,1.0/3.0};
 
    };

    //void print_meta(Interp_Meta meta){
    //    std::cout<<(int)meta.interpAlgo<<" "<<(int)meta.interpParadigm<<" "<<(int)meta.cubicSplineType<<" "<<(int)meta.interpDirection<<" "<<(int)meta.adjInterp<<std::endl;

    //}

    template<class T>
    const char *enum2Str(T e) {
        if (std::is_same<T, ALGO>::value) {
            return ALGO_STR[e];
        } else if (std::is_same<T, INTERP_ALGO>::value) {
            return INTERP_ALGO_STR[e];
        } else if (std::is_same<T, EB>::value) {
            return EB_STR[e];

        }  else if (std::is_same<T, TUNING_TARGET>::value) {
            return TUNING_TARGET_STR[e];
            
        }else {
            printf("invalid enum type for enum2Str()\n ");
            exit(0);
        }
    }

    class Config{
    public:
        template<class ... Dims>
        Config(Dims ... args) {
            dims = std::vector<size_t>{static_cast<size_t>(std::forward<Dims>(args))...};
            setDims(dims.begin(), dims.end());
        }

        template<class Iter>
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
                std::cout << "Can't load cfg file  <<" << cfgpath << std::endl;
                exit(0);
            } else {
                //std::cout << "Load cfg from " << cfgpath << std::endl;
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
            auto tuningTargetStr = cfg.Get("GlobalSettings", "tuningTarget", "");
            if (tuningTargetStr == TUNING_TARGET_STR[TUNING_TARGET_RD]) {
                tuningTarget = TUNING_TARGET_RD;
            } else if (tuningTargetStr == TUNING_TARGET_STR[TUNING_TARGET_CR]) {
                tuningTarget = TUNING_TARGET_CR;
            }
            else if (tuningTargetStr == TUNING_TARGET_STR[TUNING_TARGET_SSIM]) {
                tuningTarget = TUNING_TARGET_SSIM;
            }
            else if (tuningTargetStr == TUNING_TARGET_STR[TUNING_TARGET_AC]) {
                tuningTarget = TUNING_TARGET_AC;
            }
                


            absErrorBound = cfg.GetReal("GlobalSettings", "AbsErrorBound", absErrorBound);
            relErrorBound = cfg.GetReal("GlobalSettings", "RelErrorBound", relErrorBound);
            psnrErrorBound = cfg.GetReal("GlobalSettings", "PSNRErrorBound", psnrErrorBound);
            l2normErrorBound = cfg.GetReal("GlobalSettings", "L2NormErrorBound", l2normErrorBound);
            //prewave_absErrorBound= cfg.GetReal("GlobalSettings", "prewave_absErrorBound", prewave_absErrorBound);
            alpha = cfg.GetReal("AlgoSettings", "alpha", alpha);
            beta = cfg.GetReal("AlgoSettings", "beta", beta);
            autoTuningRate = cfg.GetReal("AlgoSettings", "autoTuningRate", autoTuningRate);
            predictorTuningRate = cfg.GetReal("AlgoSettings", "predictorTuningRate", predictorTuningRate);
           
            lorenzoBrFix = cfg.GetReal("AlgoSettings", "lorenzoBrFix", lorenzoBrFix);
            blockwiseSampleRate=cfg.GetReal("AlgoSettings", "blockwiseSampleRate", blockwiseSampleRate);
            //anchorThreshold= cfg.GetReal("AlgoSettings", "anchorThreshold", anchorThreshold);
            //sperr_eb_coeff = cfg.GetReal("AlgoSettings", "sperr_eb_coeff", sperr_eb_coeff);

            openmp = cfg.GetBoolean("GlobalSettings", "OpenMP", openmp);
            lorenzo = cfg.GetBoolean("AlgoSettings", "Lorenzo", lorenzo);
            lorenzo2 = cfg.GetBoolean("AlgoSettings", "Lorenzo2ndOrder", lorenzo2);
            regression = cfg.GetBoolean("AlgoSettings", "Regression", regression);
            regression2 = cfg.GetBoolean("AlgoSettings", "Regression2ndOrder", regression2);
            
            naturalSpline = cfg.GetBoolean("AlgoSettings", "naturalSpline", naturalSpline);
            
           
            verbose=cfg.GetBoolean("AlgoSettings", "verbose", verbose);
            var_first = cfg.GetBoolean("AlgoSettings", "var_first", var_first);
            blockwiseTuning = cfg.GetBoolean("AlgoSettings", "blockwiseTuning", blockwiseTuning);
            quadInterp = cfg.GetBoolean("AlgoSettings", "quadInterp", quadInterp);
            freezeDimTest = cfg.GetBoolean("AlgoSettings", "freezeDimTest", freezeDimTest);
            dynamicDimCoeff = cfg.GetBoolean("AlgoSettings", "dynamicDimCoeff", dynamicDimCoeff);
           // fineGrainTuning = cfg.GetBoolean("AlgoSettings", "fineGrainTuning", fineGrainTuning);
            //external_wave = cfg.GetBoolean("AlgoSettings", "external_wave", external_wave);
            
            
            
            
            auto interpAlgoStr = cfg.Get("AlgoSettings", "InterpolationAlgo", "");
            if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_LINEAR]) {
                interpMeta.interpAlgo = INTERP_ALGO_LINEAR;
            } else if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_CUBIC]) {
                interpMeta.interpAlgo = INTERP_ALGO_CUBIC;
            }
            else if (interpAlgoStr == INTERP_ALGO_STR[INTERP_ALGO_QUAD]) {
                interpMeta.interpAlgo = INTERP_ALGO_QUAD;
            }
            QoZ=cfg.GetInteger("AlgoSettings", "QoZ", QoZ);
            interpMeta.interpParadigm = cfg.GetInteger("AlgoSettings", "interpParadigm", interpMeta.interpParadigm);
            interpMeta.cubicSplineType = cfg.GetInteger("AlgoSettings", "cubicSplineType", interpMeta.cubicSplineType);
            interpMeta.interpDirection = cfg.GetInteger("AlgoSettings", "InterpDirection", interpMeta.interpDirection);
            interpMeta.adjInterp = cfg.GetInteger("AlgoSettings", "adjInterp", interpMeta.adjInterp);
            interpBlockSize = cfg.GetInteger("AlgoSettings", "InterpBlockSize", interpBlockSize);
            blockSize = cfg.GetInteger("AlgoSettings", "BlockSize", blockSize);
            quantbinCnt = cfg.GetInteger("AlgoSettings", "QuantizationBinTotal", quantbinCnt);
            maxStep=cfg.GetInteger("AlgoSettings", "maxStep", maxStep);
            sampleBlockSize=cfg.GetInteger("AlgoSettings", "sampleBlockSize", sampleBlockSize);
            levelwisePredictionSelection=cfg.GetInteger("AlgoSettings", "levelwisePredictionSelection", levelwisePredictionSelection);
            //exhaustiveTuning=cfg.GetInteger("AlgoSettings", "exhaustiveTuning", exhaustiveTuning);
            testLorenzo=cfg.GetInteger("AlgoSettings", "testLorenzo", testLorenzo);
            //linearReduce=cfg.GetBoolean("AlgoSettings", "linearReduce", linearReduce);
            //train=cfg.GetBoolean("AlgoSettings", "train", train);
            //useCoeff=cfg.GetBoolean("AlgoSettings", "useCoeff", useCoeff);
            //regSampleStep=cfg.GetInteger("AlgoSettings", "regSampleStep", regSampleStep);
            multiDimInterp=cfg.GetInteger("AlgoSettings", "multiDimInterp", multiDimInterp);
            //mdCrossInterp=cfg.GetInteger("AlgoSettings", "mdCrossInterp", mdCrossInterp);
            profiling=cfg.GetInteger("AlgoSettings", "profiling", profiling);
            SSIMBlockSize=cfg.GetInteger("AlgoSettings", "SSIMBlockSize", SSIMBlockSize);
            fixBlockSize=cfg.GetInteger("AlgoSettings", "fixBlockSize", fixBlockSize);
            
            
            //profilingFix=cfg.GetBoolean("AlgoSettings", "profilingFix",profilingFix);



            
            pdTuningRealComp=cfg.GetBoolean("AlgoSettings", "pdTuningRealComp", pdTuningRealComp);
            pdTuningAbConf=cfg.GetInteger("AlgoSettings", "pdTuningAbConf", pdTuningAbConf);
           // pdAlpha=cfg.GetReal("AlgoSettings", "pdAlpha", pdAlpha);
            //pdBeta=cfg.GetReal("AlgoSettings", "pdBeta", pdBeta);
            //lastPdTuning=cfg.GetBoolean("AlgoSettings", "lastPdTuning", lastPdTuning);
            abList=cfg.GetInteger("AlgoSettings", "abList", abList);
            crossBlock=cfg.GetInteger("AlgoSettings", "crossBlock", crossBlock);
            //sampleBlockSampleBlockSize=cfg.GetInteger("AlgoSettings", "sampleBlockSampleBlockSize", sampleBlockSampleBlockSize);
            peTracking=cfg.GetBoolean("AlgoSettings", "peTracking", peTracking);
           
            //transformation=cfg.GetInteger("AlgoSettings", "transformation", transformation);
           // trimToZero = cfg.GetInteger("AlgoSettings", "trimToZero", trimToZero);
           
            //blockOrder = cfg.GetInteger("AlgoSettings", "blockOrder", blockOrder);
            //coeffTracking = cfg.GetInteger("AlgoSettings", "coeffTracking", coeffTracking);
           
            
            profStride = cfg.GetInteger("AlgoSettings", "profStride", profStride);
           
            regressiveInterp = cfg.GetInteger("AlgoSettings", "regressiveInterp", regressiveInterp);
            naturalSpline = cfg.GetInteger("AlgoSettings", "naturalSpline", naturalSpline );
            adaptiveMultiDimStride = cfg.GetInteger("AlgoSettings", "adaptiveMultiDimStride", adaptiveMultiDimStride);
            fullAdjacentInterp = cfg.GetInteger("AlgoSettings", "fullAdjacentInterp", fullAdjacentInterp);
           // minAnchorLevel = cfg.GetInteger("AlgoSettings", "minAnchorLevel", minAnchorLevel);




        }

        static size_t size_est() {
            
            return sizeof(size_t) * 5 + sizeof(double) * 4 + sizeof(bool) * 5 + sizeof(uint8_t) * 6+ sizeof(int) * 6 + 50; //original SZ3
        }

        size_t save(unsigned char *&c) {
            auto c0 = c;

            write(sz3MagicNumber, c);
            write(sz3DataVer, c);
            write(N, c);

            write(qoz_used, c);
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

            //write(tuningTarget, c);
            //write(absErrorBound, c);
            //write(relErrorBound, c);
            //write(alpha,c);
           // write(beta,c);
            //write(autoTuningRate,c);
            //write(predictorTuningRate,c);
            //write(lorenzo, c);
            //write(lorenzo2, c);
           // write(regression, c);
            //write(regression2, c);
            /*
            write(interpMeta.interpAlgo, c);
            write(interpMeta.interpParadigm, c);
            write(interpMeta.cubicSplineType, c);
            write(interpMeta.interpDirection, c);
            write(interpMeta.adjInterp, c);
            */
            //write(interpMeta,c);
            write(interpBlockSize, c);
            uint8_t boolvals = (lorenzo & 1) << 7 | (lorenzo2 & 1) << 6 | (regression & 1) << 5 | (regression2 & 1) << 4 |
                       (openmp & 1) << 3;
            write(boolvals, c);

            write(dataType, c);
            write(lossless, c);
            write(encoder, c);
            write(quantbinCnt, c);
            write(blockSize, c);
            
           // write(levelwisePredictionSelection, c);
            //write(blockwiseTuning, c);
            write(stride, c);
           // write(maxStep,c);
            write(pred_dim, c);
            //write(openmp, c);
            //write(fixBlockSize, c);
            //write(blockwiseSampleBlockSize, c);
            //write(QoZ, c);//recently changed.
            //write(crossBlock, c);
            return c-c0;
            

            
        }

        void load(const unsigned char *&c) {

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

            read(qoz_used, c);
            //dims.resize(N);
            //read(dims.data(), N, c);
            uint8_t bitWidth;
            read(bitWidth, c);
            dims = bytes2vector<size_t>(c, bitWidth, N);

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
            //read(tuningTarget, c);
            //read(alpha,c);
            //read(beta,c);
            //read(autoTuningRate,c);
            //read(predictorTuningRate,c);
            //read(lorenzo, c);
            //read(lorenzo2, c);
            //read(regression, c);
           // read(regression2, c);
            /*
            read(interpMeta.interpAlgo, c);
            read(interpMeta.interpParadigm, c);
            read(interpMeta.cubicSplineType, c);
            read(interpMeta.interpDirection, c);
            read(interpMeta.adjInterp, c);
            */

            read(interpBlockSize, c);


            uint8_t boolvals;
            read(boolvals, c);
            lorenzo = (boolvals >> 7) & 1;
            lorenzo2 = (boolvals >> 6) & 1;
            regression = (boolvals >> 5) & 1;
            regression2 = (boolvals >> 4) & 1;
            openmp = (boolvals >> 3) & 1;

            //read(interpMeta,c);
           
            read(dataType, c);
            read(lossless, c);
            read(encoder, c);
            read(quantbinCnt, c);
            read(blockSize, c);
            //read(levelwisePredictionSelection, c);
            //read(blockwiseTuning, c);
            read(stride, c);
            //read(maxStep,c);
            read(pred_dim, c);
            //read(openmp, c);
            //read(fixBlockSize, c);
            //read(blockwiseSampleBlockSize, c);
            //read(QoZ, c);//recently changed.
            //read(crossBlock, c);
            
            //read(transformation, c);
            //read(external_wave, c);
            
          
           
            //read(blockOrder, c);
          
            //read(trimToZero, c);
            //read(prewave_absErrorBound, c);
            

            
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
            //printf("InterpolationAlgo = %s\n", enum2Str(static_cast<INTERP_ALGO>(interpAlgo))); //todo: add qoz-related interp information 
            //printf("InterpolationDirection = %d\n", interpDirection);
            printf("QuantizationBinTotal = %d\n", quantbinCnt);
            printf("BlockSize = %d\n", blockSize);
            printf("Stride = %d\n", stride);
            printf("PredDim = %d\n", pred_dim);
            printf("===================== End SZ3 Configuration =====================\n");
        }

        uint32_t sz3MagicNumber = SZ3_MAGIC_NUMBER;
        uint32_t sz3DataVer = versionInt(SZ3_DATA_VER);

        bool qoz_used =false;
        int QoZ = -1;
        char N;
        std::vector<size_t> dims;
        size_t num;
        uint8_t cmprAlgo = ALGO_INTERP_LORENZO;// move to basic
        uint8_t errorBoundMode = EB_ABS;
        uint8_t tuningTarget=TUNING_TARGET_RD;
        double absErrorBound = 1e-3;
        double relErrorBound = 0.0;
        double psnrErrorBound = 0.0;
        double l2normErrorBound = 0.0;
        //double prewave_absErrorBound;
        double rng=-1;
        double alpha=-1;
        double beta=-1;
        double autoTuningRate=0.0;
        double predictorTuningRate=0.0;
        bool lorenzo = true;
        bool lorenzo2 = false;
        bool regression = true;
        bool regression2 = false;
        bool openmp = false; //move to basic
        uint8_t dataType = SZ_FLOAT;  // dataType is only used in HDF5 filter
        uint8_t lossless = 1;         // 0-> skip lossless(use lossless_bypass); 1-> zstd
        uint8_t encoder = 1;          // 0-> skip encoder; 1->HuffmanEncoder; 2->ArithmeticEncoder
        /*
        uint8_t interpAlgo = INTERP_ALGO_CUBIC;
        uint8_t interpParadigm = 0;//1D, MD,HD
        uint8_t cubicSplineType = 0;//noknot,nat
        uint8_t interpDirection = 0;//0,N!-1
        uint8_t adjInterp=0;//0,1
        */
        Interp_Meta interpMeta;
        std::vector <Interp_Meta> interpMeta_list;

        
        int levelwisePredictionSelection=0;
        /*
        std::vector <uint8_t> interpAlgo_list;
        std::vector <uint8_t> interpDirection_list;
        std::vector <uint8_t> cubicSplineType_list;
        */

        
        
        size_t maxStep=0;
        int interpBlockSize = 32;
        int quantbinCnt = 65536;
        int blockSize = 0;
        //int exhaustiveTuning=0;
        int testLorenzo=0;
        std::vector<int> quant_bins;
        //double pred_square_error;
        //double decomp_square_error;
        std::vector<size_t> quant_bin_counts;
        int sampleBlockSize=0;
        bool blockwiseTuning=false;
        int stride = 0; //not used now
        uint8_t pred_dim = 0; // not used now
        //bool linearReduce=0;
        //bool train=0;
        //bool useCoeff=0;
        //int regSampleStep=6;
        int multiDimInterp=0;
        //int mdCrossInterp=0;
        int profiling=0;//since there may be multiple ways of profiling set it to int
        int SSIMBlockSize=8;
        int fixBlockSize=0;
        double blockwiseSampleRate=3.0;
        bool dynamicDimCoeff=false;
        bool freezeDimTest=false;
        int adaptiveMultiDimStride=8;
        //std::vector<double> lorenzo1_coeffs;
        //std::vector<double> lorenzo2_coeffs;
        bool verbose=1;
       
        bool pdTuningRealComp=0;
        int pdTuningAbConf=0;
        //double pdAlpha=-1;
        //double pdBeta=-1;
        //bool lastPdTuning=0;
        int abList=0;
        int crossBlock=0;
        //int sampleBlockSampleBlockSize=0;
        bool peTracking=false;
       //bool external_wave=0;
     
        //vint coeffTracking=0;//0 no. 1: output coeff. 2: print stats of coeff 3: both
        //int pid=0;

        //int transformation = 0; //0: no trans; 1: sigmoid 2: tanh
        //std::vector<float> predictionErrors;//for debug, to delete in final version.
        //std::vector<uint8_t> interp_ops;//for debug, to delete in final version.
        //std::vector<uint8_t> interp_dirs;//for debug, to delete in final version.
        //int trimToZero=0;//1: trim only when quantizing;2: also trim before compression.
        //double preTrim=0.0;//trim small numbers to zero before compression.
        //int blockOrder = 0;//order of blocks.
        double lorenzoBrFix = 1.0;
        bool var_first=false;
        size_t profStride=0;
        
        
       
      
       
        int frozen_dim=-1;
        int regressiveInterp=0;
        int fullAdjacentInterp=0;
        bool naturalSpline=0;
        bool quadInterp=false;


        //bool fineGrainTuning=false;
        //bool profilingFix=true;//only for test

       // double anchorThreshold=0.0;
       // size_t minAnchorLevel=3;


        

    };





}  // namespace SZ3

#endif  // SZ_CONFIG_HPP
