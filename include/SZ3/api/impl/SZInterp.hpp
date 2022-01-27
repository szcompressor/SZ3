#ifndef SZ3_SZINTERP_HPP
#define SZ3_SZINTERP_HPP

#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZFastFrontend.hpp"
#include "SZ3/frontend/SZGeneralFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/predictor/SimplePredictor.hpp"
#include "SZ3/api/impl/SZLorenzoReg.hpp"
#include <cmath>
#include <memory>


template<class T, uint N>
char *SZ_compress_Interp(SZ::Config &conf, T *data, size_t &outSize) {

//    std::cout << "****************** Interp Compression ****************" << std::endl;
//    std::cout << "Interp Op          = " << interpAlgo << std::endl
//              << "Direction          = " << direction << std::endl
//              << "SZ block size      = " << blockSize << std::endl
//              << "Interp block size  = " << interpBlockSize << std::endl;

    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_INTERP);
    SZ::calAbsErrorBound(conf, data);

    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(conf.absErrorBound),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    char *cmpData = (char *) sz.compress(conf, data, outSize);
    return cmpData;
}


template<class T, uint N>
void SZ_decompress_Interp(const SZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_INTERP);
    SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData;
    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    sz.decompress(cmpDataPos, cmpSize, decData);
}


template<class T, uint N>
double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num,
                                                  double eb, int interp_op, int direction_op, int block_size) {

    std::vector<T> data1(data, data + num);
    size_t compressed_size = 0;

    SZ::Config conf;
    conf.absErrorBound = eb;
    conf.setDims(dims.begin(), dims.end());
    conf.blockSize = block_size;
    conf.stride = conf.blockSize;
    auto sz = SZ::make_sz_block_interpolation_compressor<T, N>(
            conf,
            SZ::SimplePredictor<T, N>(eb),
            SZ::LinearQuantizer<T>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd(),
            interp_op,
            direction_op,
            -1
    );

    auto cmpData = sz.compress(data1.data(), compressed_size);
    delete[]cmpData;
    auto compression_ratio = num * sizeof(T) * 1.0 / compressed_size;
    return compression_ratio;
}

template<class T, uint N>
char *SZ_compress_Interp_lorenzo(SZ::Config &conf, T *data, size_t &outSize) {
    assert(conf.cmprAlgo == ALGO_INTERP_LORENZO);

    std::cout << "====================================== BEGIN TUNING ================================" << std::endl;
    SZ::Timer timer(true);

    SZ::calAbsErrorBound(conf, data);

    size_t sampling_num, sampling_block;
    std::vector<size_t> sample_dims(N);
    std::vector<T> sampling_data = SZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
//    printf("%lu %lu %lu %lu %lu\n", sampling_data.size(), sampling_num, sample_dims[0], sample_dims[1], sample_dims[2]);

    SZ::Config lorenzo_config = conf;
    lorenzo_config.cmprAlgo = ALGO_LORENZO_REG;
    lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
    lorenzo_config.lorenzo = true;
    lorenzo_config.lorenzo2 = true;
    lorenzo_config.regression = false;
    lorenzo_config.regression2 = false;
    lorenzo_config.openmp = false;
    lorenzo_config.blockSize = 5;
    lorenzo_config.quantbinCnt = 65536 * 2;
    size_t sampleOutSize;
    auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
    delete[]cmprData;
    double ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
//    printf("Lorenzo ratio = %.2f\n", ratio);

    double best_lorenzo_ratio = ratio, best_interp_ratio = 0;

    {
        //tune interp
        for (auto &interp_op: {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC}) {
            ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                     interp_op, conf.interpDirection, sampling_block);
            if (ratio > best_interp_ratio) {
                best_interp_ratio = ratio;
                conf.interpAlgo = interp_op;
            }
        }
        std::cout << "interp best interpAlgo = " << (conf.interpAlgo == 0 ? "LINEAR" : "CUBIC") << std::endl;

        int direction_op = SZ::factorial(N) - 1;
        ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                 conf.interpAlgo, direction_op, sampling_block);
        if (ratio > best_interp_ratio * 1.02) {
            best_interp_ratio = ratio;
            conf.interpDirection = direction_op;
        }
        std::cout << "interp best direction = " << (unsigned) conf.interpDirection << std::endl;
    }

    bool useInterp = !(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80);
//    printf("\nLorenzo compression ratio = %.2f\n", best_lorenzo_ratio);
//    printf("Interp compression ratio = %.2f\n", best_interp_ratio);
    printf("choose %s\n", useInterp ? "interp" : "Lorenzo");

    if (useInterp) {
        conf.cmprAlgo = ALGO_INTERP;
        double tuning_time = timer.stop();
//        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        std::cout << "====================================== END TUNING ======================================" << std::endl;
        return SZ_compress_Interp<T, N>(conf, data, outSize);
    } else {
        //further tune lorenzo
        if (N == 3) {
            lorenzo_config.quantbinCnt = SZ::optimize_quant_invl_3d<T>(data, conf.dims[0], conf.dims[1], conf.dims[2], conf.absErrorBound);
            lorenzo_config.pred_dim = 2;
            cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
            delete[]cmprData;
            ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
//            printf("Lorenzo, pred_dim=2, ratio = %.2f\n", ratio);
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.pred_dim = 3;
            }
        }

        if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5) {
            auto quant_num = lorenzo_config.quantbinCnt;
            lorenzo_config.quantbinCnt = 16384;
            cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
            delete[]cmprData;
            ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
//            printf("Lorenzo, quant_bin=8192, ratio = %.2f\n", ratio);
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.quantbinCnt = quant_num;
            }
        }
        lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
        conf = lorenzo_config;
        double tuning_time = timer.stop();
//        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        std::cout << "====================================== END TUNING ======================================" << std::endl;
        return SZ_compress_LorenzoReg<T, N>(conf, data, outSize);
    }


}

#endif