#ifndef SZ3_SZINTERP_HPP
#define SZ3_SZINTERP_HPP

#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZMetaFrontend.hpp"
#include "SZ3/frontend/SZGeneralFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/predictor/SimplePredictor.hpp"
#include "szLorenzoReg.hpp"
#include <cmath>
#include <memory>


template<class T, uint N>
char *SZ_compress_Interp_N(SZ::Config &conf, T *data, size_t &outSize) {

//    std::cout << "****************** Interp Compression ****************" << std::endl;
//    std::cout << "Interp Op          = " << interp_op << std::endl
//              << "Direction          = " << direction_op << std::endl
//              << "SZ block size      = " << block_size << std::endl
//              << "Interp block size  = " << interp_block_size << std::endl;

    assert(N == conf.N);
    SZ::calAbsErrorBound(conf, data);

    char *cmpData;
    if (conf.cmprMethod == METHOD_INTERP) {
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(conf.absErrorBound),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd());
        cmpData = (char *) sz.compress(conf, data, outSize);
    }
    return cmpData;
}


template<class T, uint N>
T *SZ_decompress_Interp_N(const SZ::Config &conf, char *cmpData, size_t cmpSize) {
    SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData;
    if (conf.cmprMethod == METHOD_INTERP) {
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd());
        return sz.decompress(cmpDataPos, cmpSize, conf.num);
    }
    return nullptr;
}


template<class T, uint N>
double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num, double eb, int interp_level,
                                                  int interp_op, int direction_op,
                                                  int block_size, int interp_block_size) {

    std::cout << "****************** Interp Compression ****************" << std::endl;
    std::cout << "Interp Op          = " << interp_op << std::endl
              << "Direction          = " << direction_op << std::endl
              << "SZ block size      = " << block_size << std::endl
              << "Interp block size  = " << interp_block_size << std::endl;
//              << "SZ_mode            = " << sz_op << std::endl;

    SZ::Timer timer(true);

    std::vector<T> data1(data, data + num);
    size_t compressed_size = 0;

    SZ::Config conf;
    conf.absErrorBound = eb;
    conf.update_dims(dims.begin(), dims.end());
    conf.block_size = block_size;
    conf.stride = conf.block_size;
    auto sz = SZ::make_sz_block_interpolation_compressor<T, N>(
            conf,
            SZ::SimplePredictor<T, N>(eb),
            SZ::LinearQuantizer<T>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd(),
            interp_op,
            direction_op,
            interp_level
    );

    auto cmpData = sz.compress(data1.data(), compressed_size);
    delete[]cmpData;

    double compression_time = timer.stop();

    auto compression_ratio = num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compression_ratio << std::endl;
    std::cout << "Interp compression time = " << compression_time
              << " Ratio = " << compression_ratio
              << " Params = " << interp_level
              << " " << interp_op
              << " " << direction_op
              << " " << block_size
              << " " << interp_block_size
              //              << " " << sz_op
              << std::endl;


    std::cout << "****************** Interp end ****************" << std::endl;
    return compression_ratio;
}

template<class T, uint N>
char *SZ_compress_Interp_lorenzo_N(SZ::Config &conf, T *data, size_t &outSize) {
    std::cout << "================================ BEGIN TUNING ================================" << std::endl;
    SZ::Timer timer(true);

    SZ::calAbsErrorBound(conf, data);

    size_t sampling_num, sampling_block;
    std::vector<size_t> sample_dims(N);
    std::vector<T> sampling_data = SZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
//    printf("%lu %lu %lu %lu %lu\n", sampling_data.size(), sampling_num, sample_dims[0], sample_dims[1], sample_dims[2]);

    SZ::Config lorenzo_config = conf;
    lorenzo_config.cmprMethod = METHOD_LORENZO_REG_FAST;
    lorenzo_config.update_dims(sample_dims.begin(), sample_dims.end());
    lorenzo_config.enable_lorenzo = true;
    lorenzo_config.enable_2ndlorenzo = true;
    lorenzo_config.enable_regression = false;
    lorenzo_config.block_size = 5;
    lorenzo_config.quant_state_num = 65536 * 2;
    size_t sampleOutSize;
    auto cmprData = SZ_compress_LorenzoReg_N<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
    delete[]cmprData;
    double ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
    printf("Lorenzo ratio = %.2f\n", ratio);

    double best_lorenzo_ratio = ratio, best_interp_ratio = 0;
    int interp_level = -1, interp_op, direction_op = 0, block_size = sampling_block, interp_block_size = sampling_block;

    {
        //tune interp
        for (int i = 0; i < 2; i++) {
            ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                     interp_level,
                                                                     i, direction_op,
                                                                     block_size, interp_block_size);
            if (ratio > best_interp_ratio) {
                best_interp_ratio = ratio;
                interp_op = i;
            }
        }
        std::cout << "interp best interp_op = " << interp_op << " , best ratio = " << best_interp_ratio << std::endl;

        ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound, interp_level,
                                                                 interp_op, 5,
                                                                 block_size, interp_block_size);
        if (ratio > best_interp_ratio * 1.02) {
            best_interp_ratio = ratio;
            direction_op = 5;
        }
        std::cout << "interp best direction_op = " << direction_op << " , best ratio = " << best_interp_ratio
                  << std::endl;
    }

    bool useInterp = !(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80);
    printf("\nLorenzo compression ratio = %.2f\n", best_lorenzo_ratio);
    printf("Interp compression ratio = %.2f\n", best_interp_ratio);
    printf("choose %s\n", useInterp ? "Interp" : "Lorenzo");

    if (useInterp) {
        conf.cmprMethod = METHOD_INTERP;
        conf.interp_op = interp_op;
        conf.interp_direction_op = direction_op;
        double tuning_time = timer.stop();
        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        std::cout << "====================================== END TUNING ======================================"
                  << std::endl << std::endl;
        return SZ_compress_Interp_N<T, N>(conf, data, outSize);
    } else {
        //further tune lorenzo
        if (N == 3) {
            lorenzo_config.quant_state_num = SZ::optimize_quant_invl_3d<T>(data, conf.dims[0], conf.dims[1], conf.dims[2], conf.absErrorBound);
            lorenzo_config.pred_dim = 2;
            cmprData = SZ_compress_LorenzoReg_N<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
            delete[]cmprData;
            ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
            printf("Lorenzo, pred_dim=2, ratio = %.2f\n", ratio);
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.pred_dim = 3;
            }
        }

        if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5) {
            auto quant_num = lorenzo_config.quant_state_num;
            lorenzo_config.quant_state_num = 16384;
            cmprData = SZ_compress_LorenzoReg_N<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
            delete[]cmprData;
            ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
            printf("Lorenzo, quant_bin=8192, ratio = %.2f\n", ratio);
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.quant_state_num = quant_num;
            }
        }
        lorenzo_config.update_dims(conf.dims.begin(), conf.dims.end());
        conf = lorenzo_config;
        double tuning_time = timer.stop();
        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        std::cout << "====================================== END TUNING ======================================"
                  << std::endl << std::endl;
        return SZ_compress_LorenzoReg_N<T, N>(conf, data, outSize);
    }


}

#endif