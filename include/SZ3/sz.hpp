#ifndef SZ3_SZ_HPP
#define SZ3_SZ_HPP

#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZMetaFrontend.hpp"
#include "SZ3/frontend/SZGeneralFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/predictor/PolyRegressionPredictor.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Verification.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/szInterp.hpp"
#include <cmath>
#include <memory>


template<class T, uint N, class Quantizer, class Encoder, class Lossless>
std::shared_ptr<SZ::concepts::CompressorInterface < T>>
make_lorenzo_regression_compressor(const SZ::Config &conf, Quantizer quantizer, Encoder encoder, Lossless lossless) {
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;

    int use_single_predictor =
            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return SZ::make_sz_general_compressor<T, N>(
                    SZ::make_sz_general_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return SZ::make_sz_general_compressor<T, N>(
                    SZ::make_sz_general_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return SZ::make_sz_general_compressor<T, N>(
                    SZ::make_sz_general_frontend<T, N>(conf, SZ::RegressionPredictor<T, N>(conf.block_size, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.absErrorBound));
        }
    }

    if (conf.enable_2ndregression) {
        if (use_single_predictor) {
            return SZ::make_sz_general_compressor<T, N>(
                    SZ::make_sz_general_frontend<T, N>(conf, SZ::PolyRegressionPredictor<T, N>(conf.block_size, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<SZ::PolyRegressionPredictor<T, N>>(conf.block_size, conf.absErrorBound));
        }
    }
    return SZ::make_sz_general_compressor<T, N>(
            SZ::make_sz_general_frontend<T, N>(conf, SZ::ComposedPredictor<T, N>(predictors),
                                               quantizer), encoder, lossless);
}

template<class T, uint N>
char *SZ_compress_N(SZ::Config &conf, T *data, size_t &outSize) {

//    std::cout << "****************** Interp Compression ****************" << std::endl;
//    std::cout << "Interp Op          = " << interp_op << std::endl
//              << "Direction          = " << direction_op << std::endl
//              << "SZ block size      = " << block_size << std::endl
//              << "Interp block size  = " << interp_block_size << std::endl;

    assert(N == conf.N);
    if (conf.errorBoundMode != ABS) {
        if (conf.errorBoundMode == REL) {
            conf.errorBoundMode = ABS;
            conf.absErrorBound = conf.relErrorBound * SZ::data_range(data, conf.num);
        } else {
            printf("Error, error bound mode not supported\n");
            exit(0);
        }
    }

    char *cmpData;
    if (conf.cmprMethod == METHOD_LORENZO_REG || conf.cmprMethod == METHOD_LORENZO_REG_FAST) {
        auto quantizer = SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quant_state_num / 2);
        if (N == 3 && conf.cmprMethod == METHOD_LORENZO_REG_FAST) {
            auto sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_meta_frontend<T, N>(conf, quantizer), SZ::HuffmanEncoder<int>(),
                                                           SZ::Lossless_zstd());
            cmpData = (char *) sz->compress(conf, data, outSize);
        } else {
            auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());
            cmpData = (char *) sz->compress(conf, data, outSize);
        }

    } else if (conf.cmprMethod == METHOD_INTERP) {
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(conf.absErrorBound),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd());
        cmpData = (char *) sz.compress(conf, data, outSize);
    }
    cmpData[outSize] = conf.cmprMethod;
    outSize += 1;
    return cmpData;
}


template<class T, uint N>
T *SZ_decompress_N(const SZ::Config &conf, char *cmpData, size_t cmpSize) {
    char cmprMethod = cmpData[cmpSize - 1];
    SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData;

    if (cmprMethod == METHOD_LORENZO_REG || cmprMethod == METHOD_LORENZO_REG_FAST) {
        SZ::LinearQuantizer<T> quantizer;
        if (N == 3 && cmprMethod == METHOD_LORENZO_REG_FAST) {
            auto sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_meta_frontend<T, N>(conf, quantizer),
                                                           SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());
            return sz->decompress(cmpDataPos, cmpSize, conf.num);

        } else {
            auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());
            return sz->decompress(cmpDataPos, cmpSize, conf.num);
        }

    } else if (cmprMethod == METHOD_INTERP) {
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd());
        return sz.decompress(cmpDataPos, cmpSize, conf.num);
    }
    return nullptr;
}

template<class T>
char *SZ_compress(SZ::Config &conf, T *data, size_t &outSize) {
    if (conf.N == 1) {
        return SZ_compress_N<T, 1>(conf, data, outSize);
    } else if (conf.N == 2) {
        return SZ_compress_N<T, 2>(conf, data, outSize);
    } else if (conf.N == 3) {
        return SZ_compress_N<T, 3>(conf, data, outSize);
    } else if (conf.N == 4) {
        return SZ_compress_N<T, 4>(conf, data, outSize);
    } else {
        for (int i = 4; i < conf.N; i++) {
            conf.dims[3] *= conf.dims[i];
        }
        conf.dims.resize(4);
        conf.N = 4;
        return SZ_compress_N<T, 4>(conf, data, outSize);
    }
}

template<class T>
T *SZ_decompress(const SZ::Config &conf, char *cmpData, size_t cmpSize) {
    if (conf.N == 1) {
        return SZ_decompress_N<T, 1>(conf, cmpData, cmpSize);
    } else if (conf.N == 2) {
        return SZ_decompress_N<T, 2>(conf, cmpData, cmpSize);
    } else if (conf.N == 3) {
        return SZ_decompress_N<T, 3>(conf, cmpData, cmpSize);
    } else if (conf.N == 4) {
        return SZ_decompress_N<T, 4>(conf, cmpData, cmpSize);
    } else {
        return SZ_decompress_N<T, 4>(conf, cmpData, cmpSize);
    }
}



#endif