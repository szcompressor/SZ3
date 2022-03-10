#ifndef SZ3_SZ_LORENZO_REG_HPP
#define SZ3_SZ_LORENZO_REG_HPP

#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZFastFrontend.hpp"
#include "SZ3/frontend/SZGeneralFrontend.hpp"
#include "SZ3/frontend/SZQoIFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/quantizer/QoIIntegerQuantizer.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/predictor/PolyRegressionPredictor.hpp"
#include "SZ3/encoder/QoIEncoder.hpp"
#include "SZ3/qoi/XSquare.hpp"
#include "SZ3/qoi/LogX.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/def.hpp"
#include <cmath>
#include <memory>


template<class T, SZ::uint N, class Quantizer, class Encoder, class Lossless>
std::shared_ptr<SZ::concepts::CompressorInterface<T>>
make_lorenzo_regression_compressor(const SZ::Config &conf, Quantizer quantizer, Encoder encoder, Lossless lossless) {
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;

    int methodCnt = (conf.lorenzo + conf.lorenzo2 + conf.regression + conf.regression2);
    int use_single_predictor = (methodCnt == 1);
    if (methodCnt == 0) {
        printf("All lorenzo and regression methods are disabled.\n");
        exit(0);
    }
    if (conf.lorenzo) {
        if (use_single_predictor) {
            return SZ::make_sz_general_compressor<T, N>(
                    SZ::make_sz_general_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
        }
    }
    if (conf.lorenzo2) {
        if (use_single_predictor) {
            return SZ::make_sz_general_compressor<T, N>(
                    SZ::make_sz_general_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
        }
    }
    if (conf.regression) {
        if (use_single_predictor) {
            return SZ::make_sz_general_compressor<T, N>(
                    SZ::make_sz_general_frontend<T, N>(conf, SZ::RegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }

    if (conf.regression2) {
        if (use_single_predictor) {
            return SZ::make_sz_general_compressor<T, N>(
                    SZ::make_sz_general_frontend<T, N>(conf, SZ::PolyRegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<SZ::PolyRegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }
    return SZ::make_sz_general_compressor<T, N>(
            SZ::make_sz_general_frontend<T, N>(conf, SZ::ComposedPredictor<T, N>(predictors),
                                               quantizer), encoder, lossless);
}

void qoi_check(int a, int b){
    if(a != b){
        std::cerr << "QoI number does not match" << std::endl;
        exit(-1);
    }
}

template<class T, SZ::uint N>
char *SZ_compress_LorenzoReg(SZ::Config &conf, T *data, size_t &outSize) {

    assert(N == conf.N);
    assert(conf.cmprAlgo == SZ::ALGO_LORENZO_REG);
    SZ::calAbsErrorBound(conf, data);

    char *cmpData;
    if(conf.qoi > 0){
        std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;
        auto quantizer = SZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = SZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2);
        // text x^2
        // qoi_check(conf.qoi, 1);
        // auto qoi = SZ::QoI_X_Square<T>(conf.qoiEB, conf.num, conf.absErrorBound);
        // test log(x)
        qoi_check(conf.qoi, 2);
        auto qoi = SZ::QoI_Log_X<T>(conf.qoiEB, conf.absErrorBound);
        // will not have reg since SZ3 is used
        assert(conf.regression + conf.regression2 == 0);
        // will use both two Lorenzo predictors
        assert(conf.lorenzo);
        assert(conf.lorenzo2);
        // TODO: identify which one to use
        if(conf.lorenzo){
            auto sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                        SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
            cmpData = (char *) sz->compress(conf, data, outSize);
            conf.lorenzo2 = false;
        }
        else{
            auto sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                        SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
            cmpData = (char *) sz->compress(conf, data, outSize);            
            conf.lorenzo = false;
        }
        return cmpData;
    }
    auto quantizer = SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
    if (N == 3 && !conf.regression2) {
        // use fast version for 3D
        auto sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_fast_frontend<T, N>(conf, quantizer), SZ::HuffmanEncoder<int>(),
                                                       SZ::Lossless_zstd());
        cmpData = (char *) sz->compress(conf, data, outSize);
    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());
        cmpData = (char *) sz->compress(conf, data, outSize);
    }
    return cmpData;
}


template<class T, SZ::uint N>
void SZ_decompress_LorenzoReg(const SZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == SZ::ALGO_LORENZO_REG);

    SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData;
    if(conf.qoi > 0){
        std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;
        auto quantizer = SZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = SZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2);
        // text x^2
        // qoi_check(conf.qoi, 1);
        // auto qoi = SZ::QoI_X_Square<T>(conf.qoiEB, conf.num, conf.absErrorBound);
        // test log(x)
        // qoi_check(conf.qoi, 2);
        auto qoi = SZ::QoI_Log_X<T>(conf.qoiEB, conf.absErrorBound);
        // identify which one to use
        if(conf.lorenzo){
            auto sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                        SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
            sz->decompress(cmpDataPos, cmpSize, decData);
        }
        else{
            auto sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                        SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
            sz->decompress(cmpDataPos, cmpSize, decData);           
        }
        return;
    }    
    SZ::LinearQuantizer<T> quantizer;
    if (N == 3 && !conf.regression2) {
        // use fast version for 3D
        auto sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_fast_frontend<T, N>(conf, quantizer),
                                                       SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());
        sz->decompress(cmpDataPos, cmpSize, decData);
        return;

    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());
        sz->decompress(cmpDataPos, cmpSize, decData);
        return;
    }

}

#endif