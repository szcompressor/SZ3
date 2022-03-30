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
#include "SZ3/qoi/QoIInfo.hpp"
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

template<class T, SZ::uint N, class Quantizer, class Quantizer_EB>
std::shared_ptr<SZ::concepts::CompressorInterface<T>>
make_qoi_lorenzo_compressor(const SZ::Config &conf, std::shared_ptr<SZ::concepts::QoIInterface<T, N>> qoi, Quantizer quantizer, Quantizer_EB quantizer_eb) {

    quantizer.clear();
    quantizer_eb.clear();
    std::shared_ptr<SZ::concepts::CompressorInterface<T>> sz;

    int methodCnt = (conf.lorenzo + conf.lorenzo2);
    int use_single_predictor = (methodCnt == 1);

    if(use_single_predictor){
        if(conf.lorenzo){
            sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                    SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
        }
        else if(conf.lorenzo2){
            sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                    SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
        }
    }
    else{
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
        predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
        predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
        sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::ComposedPredictor<T, N>(predictors), quantizer, quantizer_eb, qoi),
                                                SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
    }
    return sz;
}

template<class T, SZ::uint N>
char *SZ_compress_LorenzoReg(SZ::Config &conf, T *data, size_t &outSize) {

    assert(N == conf.N);
    assert(conf.cmprAlgo == SZ::ALGO_LORENZO_REG);
    SZ::calAbsErrorBound(conf, data);

    char *cmpData;
    if(conf.qoi > 0){
        //std::cout << "absErrorBound = " << conf.absErrorBound << std::endl;
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << " " << conf.qoiRegionSize << std::endl;
        auto quantizer = SZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = SZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2);
        auto qoi = SZ::GetQOI<T, N>(conf);
        if(conf.qoi == 3){
            conf.blockSize = conf.qoiRegionSize;
        }
        // use sampling to determine abs bound
        {
            auto dims = conf.dims;
            auto tmp_abs_eb = conf.absErrorBound;

            size_t sampling_num, sampling_block;
            std::vector<size_t> sample_dims(N);
            std::vector<T> samples = SZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
            conf.setDims(sample_dims.begin(), sample_dims.end());

            auto sz = make_qoi_lorenzo_compressor(conf, qoi, quantizer, quantizer_eb);
            T * sampling_data = (T *) malloc(sampling_num * sizeof(T));
            // get current ratio
            double ratio = 0;
            {
                size_t sampleOutSize;
                memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                auto cmprData = sz->compress(conf, sampling_data, sampleOutSize);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;                
                std::cout << "current_eb = " << conf.absErrorBound << ", current_ratio = " << ratio << std::endl;
            }
            double prev_ratio = 1;
            double current_ratio = ratio;
            double best_abs_eb = conf.absErrorBound;
            double best_ratio = current_ratio;
            // check smaller bounds
            while(true){
                auto prev_eb = conf.absErrorBound;
                prev_ratio = current_ratio;
                conf.absErrorBound /= 2;
                qoi->set_global_eb(conf.absErrorBound);
                size_t sampleOutSize;
                memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                auto cmprData = sz->compress(conf, sampling_data, sampleOutSize);
                delete[]cmprData;
                current_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;                
                std::cout << "current_eb = " << conf.absErrorBound << ", current_ratio = " << current_ratio << std::endl;
                if(current_ratio < prev_ratio){
                    if(prev_ratio > best_ratio){
                        best_abs_eb = prev_eb;
                        best_ratio = prev_ratio;
                    }
                    break;
                }
            }
            // set error bound
            free(sampling_data);
            //std::cout << "Best abs eb / pre-set eb: " << best_abs_eb / tmp_abs_eb << std::endl; 
            //std::cout << best_abs_eb << " " << tmp_abs_eb << std::endl;
            conf.absErrorBound = best_abs_eb;
            qoi->set_global_eb(best_abs_eb);
            conf.setDims(dims.begin(), dims.end());
        }
        auto sz = make_qoi_lorenzo_compressor(conf, qoi, quantizer, quantizer_eb);
        cmpData = (char *) sz->compress(conf, data, outSize);
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
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << " " << conf.qoiRegionSize << std::endl;
        auto quantizer = SZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = SZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2);
        auto qoi = SZ::GetQOI<T, N>(conf);
        std::shared_ptr<SZ::concepts::CompressorInterface<T>> sz;

        int methodCnt = (conf.lorenzo + conf.lorenzo2);
        int use_single_predictor = (methodCnt == 1);

        if(use_single_predictor){
            if(conf.lorenzo){
                sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                        SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
            }
            else if(conf.lorenzo2){
                sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                        SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
            }
        }
        else{
            std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
            sz = SZ::make_sz_general_compressor<T, N>(SZ::make_sz_qoi_frontend<T, N>(conf, SZ::ComposedPredictor<T, N>(predictors), quantizer, quantizer_eb, qoi),
                                                    SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
        }
        sz->decompress(cmpDataPos, cmpSize, decData);
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
