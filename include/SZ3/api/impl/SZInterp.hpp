#ifndef SZ3_SZINTERP_HPP
#define SZ3_SZINTERP_HPP

#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/SZQoIInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/quantizer/QoIIntegerQuantizer.hpp"
#include "SZ3/encoder/QoIEncoder.hpp"
#include "SZ3/qoi/QoIInfo.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/api/impl/SZLorenzoReg.hpp"
#include <cmath>
#include <memory>

template<class T, SZ::uint N>
char *SZ_compress_Interp(SZ::Config &conf, T *data, size_t &outSize) {

//    std::cout << "****************** Interp Compression ****************" << std::endl;
//    std::cout << "Interp Op          = " << interpAlgo << std::endl
//              << "Direction          = " << direction << std::endl
//              << "SZ block size      = " << blockSize << std::endl
//              << "Interp block size  = " << interpBlockSize << std::endl;

    assert(N == conf.N);
    assert(conf.cmprAlgo == SZ::ALGO_INTERP);
    SZ::calAbsErrorBound(conf, data);

    // conf.print();
    // directly use abs when qoi is regional average
    if(conf.qoi > 0){
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;
        auto quantizer = SZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = SZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2);
        auto qoi = SZ::GetQOI<T, N>(conf);
        auto sz = SZ::SZQoIInterpolationCompressor<T, N, SZ::VariableEBLinearQuantizer<T, T>, SZ::EBLogQuantizer<T>, SZ::QoIEncoder<int>, SZ::Lossless_zstd>(
                quantizer, quantizer_eb, qoi, SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
        double max_abs_eb = 0;
        // use sampling to determine abs bound
        {
            auto dims = conf.dims;
            auto tmp_abs_eb = conf.absErrorBound;

            size_t sampling_num, sampling_block;
            std::vector<size_t> sample_dims(N);
            std::vector<T> samples = SZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
            conf.setDims(sample_dims.begin(), sample_dims.end());

            T * sampling_data = (T *) malloc(sampling_num * sizeof(T));
            // reset dimensions for average of square
            if(conf.qoi == 3) qoi->set_dims(sample_dims);
            // get current ratio
            double ratio = 0;
            {
                size_t sampleOutSize;
                memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                // reset variables for average of square
                if(conf.qoi == 3) qoi->init();
                auto cmprData = sz.compress(conf, sampling_data, sampleOutSize);
                max_abs_eb = sz.get_max_eb();
                sz.clear();
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;                
                std::cout << "current_eb = " << conf.absErrorBound << ", current_ratio = " << ratio << std::endl;
            }
            double prev_ratio = 1;
            double current_ratio = ratio;
            double best_abs_eb = std::min(conf.absErrorBound, max_abs_eb);
            double best_ratio = current_ratio;
            // check smaller bounds
            int max_iter = 100; 
            int iter = 0;
            while(iter++ < max_iter){
                auto prev_eb = conf.absErrorBound;
                prev_ratio = current_ratio;
                conf.absErrorBound /= 2;
                qoi->set_global_eb(conf.absErrorBound);
                size_t sampleOutSize;
                memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                // reset variables for average of square
                if(conf.qoi == 3) qoi->init();
                auto cmprData = sz.compress(conf, sampling_data, sampleOutSize);
                sz.clear();
                delete[]cmprData;
                current_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;                
                std::cout << "current_eb = " << conf.absErrorBound << ", current_ratio = " << current_ratio << std::endl;
                if(current_ratio < prev_ratio * 0.99){
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
            // reset dimensions and variables for average of square
            if(conf.qoi == 3){
                qoi->set_dims(dims);
                qoi->init();
            }
        }
        char *cmpData = (char *) sz.compress(conf, data, outSize);
        return cmpData;
    }
    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(conf.absErrorBound),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    char *cmpData = (char *) sz.compress(conf, data, outSize);
    return cmpData;
}


template<class T, SZ::uint N>
void SZ_decompress_Interp(const SZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
    //std::cout << "SZ_decompress_Interp" << std::endl;
    assert(conf.cmprAlgo == SZ::ALGO_INTERP);
    SZ::uchar const *cmpDataPos = (SZ::uchar *) cmpData;
    // directly use abs when qoi is regional average
    if(conf.qoi > 0){
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;
        auto quantizer = SZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = SZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2);
        auto qoi = SZ::GetQOI<T, N>(conf);
        auto sz = SZ::SZQoIInterpolationCompressor<T, N, SZ::VariableEBLinearQuantizer<T, T>, SZ::EBLogQuantizer<T>, SZ::QoIEncoder<int>, SZ::Lossless_zstd>(
                quantizer, quantizer_eb, qoi, SZ::QoIEncoder<int>(), SZ::Lossless_zstd());
        sz.decompress(cmpDataPos, cmpSize, decData);
        return;
    }   
    auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    sz.decompress(cmpDataPos, cmpSize, decData);
}


template<class T, SZ::uint N>
double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num,
                                                  double eb, int interp_op, int direction_op, int block_size) {

    std::vector<T> data1(data, data + num);
    size_t outSize = 0;

    SZ::Config conf;
    conf.absErrorBound = eb;
    conf.setDims(dims.begin(), dims.end());
    conf.blockSize = block_size;
    conf.interpAlgo = interp_op;
    conf.interpDirection = direction_op;
    auto sz = SZ::SZBlockInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<T>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    char *cmpData = (char *) sz.compress(conf, data1.data(), outSize);
    delete[]cmpData;
    auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
    return compression_ratio;
}

template<class T, SZ::uint N>
char *SZ_compress_Interp_lorenzo(SZ::Config &conf, T *data, size_t &outSize) {
    assert(conf.cmprAlgo == SZ::ALGO_INTERP_LORENZO);

    //std::cout << "====================================== BEGIN TUNING ================================" << std::endl;
    SZ::Timer timer(true);

    SZ::calAbsErrorBound(conf, data);
    // overwrite qoi for parameter exploration
    int qoi = conf.qoi;
    auto tmp_abs_eb = conf.absErrorBound;
    if(qoi){
        // compute abs qoi eb
        T qoi_rel_eb = conf.qoiEB;
        std::cout << "qoi_rel_eb = " << qoi_rel_eb << std::endl;
        T max = data[0];
        T min = data[0];
        double min_abs ;
        double max_abs ;
        double min_abs_2;
        for (size_t i = 1; i < conf.num; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }
        if(qoi == 1 || qoi == 3){
            // x^2
            auto max_2 = max * max;
            auto min_2 = min * min;
            double max_abs_val = (max_2 > min_2) ? max_2 : min_2;
            double min_abs_val = (max_2 > min_2) ? min_2 : max_2;
            if(min< 0) min_abs_2 = 0;
            else min_abs_2 = min_abs_val;
            conf.qoiEB *= (max_abs_val - min_abs_2);
            // conf.qoiEB *= (max_abs_val);

            
        }
        // else if(qoi == 3){
        //     // regional average
        //     conf.qoiEB *= max - min;
        //     conf.absErrorBound = conf.qoiEB;
        // }
        else if(qoi == 4){
            // compute isovalues
            conf.isovalues.clear();
            int isonum = conf.qoiIsoNum;
            auto range = max - min;
            for(int i=0; i<isonum; i++){
                conf.isovalues.push_back(min + (i + 1) * 1.0 / (isonum + 1) * range);
            }
        }
        else if(qoi >= 5){
            // (x^2) + (log x) + (isoline)
            auto max_2 = max * max;
            auto min_2 = min * min;
            auto max_abs_val = (max_2 > min_2) ? max_2 : min_2;
            auto eb_x2 = conf.qoiEB * max_abs_val;
            auto eb_logx = conf.qoiEB * 10;
            auto eb_isoline = 0;
            conf.isovalues.clear();
            int isonum = conf.qoiIsoNum;
            auto range = max - min;
            for(int i=0; i<isonum; i++){
                conf.isovalues.push_back(min + (i + 1) * 1.0 / (isonum + 1) * range);
            }
            if(qoi == 5){
                // x^2 + log x
                conf.qoiEBs.push_back(eb_x2);
                conf.qoiEBs.push_back(eb_logx);
            }
            else if(qoi == 6){
                // x^2 + isoline
                conf.qoiEBs.push_back(eb_x2);
                conf.qoiEBs.push_back(eb_isoline);                
            }
            else if(qoi == 7){
                // log x + isoline
                conf.qoiEBs.push_back(eb_logx);
                conf.qoiEBs.push_back(eb_isoline);                
            }
            else if(qoi == 8){
                // x^2 + log x + isoline
                conf.qoiEBs.push_back(eb_x2);
                conf.qoiEBs.push_back(eb_logx);
                conf.qoiEBs.push_back(eb_isoline);                                
            }
        }
        // set eb base and log base if not set by config
        if(conf.qoiEBBase == 0) 
            conf.qoiEBBase = std::numeric_limits<T>::epsilon();
        if(conf.qoiEBLogBase == 0)
            conf.qoiEBLogBase = 2;        
        // update eb base
        if(qoi != 4 && qoi != 7) conf.qoiEBBase = (max - min) * qoi_rel_eb / 1030;
        std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;
    }
    else{
        // compute isovalues for comparison
        T max = data[0];
        T min = data[0];
        for (size_t i = 1; i < conf.num; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }
        conf.isovalues.clear();
        int num = conf.qoiIsoNum;
        auto range = max - min;
        for(int i=0; i<num; i++){
            conf.isovalues.push_back(min + range / (num + 1));
        }        
    }
    conf.qoi = 0;

    size_t sampling_num, sampling_block;
    std::vector<size_t> sample_dims(N);
    std::vector<T> sampling_data = SZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
    // printf("%lu %lu %lu %lu %lu\n", sampling_data.size(), sampling_num, sample_dims[0], sample_dims[1], sample_dims[2]);

    SZ::Config lorenzo_config = conf;
    lorenzo_config.cmprAlgo = SZ::ALGO_LORENZO_REG;
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
        for (auto &interp_op: {SZ::INTERP_ALGO_LINEAR, SZ::INTERP_ALGO_CUBIC}) {
            ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                     interp_op, conf.interpDirection, sampling_block);
            if (ratio > best_interp_ratio) {
                best_interp_ratio = ratio;
                conf.interpAlgo = interp_op;
            }
        }
        //std::cout << "interp best interpAlgo = " << (conf.interpAlgo == 0 ? "LINEAR" : "CUBIC") << std::endl;

        int direction_op = SZ::factorial(N) - 1;
        ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                 conf.interpAlgo, direction_op, sampling_block);
        if (ratio > best_interp_ratio * 1.02) {
            best_interp_ratio = ratio;
            conf.interpDirection = direction_op;
        }
        //std::cout << "interp best direction = " << (unsigned) conf.interpDirection << std::endl;

    }

    bool useInterp = !(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80);
//    printf("\nLorenzo compression ratio = %.2f\n", best_lorenzo_ratio);
//    printf("Interp compression ratio = %.2f\n", best_interp_ratio);
    printf("choose %s\n", useInterp ? "interp" : "Lorenzo");

    if (useInterp) {
        conf.cmprAlgo = SZ::ALGO_INTERP;
        double tuning_time = timer.stop();
//        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        //std::cout << "====================================== END TUNING ======================================" << std::endl;
        // assign qoi back
        conf.qoi = qoi;
        conf.lorenzo = false;
        conf.lorenzo2 = false;
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
        // if(qoi){
        //     // if qoi is enabled, sample between lorenzo and lorenzo2
        //     lorenzo_config.lorenzo = true;
        //     lorenzo_config.lorenzo2 = false;
        //     cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
        //     delete[]cmprData;
        //     lorenzo_config.lorenzo = false;
        //     lorenzo_config.lorenzo2 = true;
        //     size_t sampleOutSize2 = 0;
        //     cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize2);
        //     delete[]cmprData;
        //     printf("lorenzo size = %lu, lorenzo2 size = %lu\n", sampleOutSize, sampleOutSize2);
        //     if(sampleOutSize < sampleOutSize2){
        //         // lorenzo is better
        //         lorenzo_config.lorenzo = true;
        //         lorenzo_config.lorenzo2 = false;
        //         printf("choose lorenzo\n");
        //     }
        //     else{
        //         // lorenzo2 is better
        //         printf("choose lorenzo2\n");
        //     }
        // }
        lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
        conf = lorenzo_config;
        // assign qoi back
        conf.qoi = qoi;
        double tuning_time = timer.stop();
//        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        // printf("lorenzo = %d, lorenzo2 = %d, block size = %d\n", conf.lorenzo, conf.lorenzo2, conf.blockSize);
        //std::cout << "====================================== END TUNING ======================================" << std::endl;
        return SZ_compress_LorenzoReg<T, N>(conf, data, outSize);
    }


}

#endif
