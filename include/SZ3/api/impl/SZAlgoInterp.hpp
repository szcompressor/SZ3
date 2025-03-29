#ifndef SZ3_SZALGO_INTERP_HPP
#define SZ3_SZALGO_INTERP_HPP

#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
#include "SZ3/compressor/specialized/SZBlockInterpolationCompressor.hpp"
#include "SZ3/decomposition/InterpolationDecomposition.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_Interp(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_INTERP);
    calAbsErrorBound(conf, data);
    if (conf.interp_anchorStride <= 0){
        std::array<size_t, 4> anchor_strides = {512,64,32,16};
        conf.interp_anchorStride = anchor_strides[N - 1];
    }

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    return sz->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_Interp(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_INTERP);
    auto cmpDataPos = cmpData;
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

template <class T, uint N>
double interp_compress_test(T *data, const Config &theConf, std::vector<size_t> dims, size_t num, double eb,
                                                  int interp_op, int direction_op, int block_size, uchar *buffer,
                                                  size_t bufferCap) {
    std::vector<T> data1(data, data + num);

    Config conf = theConf;
    conf.absErrorBound = eb;
    conf.setDims(dims.begin(), dims.end());
    conf.blockSize = block_size;
    conf.interpAlgo = interp_op;
    conf.interpDirection = direction_op;
    conf.tuning = true;
    conf.interp_anchorStride = 0;
    
    auto sz = SZBlockInterpolationCompressor<T, N, LinearQuantizer<T>, HuffmanEncoder<int>, Lossless_zstd>(
        LinearQuantizer<T>(eb), HuffmanEncoder<int>(), Lossless_zstd());
    size_t outSize = sz.compress(conf, data1.data(), buffer, bufferCap);
    
    /*
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    */
    size_t outSize = sz->compress(conf, data1.data(), buffer, bufferCap);
    
    auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
    return compression_ratio;
}

template <class T, uint N>
size_t SZ_compress_Interp_lorenzo(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(conf.cmprAlgo == ALGO_INTERP_LORENZO);

    //        Timer timer(true);
    calAbsErrorBound(conf, data);

    if (conf.interp_anchorStride <= 0){
        std::array<size_t, 4> anchor_strides = {4096,128,32,16};
        conf.interp_anchorStride = anchor_strides[N - 1];
    }

    size_t sampling_num, sampling_block;
    std::vector<size_t> sample_dims(N);
    std::vector<T> sampling_data = sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
    if (sampling_num == conf.num) {
        conf.cmprAlgo = ALGO_INTERP;
        return SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
    }
    double best_lorenzo_ratio = 0, best_interp_ratio = 0, ratio;
    size_t bufferCap = conf.num * sizeof(T);
    auto buffer = static_cast<uchar *>(malloc(bufferCap));
    Config lorenzo_config = conf;
    {
        if(N < 4){
            // test lorenzo
            lorenzo_config.cmprAlgo = ALGO_LORENZO_REG;
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
            lorenzo_config.lorenzo = true;
            lorenzo_config.lorenzo2 = true;
            lorenzo_config.regression = false;
            lorenzo_config.regression2 = false;
            lorenzo_config.openmp = false;
            lorenzo_config.blockSize = 5;
            //        lorenzo_config.quantbinCnt = 65536 * 2;
            std::vector<T> data1(sampling_data);
            size_t sampleOutSize = SZ_compress_LorenzoReg<T, N>(lorenzo_config, data1.data(), buffer, bufferCap);
            //            delete[]cmprData;
            //    printf("Lorenzo ratio = %.2f\n", ratio);
            best_lorenzo_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
        }
    }
    {
        // tune interp
        for (auto &interp_op : {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC, INTERP_ALGO_CUBIC_NATURAL}) {
            ratio = interp_compress_test<T, N>(
                sampling_data.data(), conf, sample_dims, sampling_num, conf.absErrorBound, interp_op, conf.interpDirection,
                sampling_block, buffer, bufferCap);
            if (ratio > best_interp_ratio) {
                best_interp_ratio = ratio;
                conf.interpAlgo = interp_op;
            }
        }

        int direction_op = factorial(N) - 1;
        ratio = interp_compress_test<T, N>(sampling_data.data(), conf, sample_dims, sampling_num,
                                                                 conf.absErrorBound, conf.interpAlgo, direction_op,
                                                                 sampling_block, buffer, bufferCap);
        if (ratio > best_interp_ratio * 1.02) {
            best_interp_ratio = ratio;
            conf.interpDirection = direction_op;
        }
    }
    bool useInterp = !(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80);
    std::cout<<best_lorenzo_ratio<<" "<<best_interp_ratio<<std::endl;
    size_t cmpSize = 0;
    if (useInterp) {
        conf.cmprAlgo = ALGO_INTERP;
        cmpSize = SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
    } else {
        // further tune lorenzo
        if (N == 3) {
            float pred_freq, mean_freq;
            T mean_guess;
            lorenzo_config.quantbinCnt = optimize_quant_invl_3d<T>(
                data, conf.dims[0], conf.dims[1], conf.dims[2], conf.absErrorBound, pred_freq, mean_freq, mean_guess);
            lorenzo_config.pred_dim = 2;
            size_t sampleOutSize =
                SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), buffer, bufferCap);
            ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.pred_dim = 3;
            }
        }

        if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
            auto quant_num = lorenzo_config.quantbinCnt;
            lorenzo_config.quantbinCnt = 16384;
            size_t sampleOutSize =
                SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), buffer, bufferCap);
            //                delete[]cmprData;
            ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.quantbinCnt = quant_num;
            }
        }
        lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
        conf = lorenzo_config;
        //            double tuning_time = timer.stop();
        cmpSize = SZ_compress_LorenzoReg<T, N>(conf, data, cmpData, cmpCap);
    }

    free(buffer);
    return cmpSize;
}
}  // namespace SZ3
#endif
