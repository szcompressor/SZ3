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
#include "SZ3/utils/Sample.hpp"

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
    //conf.tuning = true;
    conf.interp_anchorStride = 0;
    
    auto sz = SZBlockInterpolationCompressor<T, N, LinearQuantizer<T>, HuffmanEncoder<int>, Lossless_zstd>(
        LinearQuantizer<T>(eb), HuffmanEncoder<int>(), Lossless_zstd());
    size_t outSize = sz.compress(conf, data1.data(), buffer, bufferCap);
    
    /*
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    
    size_t outSize = sz->compress(conf, data1.data(), buffer, bufferCap);
    */
    
    auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
    return compression_ratio;
}

template <class T, uint N>
double interp_compress_test_qoz(const std::vector< std::vector<T> > sampled_blocks, const Config conf,
                                                  int block_size, uchar *cmpData,
                                                  size_t cmpCap) {


    //conf.setDims(dims.begin(), dims.end());
    //conf.interpAlgo = interp_op;
    //conf.interpDirection = direction_op;
    //conf.tuning = true;

    
    
    auto sz = make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2));

    std::vector<int> total_quant_bins;
    /*
    int max_level = log2(block_size);
    std::vector<size_t> prefix(max_level,0);
    int side_length = block_size, level = 0;
    while(level < max_level){
        side_length = (block_size + 1) / 2;
        prefix[level++] = pow(side_length, N);
    }
    std::vector<std::vector<int> > block_q_bins;
    */
    for (int k = 0; k < sampled_blocks.size(); k++){
        auto cur_block = sampled_blocks[k];
    
        auto quant_bins = sz.compress(conf, cur_block.data());

        total_quant_bins.insert(total_quant_bins.end(), quant_bins.begin(), quant_bins.end());
        //block_q_bins.push_back(quant_bins);
    }
    /*
    size_t level_num = prefix.size();
    size_t last_pos = 0;
    for(int k = level_num - 1; k >= 0; k--){
        for (size_t l = 0; l < sampled_blocks.size(); l++){
            for (size_t m = last_pos; m < (k == 0 ? conf.num : prefix[k - 1]); m++){
                total_quant_bins.push_back(block_q_bins[l][m]);
            }
        }
        last_pos = prefix[k];
    }    
    std::cout<<total_quant_bins.size()<<std::endl;
    */


    auto encoder = HuffmanEncoder<int>();
    auto lossless = Lossless_zstd();
    encoder.preprocess_encode(total_quant_bins, sz.get_out_range().second);
    size_t bufferSize = std::max<size_t>(
        1000, 1.2 * (sz.size_est() + encoder.size_est() + sizeof(T) * total_quant_bins.size()));

    auto buffer = static_cast<uchar *>(malloc(bufferSize));
    uchar *buffer_pos = buffer;
    //sz.save(buffer_pos);
    //encoder.save(buffer_pos);

    //store the size of quant_inds is necessary as it is not always equal to conf.num
    write<size_t>(total_quant_bins.size(), buffer_pos);
    encoder.encode(total_quant_bins, buffer_pos);
    encoder.postprocess_encode();
    auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
    free(buffer);
    auto compression_ratio = conf.num * sampled_blocks.size() * sizeof(T) * 1.0 / cmpSize;
    return compression_ratio;
}

template <class T, uint N>
double lorenzo_compress_test_qoz(const std::vector< std::vector<T> > sampled_blocks, const Config &conf, uchar *cmpData, size_t cmpCap) {
    auto sz = make_decomposition_lorenzo_regression<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2));
    std::vector<int> total_quant_bins;

    for (int k = 0; k < sampled_blocks.size(); k++){
        auto cur_block = sampled_blocks[k];
    
        auto quant_bins = sz.compress(conf, cur_block.data());

        total_quant_bins.insert(total_quant_bins.end(), quant_bins.begin(), quant_bins.end());
    }

    auto encoder = HuffmanEncoder<int>();
    auto lossless = Lossless_zstd();
    encoder.preprocess_encode(total_quant_bins, sz.get_out_range().second);
    size_t bufferSize = std::max<size_t>(
        1000, 1.2 * (sz.size_est() + encoder.size_est() + sizeof(T) * total_quant_bins.size()));

    auto buffer = static_cast<uchar *>(malloc(bufferSize));
    uchar *buffer_pos = buffer;
    sz.save(buffer_pos);
    encoder.save(buffer_pos);

    //store the size of quant_inds is necessary as it is not always equal to conf.num
    write<size_t>(total_quant_bins.size(), buffer_pos);
    encoder.encode(total_quant_bins, buffer_pos);
    encoder.postprocess_encode();
    auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
    free(buffer);
    auto compression_ratio = conf.num * sampled_blocks.size() * sizeof(T) * 1.0 / cmpSize;
    return compression_ratio;


}


template <class T, uint N>
size_t SZ_compress_Interp_lorenzo(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(conf.cmprAlgo == ALGO_INTERP_LORENZO);

    //        Timer timer(true);
    calAbsErrorBound(conf, data);

    if (conf.interp_anchorStride <= 0){
        std::array<size_t, 4> anchor_strides = {4096, 128, 32, 16};
        conf.interp_anchorStride = anchor_strides[N - 1];
    }
    
    std::array<double,4> sample_Rates={0.01, 0.01, 0.005, 0.005};
    auto sampleRate = sample_Rates[N - 1];
    std::array<size_t,4> sampleBlock_Sizes={4096, 128, 32, 16};
    size_t sampleBlockSize = sampleBlock_Sizes[N - 1];
    size_t shortest_edge = conf.dims[0];
    for (size_t i = 0; i < N; i++){
        shortest_edge = conf.dims[i] < shortest_edge ? conf.dims[i] : shortest_edge;
    }

    while (sampleBlockSize >= shortest_edge)
        sampleBlockSize /= 2;

    while (sampleBlockSize >= 16 and (pow(sampleBlockSize + 1, N) / conf.num) > 1.5 * sampleRate)
        sampleBlockSize /= 2;

    if (sampleBlockSize < 8)
        sampleBlockSize = 8;

    std::vector< std::vector<T> > sampled_blocks;
    //size_t num_sampled_blocks;
    size_t per_block_ele_num = pow(sampleBlockSize + 1, N) ;
    size_t sampling_num;
    std::vector< std::vector<size_t> >starts;
    auto profStride = sampleBlockSize / 4;
    profiling_block<T, N>(data, conf.dims, starts, sampleBlockSize, conf.absErrorBound, profStride);
    size_t num_filtered_blocks = starts.size();
    bool profiling = num_filtered_blocks * per_block_ele_num >= 0.5 * sampleRate * conf.num;//temp. to refine
    //bool profiling = false;
    sampleBlocks<T, N>(data, conf.dims, sampleBlockSize, sampled_blocks, sampleRate, profiling, starts);
    sampling_num = sampled_blocks.size() * per_block_ele_num;


    //size_t sampling_num, sampling_block;
    //std::vector<size_t> sample_dims(N);
    //std::vector<T> sampling_data = sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
    if (sampling_num == 0 or sampling_num >= conf.num * 0.2) {
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
            std::vector<size_t> sample_dims(N, sampleBlockSize + 1);
            lorenzo_config.cmprAlgo = ALGO_LORENZO_REG;
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
            lorenzo_config.lorenzo = true;
            lorenzo_config.lorenzo2 = true;
            lorenzo_config.regression = false;
            lorenzo_config.regression2 = false;
            lorenzo_config.openmp = false;
            lorenzo_config.blockSize = 5;
            //        lorenzo_config.quantbinCnt = 65536 * 2;
            best_lorenzo_ratio = lorenzo_compress_test_qoz<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
            //            delete[]cmprData;
            //    printf("Lorenzo ratio = %.2f\n", ratio);
        }
    }
    {
        // tune interp
        conf.interpDirection = 0; 
        conf.interp_alpha = 1.25;
        conf.interp_beta = 2.0;
        auto testConfig = conf;
        std::vector<size_t> dims(N, sampleBlockSize + 1);
        testConfig.setDims(dims.begin(), dims.end());
        for (auto &interp_op : {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC, INTERP_ALGO_CUBIC_NATURAL}) {
            //ratio = interp_compress_test<T, N>(
            //    sampling_data.data(), conf, sample_dims, sampling_num, conf.absErrorBound, interp_op, conf.interpDirection,
            //    sampling_block, buffer, bufferCap);

            testConfig.interpDirection = interp_op;
            ratio = interp_compress_test_qoz<T, N>(
                sampled_blocks, testConfig, sampleBlockSize, buffer, bufferCap);
            std::cout<<static_cast<int>(interp_op)<<" "<<ratio<<std::endl;
            if (ratio > best_interp_ratio) {
                best_interp_ratio = ratio;
                conf.interpAlgo = interp_op;
            }
        }
        
        testConfig.interpAlgo = conf.interpAlgo;
        testConfig.interpDirection = factorial(N) - 1;
        //ratio = interp_compress_test<T, N>(sampling_data.data(), conf, sample_dims, sampling_num,
        //                                                         conf.absErrorBound, conf.interpAlgo, direction_op,
        //                                                         sampling_block, buffer, bufferCap);
        ratio = interp_compress_test_qoz<T, N>(sampled_blocks, testConfig, sampleBlockSize, buffer, bufferCap);
        std::cout<<"reverse "<<ratio<<std::endl;
        if (ratio > best_interp_ratio * 1.02) {
            best_interp_ratio = ratio;
            conf.interpDirection = testConfig.interpDirection;
        }
        testConfig.interpDirection = conf.interpDirection;
        auto alphalist = std::vector<double>{1.0, 1.5, 2.0};
        auto betalist = std::vector<double>{1.0, 2.5, 3.0};
        for (auto i = 0; i < alphalist.size(); i++){
            auto alpha = alphalist[i];
            auto beta = betalist[i];
            testConfig.interp_alpha = alpha;
            testConfig.interp_beta = beta;
            ratio = interp_compress_test_qoz<T, N>(sampled_blocks, testConfig, sampleBlockSize, buffer, bufferCap);
            std::cout<<alpha<<" "<<beta<<" "<<ratio<<std::endl;
            if (ratio > best_interp_ratio * 1.02) {
                best_interp_ratio = ratio;
                conf.interp_alpha = alpha;
                conf.interp_beta = beta;
            }

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
            ratio  =
                lorenzo_compress_test_qoz<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.pred_dim = 3;
            }
        }

        if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
            auto quant_num = lorenzo_config.quantbinCnt;
            lorenzo_config.quantbinCnt = 16384;
            ratio  =
                lorenzo_compress_test_qoz<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
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
