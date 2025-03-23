#ifndef SZ3_SZALGO_INTERP_QOZ_HPP
#define SZ3_SZALGO_INTERP_QOZ_HPP

#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
#include "SZ3/compressor/specialized/SZBlockInterpolationCompressor.hpp"
#include "SZ3/decomposition/InterpolationDecomposition_QoZ.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Sample.hpp"
#include "SZ3/compressor/SZLosslessCompressor.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Metrics.hpp"

namespace SZ3 {

namespace QoZ {
template <class T, uint N>
size_t SZ_compress_Interp(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_INTERP);
    calAbsErrorBound(conf, data);

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    return sz->compress(conf, data, cmpData, cmpCap);
    //        return cmpData;
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
double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num, double eb,
                                                  int interp_op, int direction_op, int block_size, uchar *buffer,
                                                  size_t bufferCap) {
    std::vector<T> data1(data, data + num);

    Config conf;
    conf.absErrorBound = eb;
    conf.setDims(dims.begin(), dims.end());
    conf.blockSize = block_size;
    conf.interpMeta.interpAlgo = interp_op;
    conf.interpMeta.interpDirection = direction_op;
    auto sz = SZBlockInterpolationCompressor<T, N, LinearQuantizer<T>, HuffmanEncoder<int>, Lossless_zstd>(
        LinearQuantizer<T>(eb), HuffmanEncoder<int>(), Lossless_zstd());

    size_t outSize = sz.compress(conf, data1.data(), buffer, bufferCap);

    auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
    return compression_ratio;
}

template<class T, uint N>
inline void init_alphalist(std::vector<double> &alpha_list,const double &rel_bound, Config &conf){

    
    /*
    if(conf.linearReduce){
        alpha_list={0,0.1,0.2,0.3,0.4,0.5};

    }
    */
    //else{
        if (conf.tuningTarget!=TUNING_TARGET_CR){
           
            if(conf.abList==0)            
                alpha_list={1,1.25,1.5,1.75,2};
            else if(conf.abList==1)
                alpha_list={1,1.25,1.5,1.75,2,2.25,2.5};
            else
                alpha_list={1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};
           
        }
        else{
            
            alpha_list={-1,1,1.25,1.5,1.75,2};
           
        }
    //}
}
template<class T, uint N>
inline void init_betalist(std::vector<double> &beta_list,const double &rel_bound, Config &conf){
   
    /*
    if(conf.linearReduce){
        beta_list={1,0.75,0.5,0.33,0.25};
    }
    */
    //else{
        if (conf.tuningTarget!=TUNING_TARGET_CR){    
            
            beta_list={1.5,2,3,4};//may remove 1.5
           
        }
        else {
          
            beta_list={-1,1.5,2,3};
           
        }
    //}
}




template<class T, uint N>
void sampleBlocks(T *data,std::vector<size_t> &dims, size_t sampleBlockSize,std::vector< std::vector<T> > & sampled_blocks,double sample_rate,int profiling ,std::vector<std::vector<size_t> > &starts,int var_first=0){
    for(int i=0;i<sampled_blocks.size();i++){
                std::vector< T >().swap(sampled_blocks[i]);                
            }
            std::vector< std::vector<T> >().swap(sampled_blocks);
    for(int i=0;i<sampled_blocks.size();i++){
        std::vector< T >().swap(sampled_blocks[i]);                  
    }
    std::vector< std::vector<T> >().swap(sampled_blocks);                               
    size_t totalblock_num=1;
    for(int i=0;i<N;i++){                        
        totalblock_num*=static_cast<int>((dims[i]-1)/sampleBlockSize);
    }               
    size_t idx=0,block_idx=0;   
    if(profiling){
        size_t num_filtered_blocks=starts.size();    
        if(var_first==0){  
            size_t sample_stride=static_cast<size_t>(num_filtered_blocks/(totalblock_num*sample_rate));
            if(sample_stride<=0)
                sample_stride=1;
            
            for(size_t i=0;i<num_filtered_blocks;i+=sample_stride){
                std::vector<T> s_block;
                sample_blocks<T,N>(data, s_block,dims, starts[i],sampleBlockSize+1);
                sampled_blocks.push_back(s_block);
                
            }
            
        }
        else{
            std::vector< std::pair<double,std::vector<size_t> > >block_heap;
            for(size_t i=0;i<num_filtered_blocks;i++){
                double mean,sigma2,range;
                blockwise_profiling<T>(data,dims, starts[i],sampleBlockSize+1, mean,sigma2,range);
                block_heap.push_back(std::pair<double,std::vector<size_t> >(sigma2,starts[i]));
                
            }
            std::make_heap(block_heap.begin(),block_heap.end());
          

            size_t sampled_block_num=totalblock_num*sample_rate;
            if(sampled_block_num>num_filtered_blocks)
                sampled_block_num=num_filtered_blocks;
            if(sampled_block_num==0)
                sampled_block_num=1;

            for(size_t i=0;i<sampled_block_num;i++){
                std::vector<T> s_block;
             
                sample_blocks<T,N>(data, s_block,dims, block_heap.front().second,sampleBlockSize+1);
              
                sampled_blocks.push_back(s_block);
                std::pop_heap(block_heap.begin(),block_heap.end());
                block_heap.pop_back();
               
            }
        }
    }               
    else{
        if(var_first==0){
            size_t sample_stride=static_cast<size_t>(1.0/sample_rate);
            if(sample_stride<=0)
                sample_stride=1;
            if (N==2){                        
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                           
                    for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                        if (idx%sample_stride==0){
                            std::vector<size_t> starts{x_start,y_start};
                            std::vector<T> s_block;
                            sample_blocks<T,N>(data, s_block,dims, starts,sampleBlockSize+1);
                            sampled_blocks.push_back(s_block);
                        }
                        idx+=1;
                    }
                }
            }
            else if (N==3){                  
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                          
                    for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                        for (size_t z_start=0;z_start<dims[2]-sampleBlockSize;z_start+=sampleBlockSize){
                            if (idx%sample_stride==0){
                                std::vector<size_t> starts{x_start,y_start,z_start};
                                std::vector<T> s_block;
                                sample_blocks<T,N>(data, s_block,dims, starts,sampleBlockSize+1);
                                sampled_blocks.push_back(s_block);
                            }
                            idx+=1;
                        }
                    }
                }
            }
        }
        else{
            std::vector <std::vector<size_t> > blocks_starts;
            if (N==2){  
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                           
                    for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                       
                            blocks_starts.push_back(std::vector<size_t>{x_start,y_start});
                    }
                }

            }
            else if (N==3){           
                for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                          
                    for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                        for (size_t z_start=0;z_start<dims[2]-sampleBlockSize;z_start+=sampleBlockSize){
                            blocks_starts.push_back(std::vector<size_t>{x_start,y_start,z_start});
                        }
                    }
                }
            

                std::vector< std::pair<double,std::vector<size_t> > >block_heap;
                for(size_t i=0;i<totalblock_num;i++){
                    double mean,sigma2,range;
                    blockwise_profiling<T>(data,dims, blocks_starts[i],sampleBlockSize+1, mean,sigma2,range);
                    block_heap.push_back(std::pair<double,std::vector<size_t> >(sigma2,blocks_starts[i]));
                }
                std::make_heap(block_heap.begin(),block_heap.end());
                size_t sampled_block_num=totalblock_num*sample_rate;
                if(sampled_block_num==0)
                    sampled_block_num=1;
                for(size_t i=0;i<sampled_block_num;i++){
                    std::vector<T> s_block;
                    sample_blocks<T,N>(data, s_block,dims, block_heap.front().second,sampleBlockSize+1);
                    sampled_blocks.push_back(s_block);
                    std::pop_heap(block_heap.begin(),block_heap.end());
                    block_heap.pop_back();
                }

            }
        }
    }
}


template<class T, uint N>
std::pair<double,double> CompressTest(const Config &conf,const std::vector< std::vector<T> > & sampled_blocks,ALGO algo = ALGO_INTERP,
                    TUNING_TARGET tuningTarget=TUNING_TARGET_RD,bool useFast=true,double profiling_coeff=1,const std::vector<double> &orig_means=std::vector<double>(),
                    const std::vector<double> &orig_sigma2s=std::vector<double>(),const std::vector<double> &orig_ranges=std::vector<double>(),const std::vector<T> &flattened_sampled_data=std::vector<T>()){
    Config testConfig(conf);
    size_t ssim_size=conf.SSIMBlockSize;    
    if(algo == ALGO_LORENZO_REG){
        testConfig.cmprAlgo = ALGO_LORENZO_REG;
        testConfig.dims=conf.dims;
        testConfig.num=conf.num;
        testConfig.lorenzo = true;
        testConfig.lorenzo2 = true;
        testConfig.regression = false;
        testConfig.regression2 = false;
        testConfig.openmp = false;
        testConfig.blockSize = 5;//why?
        testConfig.quantbinCnt = 65536 * 2;
    }
    double square_error=0.0;
    double bitrate=0.0;
    double metric=0.0;
    size_t sampleBlockSize=testConfig.sampleBlockSize;
    size_t num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(sampleBlockSize+1,N);
    size_t ele_num=num_sampled_blocks*per_block_ele_num;
    std::vector<T> cur_block(testConfig.num,0);
    std::vector<int> q_bins;
    std::vector<std::vector<int> > block_q_bins;
    std::vector<size_t> q_bin_counts;
    std::vector<T> flattened_cur_blocks;
    size_t idx=0;   
    std::shared_ptr< concepts::DecompositionInterface<T,int,N>> sz;
    //size_t totalOutSize=0;
    if(algo == ALGO_LORENZO_REG){
        auto sz = make_decomposition_lorenzo_regression<T, N, LinearQuantizer<T> >(
                        testConfig,
                        LinearQuantizer<T>(testConfig.absErrorBound, testConfig.quantbinCnt / 2)
                        );
        for (int k=0;k<num_sampled_blocks;k++){
            size_t sampleOutSize;
            std::vector<T> cur_block(testConfig.num);
            std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());
            
             
            double decomp_square_error;
            auto quant_bins = sz->compress(testConfig, cur_block.data());

            
          
            q_bins.insert(q_bins.end(),quant_bins.start(),quant_bins.end());

            if(tuningTarget==TUNING_TARGET_RD){
                if(algo==ALGO_INTERP)
                    square_error+=decomp_square_error;
                else{
                   
                    for(size_t j=0;j<per_block_ele_num;j++){
                        T value=sampled_blocks[k][j]-cur_block[j];
                        square_error+=value*value;
                    
                    }
                }
            }
            else if (tuningTarget==TUNING_TARGET_SSIM){
                size_t ssim_block_num=orig_means.size();                       
                double mean=0,sigma2=0,cov=0,range=0;
                double orig_mean=0,orig_sigma2=0,orig_range=0;  
                std::vector<size_t>block_dims(N,sampleBlockSize+1);                      
                if(N==2){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            orig_mean=orig_means[idx];
                            orig_sigma2=orig_sigma2s[idx];
                            orig_range=orig_ranges[idx];
                            std::vector<size_t> starts{i,j};
                            blockwise_profiling<T>(cur_block.data(),block_dims,starts,ssim_size,mean,sigma2,range);
                            cov=blockwise_cov<T>(sampled_blocks[k].data(),cur_block.data(),block_dims,starts,ssim_size,orig_mean,mean);
                            metric+=SSIM(orig_range,orig_mean,orig_sigma2,mean,sigma2,cov)/ssim_block_num;
                            idx++;


                        }
                    }
                }
                else if(N==3){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            for (size_t kk=0;kk+ssim_size<sampleBlockSize+1;kk+=ssim_size){
                                orig_mean=orig_means[idx];
                                orig_sigma2=orig_sigma2s[idx];
                                orig_range=orig_ranges[idx];
                                std::vector<size_t> starts{i,j,kk};
                                blockwise_profiling<T>(cur_block.data(),block_dims,starts,ssim_size,mean,sigma2,range);
                                cov=blockwise_cov<T>(sampled_blocks[k].data(),cur_block.data(),block_dims,starts,ssim_size,orig_mean,mean);
                                //printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",orig_range,orig_sigma2,orig_mean,range,sigma2,mean,cov);
                                metric+=SSIM(orig_range,orig_mean,orig_sigma2,mean,sigma2,cov)/ssim_block_num;
                         
                                idx++;
                            }
                        }
                    }
                }
            }
            else if (tuningTarget==TUNING_TARGET_AC){
                flattened_cur_blocks.insert(flattened_cur_blocks.end(),cur_block.begin(),cur_block.end());
            }                      
        }


    }
    else if(algo == ALGO_INTERP){

        auto sz = make_decomposition_interpolation<T, N, LinearQuantizer<T> >(
                        testConfig,
                        LinearQuantizer<T>(testConfig.absErrorBound, testConfig.quantbinCnt / 2)
                        );


        
        for (int k=0;k<num_sampled_blocks;k++){
            size_t sampleOutSize;
            std::vector<T> cur_block(testConfig.num);
            std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());
            
             
            double decomp_square_error;
            auto quant_bins = sz->compress(testConfig, cur_block.data(), 1,decomp_square_error,q_bin_counts);

            
            

            
            block_q_bins.push_back(quant_bins);
            
        

            if(tuningTarget==TUNING_TARGET_RD){
                if(algo==ALGO_INTERP)
                    square_error+=decomp_square_error;
                else{
                   
                    for(size_t j=0;j<per_block_ele_num;j++){
                        T value=sampled_blocks[k][j]-cur_block[j];
                        square_error+=value*value;
                    
                    }
                }
            }
            else if (tuningTarget==TUNING_TARGET_SSIM){
                size_t ssim_block_num=orig_means.size();                       
                double mean=0,sigma2=0,cov=0,range=0;
                double orig_mean=0,orig_sigma2=0,orig_range=0;  
                std::vector<size_t>block_dims(N,sampleBlockSize+1);                      
                if(N==2){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            orig_mean=orig_means[idx];
                            orig_sigma2=orig_sigma2s[idx];
                            orig_range=orig_ranges[idx];
                            std::vector<size_t> starts{i,j};
                            blockwise_profiling<T>(cur_block.data(),block_dims,starts,ssim_size,mean,sigma2,range);
                            cov=blockwise_cov<T>(sampled_blocks[k].data(),cur_block.data(),block_dims,starts,ssim_size,orig_mean,mean);
                            metric+=SSIM(orig_range,orig_mean,orig_sigma2,mean,sigma2,cov)/ssim_block_num;
                            idx++;


                        }
                    }
                }
                else if(N==3){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            for (size_t kk=0;kk+ssim_size<sampleBlockSize+1;kk+=ssim_size){
                                orig_mean=orig_means[idx];
                                orig_sigma2=orig_sigma2s[idx];
                                orig_range=orig_ranges[idx];
                                std::vector<size_t> starts{i,j,kk};
                                blockwise_profiling<T>(cur_block.data(),block_dims,starts,ssim_size,mean,sigma2,range);
                                cov=blockwise_cov<T>(sampled_blocks[k].data(),cur_block.data(),block_dims,starts,ssim_size,orig_mean,mean);
                                //printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",orig_range,orig_sigma2,orig_mean,range,sigma2,mean,cov);
                                metric+=SSIM(orig_range,orig_mean,orig_sigma2,mean,sigma2,cov)/ssim_block_num;
                         
                                idx++;
                            }
                        }
                    }
                }
            }
            else if (tuningTarget==TUNING_TARGET_AC){
                flattened_cur_blocks.insert(flattened_cur_blocks.end(),cur_block.begin(),cur_block.end());
            }                      
        }
        

    }
    else{
        if(conf.verbose)
            std::cout<<"algo type error!"<<std::endl;
        return std::pair<double,double>(0,0);
    }

    if(algo==ALGO_INTERP ){
        size_t level_num=q_bin_counts.size();
        size_t last_pos=0;
        for(int k=level_num-1;k>=0;k--){
            for (size_t l =0;l<num_sampled_blocks;l++){
                for (size_t m=last_pos;m<q_bin_counts[k];m++){
                    q_bins.push_back(block_q_bins[l][m]);
                }
            }
            last_pos=q_bin_counts[k];
        }      
    }
    size_t sampleOutSize;
    
    //auto cmprData=sz->encoding_lossless(totalOutSize,q_bins);             
    //delete[]cmprData;
    size_t estimated_size = 1.2*sizeof(int)*q_bins.size();

    uchar * cmpData = new uchar[estimated_size];
    auto lossless = make_compressor_sz_encodinglossless<HuffmanEncoder<int>, Lossless_zstd>(HuffmanEncoder<int>(), Lossless_zstd());

    auto totalOutSize = lossless.compress(q_bins,cmpData,estimated_size);

    delete [] cmpData;


  
    
    bitrate=8*double(totalOutSize)/ele_num;
    
    bitrate*=profiling_coeff;
    if(tuningTarget==TUNING_TARGET_RD){
        double mse=square_error/ele_num;
        mse*=profiling_coeff;      
        
        metric=PSNR(testConfig.rng,mse);
    }
    else if (tuningTarget==TUNING_TARGET_AC){                       
        metric=1.0-autocorrelation<T>(flattened_sampled_data.data(),flattened_cur_blocks.data(),ele_num);                        
    }                    
   

    if(algo==ALGO_LORENZO_REG)    {
        bitrate*=testConfig.lorenzoBrFix;
    }
    //delete sz;
    return std::pair<double,double>(bitrate,metric);
}

std::pair <double,double> setABwithRelBound(double rel_bound,int configuration=0){

    double cur_alpha=-1,cur_beta=-1;
    if(configuration==0){              
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=2;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=2;
        }                 
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }                   
        else if (rel_bound>0.001){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
            cur_alpha=1;
            cur_beta=1;
        }
    }
    else if(configuration==1){                
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=4;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=3;
        }
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }           
        else if (rel_bound>0.001){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
                cur_alpha=1;
                cur_beta=1;
            }
    }
    else if(configuration==2){                
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=4;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=3;
        }                    
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }
        else if (rel_bound>0.001){
            cur_alpha=1.5;
            cur_beta=1.5;
        }
        else if (rel_bound>0.0005){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
            cur_alpha=1;
            cur_beta=1;
        }
    }
    return std::pair<double,double>(cur_alpha,cur_beta);
}

void setLorenzoFixRates(Config &conf,double rel_bound){
    double e1=1e-5;
    double e2=1e-4;
    double e3=1e-3;
    //double e4=1e-1;
    /*
    double f1=conf.sampleBlockSize>=64?2: 1;
    // double f2=1.1;old
    double f2=conf.sampleBlockSize>=64?3:1.2; 
    // double f3=1.2;//old need to raise 
    //double f3=1.3;
    double f3=conf.sampleBlockSize>=64?4:1.3;
    */
    double f1=conf.sampleBlockSize>=64?2: 1.1;//modified from 1.0. need further test
    // double f2=1.1;old
    double f2=conf.sampleBlockSize>=64?3:1.2; 
    // double f3=1.2;//old need to raise 
    //double f3=1.3;
    double f3=conf.sampleBlockSize>=64?4:1.3;
    if(rel_bound<=e1)
        conf.lorenzoBrFix=f1;
    else if(rel_bound<=e2)
        conf.lorenzoBrFix=f1-(f1-f2)*(rel_bound-e1)/(e2-e1);
    else if (rel_bound<=e3)
        conf.lorenzoBrFix=f2-(f2-f3)*(rel_bound-e2)/(e3-e2);
    else 
        conf.lorenzoBrFix=f3;
}

template<class T, uint N>
double Tuning(Config &conf, T *data){
   
    
    
   // Timer timer(true);
    //timer.stop("")
    if(conf.QoZ>0){
        
        //testLorenzo?
        

        //activate
        //conf.testLorenzo=0;//temporarily deactivated Lorenzo. Need further revision
        conf.profiling=1;
        if(conf.autoTuningRate<=0)
            conf.autoTuningRate = (N<=2?0.01:0.005);
        if(conf.predictorTuningRate<=0)
            conf.predictorTuningRate = (N<=2?0.01:0.005);
        if (conf.maxStep<=0){
            std::array<size_t,4> anchor_strides={256,64,32,16};
            conf.maxStep = anchor_strides[N-1];
        }
        if (conf.levelwisePredictionSelection<=0)
            conf.levelwisePredictionSelection = (N<=2?5:4);
        if (conf.sampleBlockSize<=0){
            
            conf.sampleBlockSize = (N<=2?64:32);
        }

        if(conf.QoZ>=2){
            //conf.testLorenzo=1;
            conf.multiDimInterp=1;
            conf.naturalSpline=1;
            conf.fullAdjacentInterp=1;
            conf.freezeDimTest=1;
            
        }
        if(conf.QoZ>=3){
            conf.dynamicDimCoeff=1;
        }
        if(conf.QoZ>=4){
            conf.crossBlock=1;
            conf.blockwiseTuning=1;
            if(conf.blockwiseSampleRate<1.0)
                conf.blockwiseSampleRate=3.0;
        }
    }   
    //Deactivate several features when dim not fit
    if(N!=2 and N!=3){
       // conf.QoZ=0; comment it for level-wise eb
        conf.autoTuningRate=0;
        conf.predictorTuningRate=0;
        conf.levelwisePredictionSelection=0;
        conf.multiDimInterp=0;
        conf.naturalSpline=0;
        conf.fullAdjacentInterp=0;
        conf.freezeDimTest=0;
        conf.dynamicDimCoeff=0;
        conf.blockwiseTuning=0;
    }
    

    if(conf.multiDimInterp==0)
        conf.dynamicDimCoeff=0;
    size_t sampling_num, sampling_block;
    double best_interp_cr=0.0;
    double best_lorenzo_ratio=0.0;
    bool useInterp=true;        
    std::vector<size_t> sample_dims(N);
    std::vector<T> sampling_data;
    double anchor_rate=0;
    int max_interp_level = -1;

    for (size_t i = 0; i < N; i++) {
        if ( max_interp_level < ceil(log2(conf.dims[i]))) {
             max_interp_level = static_cast<uint> (ceil(log2(conf.dims[i])));
        }
                
    }
    
    if (conf.maxStep>0){
        anchor_rate=1/(pow(conf.maxStep,N));   
        int temp_max_interp_level=static_cast<uint>(log2(conf.maxStep));//to be catious: the max_interp_level is different from the ones in szinterpcompressor, which includes the level of anchor grid.
        if (temp_max_interp_level<=max_interp_level){                  
            max_interp_level=temp_max_interp_level;
        }
        if (conf.levelwisePredictionSelection>max_interp_level)
            conf.levelwisePredictionSelection=max_interp_level;
    }
            
    

    std::vector<Interp_Meta> bestInterpMeta_list(conf.levelwisePredictionSelection);
    Interp_Meta bestInterpMeta;


    size_t shortest_edge=conf.dims[0];
    for (size_t i=0;i<N;i++){
        shortest_edge=conf.dims[i]<shortest_edge?conf.dims[i]:shortest_edge;
    }
 

    
    if (conf.sampleBlockSize<=0){
        conf.sampleBlockSize = (N==2?64:32);
            
    }

    
    
    size_t minimum_sbs=8;
    if (conf.sampleBlockSize<minimum_sbs)
        conf.sampleBlockSize=minimum_sbs;


    while(conf.sampleBlockSize>=shortest_edge)
        conf.sampleBlockSize/=2;


    

    while(conf.autoTuningRate>0 and conf.sampleBlockSize>=2*minimum_sbs and (pow(conf.sampleBlockSize+1,N)/conf.num)>1.5*conf.autoTuningRate)
        conf.sampleBlockSize/=2;

    if (conf.sampleBlockSize<minimum_sbs){
        conf.predictorTuningRate=0.0;
        conf.autoTuningRate=0.0;
    }
    else{
        int max_lps_level=static_cast<uint>(log2(conf.sampleBlockSize));//to be catious: the max_interp_level is different from the ones in szinterpcompressor, which includes the level of anchor grid.

        if (conf.levelwisePredictionSelection>max_lps_level)
            conf.levelwisePredictionSelection=max_lps_level;
    }

    std::vector< std::vector<T> > sampled_blocks;
    size_t sampleBlockSize=conf.sampleBlockSize;
    size_t num_sampled_blocks;
    size_t per_block_ele_num;
    size_t ele_num;

    
           
    size_t totalblock_num=1;  
    for(int i=0;i<N;i++){                      
        totalblock_num*=static_cast<size_t>((conf.dims[i]-1)/sampleBlockSize);
    }

    std::vector<std::vector<size_t> >starts;
    if((conf.autoTuningRate>0 or conf.predictorTuningRate>0) and conf.profiling){      
        conf.profStride=conf.sampleBlockSize/4;
        if(N==2){
            profiling_block_2d<T,N>(data,conf.dims,starts,sampleBlockSize,conf.absErrorBound,conf.profStride);
        }
        else if (N==3){
            profiling_block_3d<T,N>(data,conf.dims,starts,sampleBlockSize,conf.absErrorBound,conf.profStride);
        }
       
    }


    size_t num_filtered_blocks=starts.size();
    if(num_filtered_blocks<=static_cast<int>(0.3*conf.predictorTuningRate))//temp. to refine
        conf.profiling=0;
    double profiling_coeff=1;//It seems that this coefficent is useless. Need further test
  
    if(conf.profiling){
        profiling_coeff=static_cast<double>(num_filtered_blocks)/totalblock_num;
    }
    std::vector<size_t> global_dims=conf.dims;
    size_t global_num=conf.num;


    
    double rel_bound = 0.0;
    
    if (conf.autoTuningRate>0){

        if (conf.rng<0)
            conf.rng=data_range<T>(data,conf.num);
        if(conf.relErrorBound<=0)
            conf.relErrorBound=conf.absErrorBound/conf.rng;
        rel_bound = conf.relErrorBound;
        if(rel_bound>=3e-4 or conf.tuningTarget==TUNING_TARGET_SSIM)//rencently changed, need to fix later
            conf.testLorenzo=0;
        if (conf.testLorenzo>0)
            setLorenzoFixRates(conf,rel_bound);
    }


    bool blockwiseTuning=conf.blockwiseTuning;
    conf.blockwiseTuning=false;


    if (conf.predictorTuningRate>0 and conf.predictorTuningRate<1){
        if (conf.verbose)
            std::cout<<"Predictor tuning started."<<std::endl;
        double o_alpha=conf.alpha;
        double o_beta=conf.beta;
                    
        
        sampleBlocks<T,N>(data,conf.dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
              
        num_sampled_blocks=sampled_blocks.size();
        per_block_ele_num=pow(sampleBlockSize+1,N);
        ele_num=num_sampled_blocks*per_block_ele_num;
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;
        std::vector<T> cur_block(per_block_ele_num,0);
        
        
        
            
        double ori_eb=conf.absErrorBound;
        std::vector<size_t> coeffs_size;
        
        if(conf.testLorenzo and conf.autoTuningRate==0 ){

            std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,ALGO_LORENZO_REG,TUNING_TARGET_CR,false);
            best_lorenzo_ratio=sizeof(T)*8.0/results.first;
            
            if(conf.verbose)
                std::cout << "lorenzo best cr = " << best_lorenzo_ratio << std::endl;

        }

        if (conf.autoTuningRate>0){

           
            std::pair<double,double> ab=setABwithRelBound(rel_bound,0);//ori pdtuningqabconf
            conf.alpha=ab.first;
            conf.beta=ab.second;
       
        }
        std::vector<uint8_t> interpAlgo_Candidates={INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC};
        //if(conf.quadInterp){//deprecated
         //   interpAlgo_Candidates.push_back(INTERP_ALGO_QUAD);
        //}

        
        std::vector<uint8_t> interpParadigm_Candidates={0};//
        std::vector<uint8_t> cubicSplineType_Candidates={0};
        std::vector<uint8_t> interpDirection_Candidates={0, (uint8_t)(factorial(N) -1)};
       
        std::vector<uint8_t> adjInterp_Candidates={0};//

        
        if(conf.multiDimInterp>0){
            for(size_t i=1;i<=conf.multiDimInterp;i++)
                interpParadigm_Candidates.push_back(i);
        }

        if (conf.naturalSpline){
            cubicSplineType_Candidates.push_back(1);
        }


    
        
        if(conf.fullAdjacentInterp){
            adjInterp_Candidates.push_back(1);
        }
        
        if(conf.levelwisePredictionSelection>0 and (N==2 or N==3)){
            std::vector<Interp_Meta> interpMeta_list(conf.levelwisePredictionSelection);
            auto sz = make_decomposition_interpolation(conf, LinearQuantizer<T>(conf.absErrorBound,conf.quantbinCnt / 2));  
            double best_accumulated_interp_loss_1=0;
            double best_accumulated_interp_loss_2=0;
            std::vector<std::vector<double> > linear_interp_vars(conf.levelwisePredictionSelection),cubic_noknot_vars(conf.levelwisePredictionSelection),cubic_nat_vars(conf.levelwisePredictionSelection);
            
            for(int level=conf.levelwisePredictionSelection;level>0;level--){
              
                int start_level=(level==conf.levelwisePredictionSelection?9999:level);
                int end_level=level-1;
                
                if((conf.multiDimInterp>0 and conf.dynamicDimCoeff>0) or (conf.freezeDimTest>0 and level==1 and N>=3) ){
                    
                    
                    size_t interp_stride=pow(2,level-1);
                    size_t stride;
                    double cur_eb=level>=3?conf.absErrorBound/2:conf.absErrorBound;
                    if(level>=3)
                        stride=2;
                    else if(level>=2)
                        stride= conf.adaptiveMultiDimStride<=4?2:conf.adaptiveMultiDimStride/2;
                    else
                        stride= conf.adaptiveMultiDimStride;
                    if(conf.multiDimInterp>0 and conf.dynamicDimCoeff){
                        calculate_interp_error_vars<T,N>(data, global_dims,linear_interp_vars[level-1],0,0,stride,interp_stride,cur_eb);
                        preprocess_vars<N>(linear_interp_vars[level-1]);
                       
                        if (conf.naturalSpline){
                            calculate_interp_error_vars<T,N>(data, global_dims,cubic_nat_vars[level-1],1,1,stride,interp_stride,cur_eb);
                            preprocess_vars<N>(cubic_nat_vars[level-1]);
                         
                        }
                    }
                    calculate_interp_error_vars<T,N>(data, global_dims,cubic_noknot_vars[level-1],1,0,stride,interp_stride,cur_eb);
                    preprocess_vars<N>(cubic_noknot_vars[level-1]);
                   
                }


                double best_interp_absloss=std::numeric_limits<double>::max();
                //conf.cmprAlgo = ALGO_INTERP;    
                Interp_Meta cur_meta;
                Interp_Meta best_meta;

                for (auto &interp_op: interpAlgo_Candidates) {
                    cur_meta.interpAlgo=interp_op;
                    for (auto &interp_pd: interpParadigm_Candidates) {
                        cur_meta.interpParadigm=interp_pd;

                        

                        for (auto &interp_direction: interpDirection_Candidates) {
                            if ((interp_pd==1 or  (interp_pd==2 and N<=2)) and interp_direction!=0)
                                continue;
                            cur_meta.interpDirection=interp_direction;
                            for(auto &cubic_spline_type:cubicSplineType_Candidates){
                                if (interp_op!=INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                    break;
                                cur_meta.cubicSplineType=cubic_spline_type;
                                for(auto adj_interp:adjInterp_Candidates){
                                    if (interp_op!=INTERP_ALGO_CUBIC and adj_interp!=0)
                                        break;
                                    
                                    cur_meta.adjInterp=adj_interp;
                                    
                                    if(conf.dynamicDimCoeff>0 and interp_pd>0){
                                        if(interp_op==0){
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=linear_interp_vars[level-1][i];
                                        }
                                        else if (cubic_spline_type==0){
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=cubic_noknot_vars[level-1][i];
                                        }
                                        else{
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=cubic_nat_vars[level-1][i];
                                        }

                                    }
                                    
                                    conf.interpMeta=cur_meta;


                                    double cur_absloss=0;
                                    for (int i=0;i<num_sampled_blocks;i++){
                                        cur_block=sampled_blocks[i];  //not so efficient              
                                        double decomp_square_error;    
                                        std::vector<size_t> quant_bin_counts;                  
                                        sz.compress(conf, cur_block.data(), 2,start_level,end_level,decomp_square_error,quant_bin_counts);
                                        //delete []cmprData;                              
                                        cur_absloss+=decomp_square_error;
                                    }
                                    
                                    if (cur_absloss<best_interp_absloss){
                                        best_meta=cur_meta;
                                        best_interp_absloss=cur_absloss;
                                    }
                                    cur_meta.dimCoeffs={1.0/3.0,1.0/3.0,1.0/3.0};
                    
                                }
                            }
                        }   
                    }
                }
                best_accumulated_interp_loss_1+=best_interp_absloss;
      
                interpMeta_list[level-1]=best_meta;
                /*
                if(conf.pdTuningRealComp){
                    //place to add real compression,need to deal the problem that the sampled_blocks are changed. 
                    
                    conf.interpMeta=best_meta;
                    for (int i=0;i<num_sampled_blocks;i++){

                        size_t outSize=0;
                        double decomp_square_error;  
                        sz.compress(conf, sampled_blocks[i].data(),2,start_level,end_level,decomp_square_error);
                        //delete []cmprData;
                    }
                    
                } */
                

            }
            
            bestInterpMeta_list=interpMeta_list;

           
            /*
            
            if(conf.pdTuningRealComp and ((conf.autoTuningRate>0 and conf.autoTuningRate==conf.predictorTuningRate) or (conf.adaptiveMultiDimStride>0 and N==3))){
                    //recover sample if real compression used                  
                sampleBlocks<T,N>(data,global_dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
            }
            */
                

            //frozendim
            
            if(conf.freezeDimTest and N==3 ){

                std::vector<Interp_Meta> tempmeta_list=conf.interpMeta_list;
                conf.interpMeta_list=interpMeta_list;      
                std::pair<double,double> results=CompressTest<T,N>(conf,sampled_blocks,ALGO_INTERP,TUNING_TARGET_CR,false);
                double best_interp_cr_1=sizeof(T)*8.0/results.first;     
                conf.interpMeta_list=tempmeta_list;
                




                size_t frozen_dim=0;
                double cur_weight=cubic_noknot_vars[0][0];
                for(size_t i=1;i<N;i++){
                    if(cubic_noknot_vars[0][i]<cur_weight){
                        frozen_dim=i;
                        cur_weight=cubic_noknot_vars[0][i];
                    }
                }
                if(frozen_dim==0)
                    interpDirection_Candidates={6,7};
                else if (frozen_dim==1)
                    interpDirection_Candidates={8,9};
                else
                    interpDirection_Candidates={10,11};


                for(int level=conf.levelwisePredictionSelection;level>0;level--){
                    int start_level=(level==conf.levelwisePredictionSelection?9999:level);
                    int end_level=level-1;
                  
                    double best_interp_absloss=std::numeric_limits<double>::max();
        
                   Interp_Meta cur_meta;
                   Interp_Meta best_meta;

                    for (auto &interp_op: interpAlgo_Candidates) {
                        cur_meta.interpAlgo=interp_op;
                        for (auto &interp_pd: interpParadigm_Candidates) {
                            if (interp_pd>1)
                                continue;
                            cur_meta.interpParadigm=interp_pd;

                            

                            for (auto &interp_direction: interpDirection_Candidates) {
                                if (interp_pd>=1  and interp_direction%2!=0)//only dims[0] matters.
                                    continue;
                                cur_meta.interpDirection=interp_direction;
                                for(auto &cubic_spline_type:cubicSplineType_Candidates){
                                    if (interp_op!=INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                        break;
                                    cur_meta.cubicSplineType=cubic_spline_type;
                                    for(auto adj_interp:adjInterp_Candidates){
                                        if (interp_op!=INTERP_ALGO_CUBIC and adj_interp!=0)
                                            break;
                                       
                                        cur_meta.adjInterp=adj_interp;

                                        if(conf.dynamicDimCoeff>0 and interp_pd>0){
                                            if(interp_op==0){
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=linear_interp_vars[level-1][i];
                                            }
                                            else if (cubic_spline_type==0){
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=cubic_noknot_vars[level-1][i];
                                            }
                                            else{
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=cubic_nat_vars[level-1][i];
                                            }

                                        }
                                        
                                        conf.interpMeta=cur_meta;

                                        
                                        double cur_absloss=0;
                                        for (int i=0;i<num_sampled_blocks;i++){
                                            cur_block=sampled_blocks[i];  //not so efficient              
                                           double decomp_square_error;           
                                           std::vector<size_t> quant_bin_counts;            
                                           sz.compress(conf, cur_block.data(), 2,start_level,end_level,decomp_square_error,quant_bin_counts);                          
                                            cur_absloss+=decomp_square_error;
                                        }
                                        if (cur_absloss<best_interp_absloss){
                                            best_meta=cur_meta;
                                            best_interp_absloss=cur_absloss;
                                        }
                                        cur_meta.dimCoeffs={1.0/3.0,1.0/3.0,1.0/3.0};
                                        
                                    }
                                }
                            }   
                        }
                    }
                    best_accumulated_interp_loss_2+=best_interp_absloss;
          
                 

                    interpMeta_list[level-1]=best_meta;
                        /*
                    if(conf.pdTuningRealComp){
                        //place to add real compression,need to deal the problem that the sampled_blocks are changed. 
                  
                        conf.interpMeta=best_meta;
                        for (int i=0;i<num_sampled_blocks;i++){

                            size_t outSize=0;
                                       
                            auto cmprData =sz.compress(conf, sampled_blocks[i].data(), outSize,2,start_level,end_level);
                            delete []cmprData;
                        }
                        
                    } */
                }
                

                tempmeta_list=conf.interpMeta_list;
                conf.interpMeta_list=interpMeta_list;      
                results=CompressTest<T,N>(conf,sampled_blocks,ALGO_INTERP,TUNING_TARGET_CR,false);
                double best_interp_cr_2=sizeof(T)*8.0/results.first;     
                conf.interpMeta_list=tempmeta_list;

                if(best_interp_cr_2>best_interp_cr_1*1.05 and conf.verbose){
                    conf.frozen_dim=frozen_dim;
                    bestInterpMeta_list=interpMeta_list;
                    if(conf.verbose)
                        std::cout<<"Dim "<<frozen_dim<<" frozen"<<std::endl;
                }
            

                if(conf.pdTuningRealComp and conf.autoTuningRate>0 and conf.autoTuningRate==conf.predictorTuningRate){
                        //recover sample if real compression used                  
                    sampleBlocks<T,N>(data,global_dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
                }



            }
            
            



            
            conf.interpMeta_list=bestInterpMeta_list;
            if(conf.autoTuningRate==0){ //Qustionable part.  //when adaptivemdtride>0 there's a duplication of work. To fix.
              
                      
                std::pair<double,double> results=CompressTest<T,N>(conf,sampled_blocks,ALGO_INTERP,TUNING_TARGET_CR,false);
                double cur_best_interp_cr=sizeof(T)*8.0/results.first;     
               
                best_interp_cr=cur_best_interp_cr;
                   
                   
            }
           
        }

        else{
            
            Interp_Meta best_meta,cur_meta;
            double cur_best_interp_cr=0.0;
            for (auto &interp_op: interpAlgo_Candidates) {
                cur_meta.interpAlgo=interp_op;
                for (auto &interp_pd: interpParadigm_Candidates) {
                    cur_meta.interpParadigm=interp_pd;
                    for (auto &interp_direction: interpDirection_Candidates) {
                        if (interp_pd==1 or  (interp_pd==2 and N<=2) and interp_direction!=0)
                            continue;
                        cur_meta.interpDirection=interp_direction;
                        for(auto &cubic_spline_type:cubicSplineType_Candidates){
                            if (interp_op!=INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                break;
                            cur_meta.cubicSplineType=cubic_spline_type;
                            for(auto adj_interp:adjInterp_Candidates){
                                if (interp_op!=INTERP_ALGO_CUBIC and adj_interp!=0)
                                    break;
                                cur_meta.adjInterp=adj_interp;       
                                conf.interpMeta=cur_meta;
                                double cur_ratio=0;
                                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,ALGO_INTERP,TUNING_TARGET_CR,false);
                                cur_ratio=sizeof(T)*8.0/results.first;
                                
                                if (cur_ratio>cur_best_interp_cr){
                                    cur_best_interp_cr=cur_ratio;
                                    
                                    best_meta=cur_meta;

                                }
                            }
                        }
                    }
                }
            }
           
            bestInterpMeta=best_meta;
            conf.interpMeta=best_meta;
            if(conf.autoTuningRate==0){
                if(cur_best_interp_cr>best_interp_cr){
                    
                    best_interp_cr=cur_best_interp_cr;
                    conf.interpMeta=best_meta;
                
                }
            }
        }
        conf.absErrorBound=ori_eb;


        

        if(conf.verbose)           
            printf("Predictor tuning finished.\n");           
        conf.alpha=o_alpha;
        conf.beta=o_beta;
        conf.dims=global_dims;
        conf.num=global_num;
        useInterp= (best_interp_cr>=best_lorenzo_ratio) or best_lorenzo_ratio>=80 or best_interp_cr>=80;//orig 0.95*lorenzo_ratio
        if(conf.verbose){
            if (conf.levelwisePredictionSelection<=0){
                std::cout << "interp best interpAlgo = " << (bestInterpMeta.interpAlgo == 0 ? "LINEAR" : (bestInterpMeta.interpAlgo == 1?"CUBIC":"QUAD")) << std::endl;
                std::cout << "interp best interpParadigm = " << (bestInterpMeta.interpParadigm == 0 ? "1D" : (bestInterpMeta.interpParadigm == 1 ? "MD" : "HD") ) << std::endl;
                if(bestInterpMeta.interpParadigm!=1)
                    std::cout << "interp best direction = " << (unsigned) bestInterpMeta.interpDirection << std::endl;
                if(bestInterpMeta.interpAlgo!=0){
                    std::cout << "interp best cubic spline = " << (unsigned) bestInterpMeta.cubicSplineType << std::endl;
                    std::cout << "interp best adj = " << (unsigned) bestInterpMeta.adjInterp << std::endl;

                }
            }
            else{
                for(int level=conf.levelwisePredictionSelection;level>0;level--){
                    std::cout << "Level: " << (unsigned) level<<std::endl;
                    std::cout << "\tinterp best interpAlgo = " << (bestInterpMeta_list[level-1].interpAlgo == 0 ? "LINEAR" : (bestInterpMeta_list[level-1].interpAlgo == 1 ? "CUBIC" : "QUAD")) << std::endl;
                    std::cout << "\tinterp best interpParadigm = " << (bestInterpMeta_list[level-1].interpParadigm == 0 ? "1D" : (bestInterpMeta_list[level-1].interpParadigm == 1 ? "MD" : "HD") ) << std::endl;
                    if(bestInterpMeta_list[level-1].interpParadigm!=1)
                        std::cout << "\tinterp best direction = " << (unsigned) bestInterpMeta_list[level-1].interpDirection << std::endl;
                    if(bestInterpMeta_list[level-1].interpAlgo!=0){
                        std::cout << "\tinterp best cubic spline = " << (unsigned) bestInterpMeta_list[level-1].cubicSplineType << std::endl;
                        std::cout << "\tinterp best adj = " << (unsigned) bestInterpMeta_list[level-1].adjInterp << std::endl;

                    }
                }
            }
            if(conf.autoTuningRate==0){
                std::cout << "interp best cr = " << best_interp_cr << std::endl;
                printf("choose %s\n", useInterp ? "interp" : "Lorenzo");
            }
        }

    }
    
    else{// if (!conf.blockwiseTuning){ //recently modified. not sure.
        //Timer timer(true);
    
        sampling_data = sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
        if (sampling_num == conf.num) {
            conf.cmprAlgo = ALGO_INTERP;
       
        }

        else{

            


            Config lorenzo_config = conf;
            lorenzo_config.cmprAlgo = ALGO_LORENZO_REG;
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
            lorenzo_config.lorenzo = true;
            lorenzo_config.lorenzo2 = true;
            lorenzo_config.regression = false;
            lorenzo_config.regression2 = false;
            lorenzo_config.openmp = false;
            lorenzo_config.blockSize = 5;//why?

            size_t bufferCap = 1.2 * lorenzo_config * sizeof(T);
            auto buffer = static_cast<uchar *>(malloc(bufferCap));


            size_t sampleOutSize;
            std::vector<T> cur_sampling_data=sampling_data;
            size_t sampleOutSize = SZ_compress_LorenzoReg<T, N>(lorenzo_config, cur_sampling_data.data(), buffer, bufferCap);
                
            double ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
            if(conf.verbose)
                printf("Lorenzo ratio = %.4f\n", ratio);

            best_lorenzo_ratio = ratio;
            double best_interp_ratio = 0;


            for (auto &interp_op: {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC}) {
                ratio = do_not_use_this_interp_compress_block_test<T, N>(
                sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound, interp_op, conf.interpMeta.interpDirection,
                sampling_block, buffer, bufferCap);
                if (ratio > best_interp_ratio) {
                    best_interp_ratio = ratio;
                    conf.interpMeta.interpAlgo = interp_op;
                }
            }
            if(conf.verbose)
                std::cout << "interp best interpAlgo = " << (conf.interpMeta.interpAlgo == 0 ? "LINEAR" : "CUBIC") << std::endl;
                
            int direction_op = factorial(N) - 1;
            ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num,
                                                                 conf.absErrorBound, conf.interpMeta.interpAlgo, direction_op,
                                                                 sampling_block, buffer, bufferCap);
            if (ratio > best_interp_ratio * 1.02) {
                best_interp_ratio = ratio;
                conf.interpMeta.interpDirection = direction_op;
            }
            useInterp=!(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80);
            if(conf.verbose){
                std::cout << "interp best direction = " << (unsigned) conf.interpMeta.interpDirection << std::endl;
                
                printf("Interp ratio = %.4f\n", best_interp_ratio);
                    
                printf("choose %s\n", useInterp ? "interp" : "Lorenzo");
            }
            if (useInterp){
                conf.cmprAlgo=ALGO_INTERP;
            }
            else{
                conf.cmprAlgo=ALGO_LORENZO_REG;
            }

            free(buffer);
        }
        //if(conf.verbose)
        //    timer.stop("sz3 tuning");
    }
    /*
    if(conf.verbose){
        timer.stop("PredTuning");
        timer.start();
    }
    */
    if (useInterp and conf.autoTuningRate>0){
            
        if(conf.verbose)
            std::cout<<"B-M tuning started."<<std::endl;
       
        if (conf.autoTuningRate!=conf.predictorTuningRate){//} and (conf.predictorTuningRate!=0 or conf.autoTuningRate!=conf.waveletTuningRate)){
              
            sampleBlocks<T,N>(data,conf.dims,sampleBlockSize,sampled_blocks,conf.autoTuningRate,conf.profiling,starts,conf.var_first);
        }

        
        double bestalpha=1;
        double bestbeta=1;
        
        double bestb=9999;

        double bestm=0;
        size_t num_sampled_blocks=sampled_blocks.size();
        size_t per_block_ele_num=pow(sampleBlockSize+1,N);
        size_t ele_num=num_sampled_blocks*per_block_ele_num;
        std::vector<double> orig_means;//(num_sampled_blocks,0);
        std::vector<double> orig_sigma2s;//(num_sampled_blocks,0);
        std::vector<double> orig_ranges;//(num_sampled_blocks,0);
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;

        if(conf.tuningTarget==TUNING_TARGET_SSIM){
            size_t ssim_size=conf.SSIMBlockSize;
            for (size_t k =0;k<num_sampled_blocks;k++){

                double orig_mean=0,orig_sigma2=0,orig_range=0;       
                if(N==2){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            std::vector<size_t> starts{i,j};
                            blockwise_profiling<T>(sampled_blocks[k].data(),conf.dims,starts,ssim_size,orig_mean,orig_sigma2,orig_range);
                            orig_means.push_back(orig_mean);
                            orig_sigma2s.push_back(orig_sigma2);
                            orig_ranges.push_back(orig_range);


                        }
                    }
                }
                else if(N==3){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            for (size_t kk=0;kk+ssim_size<sampleBlockSize+1;kk+=ssim_size){
                                std::vector<size_t> starts{i,j,kk};
                                blockwise_profiling<T>(sampled_blocks[k].data(),conf.dims,starts,ssim_size,orig_mean,orig_sigma2,orig_range);
                                orig_means.push_back(orig_mean);
                                orig_sigma2s.push_back(orig_sigma2);
                                orig_ranges.push_back(orig_range);
                            }
                        }
                    }
                }                      
            }

        }
        std::vector<T> flattened_sampled_data;
           
        if (conf.tuningTarget==TUNING_TARGET_AC){
            for(int i=0;i<num_sampled_blocks;i++)
                flattened_sampled_data.insert(flattened_sampled_data.end(),sampled_blocks[i].begin(),sampled_blocks[i].end());

        }
        double oriabseb=conf.absErrorBound;
        

        
        
       
        //std::vector<double> flattened_cur_blocks;

        
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;

        std::vector<double>alpha_list;
        init_alphalist<T,N>(alpha_list,rel_bound,conf);
        size_t alpha_nums=alpha_list.size();
        std::vector<double>beta_list;
        init_betalist<T,N>(beta_list,rel_bound,conf);
        size_t beta_nums=beta_list.size();  
        
        for (size_t i=0;i<alpha_nums;i++){
            for (size_t j=0;j<beta_nums;j++){
                conf.absErrorBound=oriabseb;
                double alpha=alpha_list[i];
                double beta=beta_list[j];
                if (( (alpha>=1 and alpha>beta) or (alpha<0 and beta!=-1) ) )
                    continue;
                conf.alpha=alpha;
                conf.beta=beta; 
                
                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,ALGO_INTERP,(TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                    orig_sigma2s,orig_ranges,flattened_sampled_data);
                double bitrate=results.first;
                double metric=results.second;
                if ( (conf.tuningTarget!=TUNING_TARGET_CR and metric>=bestm and bitrate<=bestb) or (conf.tuningTarget==TUNING_TARGET_CR and bitrate<=bestb ) ){
                    bestalpha=alpha;
                    bestbeta=beta;
                    bestb=bitrate;
                    bestm=metric;
                    useInterp=true;
                }
                else if ( (conf.tuningTarget!=TUNING_TARGET_CR and metric<=bestm and bitrate>=bestb) or (conf.tuningTarget==TUNING_TARGET_CR and bitrate>bestb) ){
                    if ( ((alpha>=1 and pow(alpha,max_interp_level-1)<=beta) or (alpha<1 and alpha*(max_interp_level-1)<=beta)))
                        break;

                    continue;
                }
                else{
                    double eb_fixrate;

                    eb_fixrate=bitrate/bestb;
                    double orieb=conf.absErrorBound;
                    conf.absErrorBound*=eb_fixrate;
                        
                    std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,ALGO_INTERP,(TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                        orig_sigma2s,orig_ranges,flattened_sampled_data);
                    conf.absErrorBound=orieb;

                    double bitrate_r=results.first;
                    double metric_r=results.second;
                    double a=(metric-metric_r)/(bitrate-bitrate_r);
                    double b=metric-a*bitrate;
                    double reg=a*bestb+b;

                    if (reg>bestm){
                        bestalpha=alpha;
                        bestbeta=beta;
                        bestb=bitrate;
                        bestm=metric;
                        useInterp=true;
                    }
                }
                if ( ( (alpha>=1 and pow(alpha,max_interp_level-1)<=beta) or (alpha<1 and alpha*(max_interp_level-1)<=beta)) )
                    break;

            }
        }

        //add lorenzo
        conf.absErrorBound=oriabseb;
        if(conf.testLorenzo){    


            std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,ALGO_LORENZO_REG,(TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                    orig_sigma2s,orig_ranges,flattened_sampled_data);

            double bitrate=results.first;
            double metric=results.second;

            
            if ( (conf.tuningTarget!=TUNING_TARGET_CR and metric>=bestm and bitrate<=bestb) or (conf.tuningTarget==TUNING_TARGET_CR and bitrate<=bestb ) ){
                
                bestb=bitrate;
                bestm=metric;
                bestalpha=-1;
                bestbeta=-1;

                useInterp=false;
                   
            }
            else if ( (conf.tuningTarget!=TUNING_TARGET_CR and metric<=bestm and bitrate>=bestb) or (conf.tuningTarget==TUNING_TARGET_CR and bitrate>bestb) ){
                useInterp=true;
            }
            else{
                double eb_fixrate;

                eb_fixrate=bitrate/bestb;
                double orieb=conf.absErrorBound;
                conf.absErrorBound*=eb_fixrate;                        
                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,ALGO_LORENZO_REG,(TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                    orig_sigma2s,orig_ranges,flattened_sampled_data);
                conf.absErrorBound=orieb;
                double bitrate_r=results.first;
                double metric_r=results.second;
                double a=(metric-metric_r)/(bitrate-bitrate_r);
                double b=metric-a*bitrate;
                double reg=a*bestb+b;

                if (a>0 and reg>bestm){
                    bestb=bitrate;
                    bestm=metric;
                    bestalpha=-1;
                    bestbeta=-1;
                    useInterp=false;
                }
            }          
        }
        conf.absErrorBound=oriabseb;
        /*
        if(conf.verbose){
            timer.stop("B-M step");
            timer.start();
        }
        */
        
        if(conf.tuningTarget==TUNING_TARGET_AC){
            bestm=1-bestm;
        }
        std::string metric_name="Quality";
        if (conf.tuningTarget==TUNING_TARGET_RD ){
            metric_name="PSNR";
        }
        else if (conf.tuningTarget==TUNING_TARGET_SSIM ){
            metric_name="SSIM";
        }
        else if (conf.tuningTarget==TUNING_TARGET_AC ){
            metric_name="AutoCorrelation";
        }
        if(conf.verbose){
            printf("Autotuning finished.\n");
            if (useInterp)
                printf("Interp selected. Selected alpha: %f. Selected beta: %f. Best bitrate: %f. Best %s: %f.\n",bestalpha,bestbeta,bestb, const_cast<char*>(metric_name.c_str()),bestm);
            else
                printf("Lorenzo selected. Best bitrate: %f. Best %s: %f.\n",bestb, const_cast<char*>(metric_name.c_str()),bestm);

        }
        conf.alpha=bestalpha;
        conf.beta=bestbeta;
        conf.dims=global_dims;
        conf.num=global_num;  

    }
    else if(useInterp and conf.QoZ){
        std::pair<double,double> ab=setABwithRelBound(rel_bound,2);
        conf.alpha=ab.first;
        conf.beta=ab.second;


    }
    conf.blockwiseTuning=blockwiseTuning;
    if (useInterp){
        conf.cmprAlgo=ALGO_INTERP;
    }
    else{
         conf.cmprAlgo=ALGO_LORENZO_REG;
    } 
        
    for(int i=0;i<sampled_blocks.size();i++){
        std::vector< T >().swap(sampled_blocks[i]);              
    }
    std::vector< std::vector<T> >().swap(sampled_blocks);     
    return best_lorenzo_ratio;    
}



template<class T, QoZ::uint N>
size_t SZ_compress_Interp_lorenzo(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(conf.cmprAlgo == QoZ::ALGO_INTERP_LORENZO);
    QoZ::calAbsErrorBound(conf, data);

    if(N!=2&&N!=3){
        conf.autoTuningRate=0;
        conf.predictorTuningRate=0;
        conf.maxStep=0;
        conf.levelwisePredictionSelection=0;
        conf.multiDimInterp=0;
        conf.QoZ=0;
        

    }//merge this part to other branches
   
   
    
    
  
    


    if(conf.verbose)
        std::cout << "====================================== BEGIN TUNING ================================" << std::endl;
    QoZ::Timer timer(true);
    
     
    double best_lorenzo_ratio=Tuning<T,N>(conf,data);
    


    if (conf.cmprAlgo == ALGO_INTERP) {
    
        double tuning_time = timer.stop();
        if(conf.verbose){
            std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
            std::cout << "====================================== END TUNING ======================================" << std::endl;
        }

            
        return SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);        

    } 

    else {
        Config lorenzo_config = conf;
        size_t sampling_num, sampling_block;        
        std::vector<size_t> sample_dims(N);
        std::vector<T> sampling_data;

        size_t sampleOutSize;
        double ratio;  
            
        sampling_data = sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);    
        lorenzo_config.cmprAlgo = ALGO_LORENZO_REG;
        
        lorenzo_config.lorenzo = true;
        lorenzo_config.lorenzo2 = true;
        lorenzo_config.regression = false;
        lorenzo_config.regression2 = false;
        lorenzo_config.openmp = false;
        lorenzo_config.blockSize = 5;//why?
        if (sampling_num != conf.num) {
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());

             size_t bufferCap = 1.2 * lorenzo_config * sizeof(T);
            auto buffer = static_cast<uchar *>(malloc(bufferCap));
       
        //lorenzo_config.quantbinCnt = 65536 * 2;
                    
            if(conf.autoTuningRate>0 or conf.predictorTuningRate>0){
                auto sampleOutSize = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), buffer, bufferCap);
                //delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                //printf("Lorenzo ratio = %.2f\n", ratio);
                best_lorenzo_ratio = ratio;
            }
          
            //further tune lorenzo
            if (N == 3 ) {
                lorenzo_config.quantbinCnt = optimize_quant_invl_3d<T>(data, conf.dims[0], conf.dims[1], conf.dims[2], conf.absErrorBound);
                lorenzo_config.pred_dim = 2;
                auto sampleOutSize= SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), buffer, bufferCap);
                //delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                //printf("Lorenzo, pred_dim=2, ratio = %.4f\n", ratio);
                if (ratio > best_lorenzo_ratio * 1.02) {
                    best_lorenzo_ratio = ratio;
                } else {
                    lorenzo_config.pred_dim = 3;
                }
            }
            if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
                auto quant_num = lorenzo_config.quantbinCnt;
                lorenzo_config.quantbinCnt = 16384;
                auto sampleOutSize = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), buffer, bufferCap);
                //delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
    //            printf("Lorenzo, quant_bin=8192, ratio = %.2f\n", ratio);
                if (ratio > best_lorenzo_ratio * 1.02) {
                    best_lorenzo_ratio = ratio;
                } else {
                    lorenzo_config.quantbinCnt = quant_num;
                }
            }
            lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
            free(buffer);
        }
     
        
        conf = lorenzo_config;
        double tuning_time = timer.stop();
        if(conf.verbose){
            std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
            std::cout << "====================================== END TUNING ======================================" << std::endl;
        }
        return SZ_compress_LorenzoReg<T, N>(conf, data, cmpData.cmpCap);
    }
  
}



}  // namespace QoZ
}  // namespace SZ3
#endif
