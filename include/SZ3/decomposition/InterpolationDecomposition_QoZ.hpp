#ifndef _SZ_INTERPOLATION_DECOMPOSITION_QOZ_HPP
#define _SZ_INTERPOLATION_DECOMPOSITION_QOZ_HPP

#include <cmath>
#include <cstring>
#include <limits>
#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/Sample.hpp"

namespace SZ3 {

namespace QoZ {
    template<class T, uint N, class Quantizer>
    class QoZInterpolationDecomposition : public SZ3::concepts::DecompositionInterface<T, int, N> {//added heritage
    public:


        QoZInterpolationDecomposition(const Config &conf, Quantizer quantizer) :
                quantizer(quantizer) {

            static_assert(std::is_base_of<concepts::QuantizerInterface<T, int>, Quantizer>::value,
                      "must implement the quantizer interface");
        }


        T *decompress(const Config &conf, std::vector<int> &quant_inds, T *decData) override  {
            
            init();   
          
            //Timer timer(true);
            this->quant_inds = quant_inds.data();
            //timer.stop("decode");
            //timer.start();
            double eb = quantizer.get_eb();
            if(!anchor){
                *decData = quantizer.recover(0, quant_inds[quant_index++]);
            }
            
            else{
                recover_grid(decData,global_dimensions,maxStep,frozen_dim);                   
                interpolation_level--;           
            }
            size_t meta_index=0;//,coeff_idx=0;
            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {

                if (alpha<0) {
                    if (level >= 3) {
                        quantizer.set_eb(eb * eb_ratio);
                    } else {
                        quantizer.set_eb(eb);
                    }
                }
                
                else if (alpha>=1){
                    
                    
                    double cur_ratio=pow(alpha,level-1);
                    if (cur_ratio>beta){
                        cur_ratio=beta;
                    }
                    
                    quantizer.set_eb(eb/cur_ratio);
                }
                else{
                    
                    
                    double cur_ratio=1-(level-1)*alpha;
                    if (cur_ratio<beta){
                        cur_ratio=beta;
                    }
                   
                    quantizer.set_eb(eb*cur_ratio);
                }

                Interp_Meta cur_meta;
                if(!blockwiseTuning){
                    if (levelwise_predictor_levels==0){
                        cur_meta=interp_meta;
                    }
                    else{

                        if (level-1<levelwise_predictor_levels){

                            cur_meta=interpMeta_list[level-1];
                        }
                        else{

                            cur_meta=interpMeta_list[levelwise_predictor_levels-1];
                        }
                    }
                }
                     
                size_t stride = 1U << (level - 1);
                size_t cur_blocksize;
                
                if (fixBlockSize>0){
                    cur_blocksize=fixBlockSize;
                }
                else{
                    cur_blocksize=blocksize*stride;
                }
                
                auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(decData,std::begin(global_dimensions), std::end(global_dimensions),
                                                           cur_blocksize, 0);//,0);//blockOrder);
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                
                
                for (auto block = inter_begin; block != inter_end; ++block) {
                    if(blockwiseTuning){
                        cur_meta=interpMeta_list[meta_index++];

                    }

                    auto start_idx=block.get_global_index();
                    auto end_idx = start_idx;


                    for (int i = 0; i < N; i++) {
                        end_idx[i] += cur_blocksize;
                        if (end_idx[i] > global_dimensions[i] - 1) {
                            end_idx[i] = global_dimensions[i] - 1;
                        }
                    }

                     block_interpolation(decData, block.get_global_index(), end_idx, PB_recover,
                                        interpolators[cur_meta.interpAlgo], cur_meta, stride,0,cross_block);//,cross_block,regressiveInterp);
                
                        


                

                }
               
            }
            quantizer.postdecompress_data();
            return decData;
        }
        std::vector<int> compress(const Config &conf, T *data) override 
        {
            double temp;
            std::vector<size_t> quant_bin_counts;
            return compress_with_tuning(conf, data,0,0,0,temp,quant_bin_counts);
        }

        std::vector<int> compress_with_tuning(const Config &conf, T *data,int tuning,double & predict_error){
            //double temp;
            std::vector<size_t> quant_bin_counts;
            return compress_with_tuning(conf, data,tuning,0,0,predict_error,quant_bin_counts);
        }

        std::vector<int> compress_with_tuning(const Config &conf, T *data,int tuning,double & predict_error,std::vector<size_t> &quant_bin_counts){//compresstest
            //double temp;
            return compress_with_tuning(conf, data,tuning,0,0,predict_error,quant_bin_counts);
        }

       

        // compress given the error bound
        std::vector<int> compress_with_tuning(const Config &conf, T *data,int tuning,int start_level,int end_level,double & predict_error,std::vector<size_t> &quant_bin_counts) {
            //tuning 0: normal compress 1:tuning to return qbins and psnr 2: tuning to return prediction loss
            //Timer timer;
            //timer.start();
            
            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
            blocksize = conf.interpBlockSize;
            maxStep=conf.maxStep;
            blockwiseTuning=conf.blockwiseTuning;
            fixBlockSize=conf.fixBlockSize;
            frozen_dim = conf.frozen_dim;

            interp_meta=conf.interpMeta;

            alpha=conf.alpha;
            beta=conf.beta;
            std::vector<Interp_Meta>interp_metas;
            cross_block=conf.crossBlock;

            interpMeta_list = conf.interpMeta_list;

            
            //std::cout<<"ap1"<<std::endl;

            init();
            std::vector<int> quant_inds_vec(num_elements);
            quant_inds = quant_inds_vec.data();
            //std::cout<<"ap2"<<std::endl;

            if (tuning){
                predict_error=0.0;
                if(tuning == 1)
                    quant_bin_counts=std::vector<size_t>(interpolation_level,0);

            }
            double eb = quantizer.get_eb();
            //std::cout<<eb<<std::endl;
            if (start_level<=0 or start_level>interpolation_level ){

                start_level=interpolation_level;

                
            } 
            if(end_level>=start_level or end_level<0){
                end_level=0;
            }
            //std::cout<<start_level<<" "<<end_level<< std::endl;

            if(!anchor){
                quant_inds[quant_index++] = quantizer.quantize_and_overwrite(*data, 0);
            }
            else if (start_level==interpolation_level){
                
                build_grid(conf,data,maxStep,tuning);
                if(tuning==1){
                    quant_bin_counts[start_level-1]=quant_index;
                }
                start_level--;
            }
            
           // double predict_error=0.0;
            levelwise_predictor_levels=interpMeta_list.size();
            //std::cout<<"ap3"<<std::endl;

            for (uint level = start_level; level > end_level && level <= start_level; level--) {
                //std::cout<<"ap4 "<<level<<std::endl;

                cur_level=level;
                double cur_eb;
                if (alpha<0) {
                    if (level >= 3) {
                        cur_eb=eb * eb_ratio;
                    } else {
                        cur_eb=eb;
                    }
                }
                else if (alpha>=1){              
                    double cur_ratio=pow(alpha,level-1);
                    if (cur_ratio>beta){
                        cur_ratio=beta;
                    }            
                    cur_eb=eb/cur_ratio;
                }
                else{              
                    double cur_ratio=1-(level-1)*alpha;
                    if (cur_ratio<beta){
                        cur_ratio=beta;
                    }             
                    cur_eb=eb*cur_ratio;
                }
                quantizer.set_eb(cur_eb);
                //std::cout<<"ap4.1 "<<level<<std::endl;

                Interp_Meta cur_meta;
                if (levelwise_predictor_levels==0){
                    cur_meta=interp_meta;
                    
                }
                else{
                    if (level-1<levelwise_predictor_levels){
                        cur_meta=interpMeta_list[level-1];
                    }
                    else{
                        cur_meta=interpMeta_list[levelwise_predictor_levels-1];
                    }
                }
                //std::cout<<"ap4.2 "<<level<<std::endl;
                Interp_Meta cur_level_meta;
                if(blockwiseTuning)
                    cur_level_meta=cur_meta;
                size_t stride = 1U << (level - 1);
                size_t cur_blocksize;
                if (fixBlockSize>0){
                    cur_blocksize=fixBlockSize;
                }
                else{
                    cur_blocksize=blocksize*stride;
                }       
                auto inter_block_range = std::make_shared<
                        multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                           std::end(global_dimensions),
                                                           cur_blocksize, 0);//,0);//conf.blockOrder);
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                //std::cout<<"ap4.3 "<<level<<std::endl;
                for (auto block = inter_begin; block != inter_end; ++block) {
                    auto start_idx=block.get_global_index();
                    auto end_idx = start_idx;
                    for (int i = 0; i < N; i++) {
                        end_idx[i] += cur_blocksize ;
                        if (end_idx[i] > global_dimensions[i] - 1) {
                            end_idx[i] = global_dimensions[i] - 1;
                        }
                    }

                    if(!blockwiseTuning or (N!=2 and N!=3)){
                        predict_error+=block_interpolation(data, start_idx, end_idx, PB_predict_overwrite,
                                    interpolators[cur_meta.interpAlgo],cur_meta, stride,tuning,cross_block);//,cross_block,regressiveInterp);

                    }

                    else{

                        size_t min_len=8;
                        auto start_idx=block.get_global_index();
                        auto end_idx = start_idx;
                        std::array<size_t,N> sample_starts,sample_ends;

                        for (int i = 0; i < N; i++) {
                            end_idx[i] += cur_blocksize ;
                            if (end_idx[i] > global_dimensions[i] - 1) {
                                end_idx[i] = global_dimensions[i] - 1;
                            }

                            double cur_rate=level>=3?1.0:conf.blockwiseSampleRate;//to finetuning
                            size_t  cur_length=(end_idx[i]-start_idx[i])+1,cur_stride=stride*cur_rate;
                            while(cur_stride>stride){
                                if(cur_length/cur_stride>=min_len)
                                    break;
                                cur_stride/=2;
                                cur_rate/=2;
                                if(cur_stride<stride){
                                    cur_stride=stride;
                                    cur_rate=1;
                                }
                            }
            
                            double temp1=0.5-0.5/cur_rate,temp2=0.5+0.5/cur_rate;
                            sample_starts[i]=static_cast<size_t>((temp1*cur_length)/(2*stride))*2*stride+start_idx[i];
                            sample_ends[i]=static_cast<size_t>((temp2*cur_length)/(2*stride))*2*stride+start_idx[i];
                            if(sample_ends[i]>end_idx[i])
                                sample_ends[i]=end_idx[i];

                        }
                        std::vector<T> orig_sampled_block;
                        std::array<size_t,N>sample_strides;
                        for(size_t i=0;i<N;i++)
                            sample_strides[i]=stride;
                        if(frozen_dim>=0)
                            sample_strides[frozen_dim]=1;
                        if(N==2){
                            for(size_t x=sample_starts[0];x<=sample_ends[0] ;x+=sample_strides[0]){
                                //sb_ends[0]++;
                                for(size_t y=sample_starts[1];y<=sample_ends[1];y+=sample_strides[1]){
                                    //sb_ends[1]++;
                                    size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1];
                                    orig_sampled_block.push_back(data[global_idx]);
                                    
                                }
                            }
                        }
                        else if(N==3){
                            for(size_t x=sample_starts[0];x<=sample_ends[0]  ;x+=sample_strides[0]){
                               
                                for(size_t y=sample_starts[1];y<=sample_ends[1] ;y+=sample_strides[1]){
                                    
                                    for(size_t z=sample_starts[2];z<=sample_ends[2] ;z+=sample_strides[2]){
                                       
                                        size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1]+z*dimension_offsets[2];
                                        orig_sampled_block.push_back(data[global_idx]);
                                    }
                                }
                            }
                        } 

                        Interp_Meta best_meta,cur_meta;
                        double best_loss=std::numeric_limits<double>::max();
                        std::vector<uint8_t> interpAlgo_Candidates={cur_level_meta.interpAlgo};
                        std::vector<uint8_t> interpParadigm_Candidates={0};
                        std::vector<uint8_t> cubicSplineType_Candidates={cur_level_meta.cubicSplineType};
                        std::vector<uint8_t> interpDirection_Candidates={0, static_cast<uint8_t>(factorial(N) -1)};
                        if(frozen_dim>=0){
                            if(frozen_dim==0)
                                interpDirection_Candidates={6,7};
                            else if (frozen_dim==1)
                                interpDirection_Candidates={8,9};
                            else
                                interpDirection_Candidates={10,11};
                        }
                        std::vector<uint8_t> adjInterp_Candidates={cur_level_meta.adjInterp};

                        std::vector<double>interp_vars;

                        std::vector<size_t>block_dims(N,0);
                        for (size_t i=0;i<N;i++)
                            block_dims[i]=(sample_ends[i]-sample_starts[i])/stride+1;
                        if(conf.multiDimInterp>0){
                            for(size_t i=1;i<=conf.multiDimInterp;i++)
                                interpParadigm_Candidates.push_back(i);
                            if(conf.dynamicDimCoeff){
                               
                                calculate_interp_error_vars<T,N>(orig_sampled_block.data(),block_dims,interp_vars,cur_level_meta.interpAlgo,cur_level_meta.cubicSplineType,2,1,cur_eb);//cur_eb or 0?
                                preprocess_vars<N>(interp_vars);

                            }
                       
                        }   
                        for (auto &interp_op: interpAlgo_Candidates) {
                            cur_meta.interpAlgo=interp_op;
                            for (auto &interp_pd: interpParadigm_Candidates) {
                                if(frozen_dim>=0 and interp_pd>1)
                                    continue;
                                cur_meta.interpParadigm=interp_pd;

                                for (auto &interp_direction: interpDirection_Candidates) {
                                    if (frozen_dim<0 and (interp_pd==1 or  (interp_pd==2 and N<=2)) and interp_direction!=0)
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

                                            if(conf.dynamicDimCoeff){
                                                for(size_t i=0;i < N;i++){
                                                    cur_meta.dimCoeffs[i]=interp_vars[i];
                                                }

                                            }
                                            double cur_loss=std::numeric_limits<double>::max();
                                            cur_loss=block_interpolation(data, sample_starts, sample_ends, PB_predict_overwrite,
                                                                interpolators[cur_meta.interpAlgo],cur_meta, stride,2,cross_block);//,cross_block,regressiveInterp);

                                            if(cur_loss<best_loss){
                                                best_loss=cur_loss;
                                                best_meta=cur_meta;
                                            }
                                            size_t local_idx=0;
                                            
                                            if(N==2){
                                                for(size_t x=sample_starts[0];x<=sample_ends[0] ;x+=stride){
                                                    
                                                    for(size_t y=sample_starts[1];y<=sample_ends[1];y+=stride){
                                                       
                                                        size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1];
                                                        data[global_idx]=orig_sampled_block[local_idx++];
                                                        
                                                    }
                                                }
                                            }
                                            else if(N==3){
                                                
                                                for(size_t x=sample_starts[0];x<=sample_ends[0]  ;x+=sample_strides[0]){
                                                   
                                                    for(size_t y=sample_starts[1];y<=sample_ends[1] ;y+=sample_strides[1]){
                                                      
                                                        for(size_t z=sample_starts[2];z<=sample_ends[2] ;z+=sample_strides[2]){
                                                          
                                                            size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1]+z*dimension_offsets[2];
                                                            data[global_idx]=orig_sampled_block[local_idx++];
                                                        }
                                                    }
                                                }
                                            } 
                                            cur_meta.dimCoeffs={1.0/3.0,1.0/3.0,1.0/3.0};
                                            

                                            
                                        }
                                    }
                                }
                            }
                        }
                        interp_metas.push_back(best_meta);

                        predict_error+=block_interpolation(data, start_idx, end_idx, PB_predict_overwrite,
                                        interpolators[best_meta.interpAlgo],best_meta, stride,tuning,cross_block);//,cross_block,regressiveInterp);
                    }

                }
                //std::cout<<"ap4.4 "<<level<<std::endl;
                if(tuning==1){
                    quant_bin_counts[level-1]=quant_index;
                }
            }

            if (!interp_metas.empty()){
                interpMeta_list = interp_metas;

            }                    
            //timer.start();

            quantizer.set_eb(eb);

            //if (tuning){
                //conf.quant_bins=quant_inds;
                //std::vector<int>().swap(quant_inds);
                //decomp_square_error=predict_error;
                //size_t bufferSize = 1;
                //uchar *buffer = new uchar[bufferSize];
                //buffer[0]=0;
                //return buffer;
           // }

            //if(conf.verbose)
            //    timer.stop("prediction");//can remove later
            
            //timer.start();
            return quant_inds_vec;
        }

        
        
        
        /*
        uchar *encoding_lossless(size_t &compressed_size,const std::vector<int> &q_inds=std::vector<int>()) {

            if(q_inds.size()>0)
                quant_inds=q_inds;
            size_t bufferSize = 2.5 * (quant_inds.size() * sizeof(T) + quantizer.size_est());//original is 3
            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;
            quantizer.save(buffer_pos);
            quantizer.clear();
            quantizer.postcompress_data();
            //timer.start();
            encoder.preprocess_encode(quant_inds, 0);
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();
//            timer.stop("Coding");
            assert(buffer_pos - buffer < bufferSize);
            //timer.start();
            uchar *lossless_data = lossless.compress(buffer,
                                                     buffer_pos - buffer,
                                                     compressed_size);
            lossless.postcompress_data(buffer);
//            timer.stop("Lossless");
            return lossless_data;

        }
        void set_eb(double eb){
            quantizer.set_eb(eb);
        }
        */ 


       

        void save(uchar *&c) override {

            write(global_dimensions.data(), N, c);
            write(blocksize, c);
            write(interp_meta, c);
            write(alpha,c);
            write(beta,c);
            write(maxStep,c);
            write(levelwise_predictor_levels,c);
            write(blockwiseTuning,c);
            write(fixBlockSize,c);
            write(frozen_dim,c);
            write(cross_block,c);
            //write(conf.regressiveInterp,buffer_pos);
            if(blockwiseTuning){
                size_t meta_num=interpMeta_list.size();
                write(meta_num,c);
                write(interpMeta_list.data(),meta_num,c);
            }
            else if(levelwise_predictor_levels>0){
                write(interpMeta_list.data(),levelwise_predictor_levels,c);
               
            }

            quantizer.save(c);
        }

        void load(const uchar *&c, size_t &remaining_length) override {
            read(global_dimensions.data(), N, c, remaining_length);        
            read(blocksize, c, remaining_length);
            read(interp_meta, c, remaining_length);    
            read(alpha,c,remaining_length);
            read(beta,c,remaining_length);
            read(maxStep,c,remaining_length);
            read(levelwise_predictor_levels,c, remaining_length);
            read(blockwiseTuning,c, remaining_length);
            read(fixBlockSize,c, remaining_length);
            //int frozen_dim=-1;
            read(frozen_dim,c, remaining_length);
            //int cross_block=0;
            read(cross_block,c, remaining_length);
            //int regressiveInterp;   
            //read(regressiveInterp,c, remaining_length);     
            if(blockwiseTuning){
                size_t meta_num;
                read(meta_num,c, remaining_length);
                interpMeta_list.resize(meta_num);
                read(interpMeta_list.data(),meta_num,c, remaining_length);
            }

            else if(levelwise_predictor_levels>0){
                  interpMeta_list.resize(levelwise_predictor_levels);
                read(interpMeta_list.data(),levelwise_predictor_levels,c, remaining_length);
            }   

            quantizer.load(c, remaining_length);
        }

        std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }        

    private:

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover
        };
        
        void init() {
            quant_index = 0;
            assert(blocksize % 2 == 0 && "Interpolation block size should be even numbers");
            num_elements = 1;

            interpolation_level = -1;

            for (int i = 0; i < N; i++) {
                if (interpolation_level < ceil(log2(global_dimensions[i]))) {
                    interpolation_level = static_cast<uint> (ceil(log2(global_dimensions[i])) );
                }
                num_elements *= global_dimensions[i];
            }
            if (maxStep>0){
                anchor=true;//recently moved out of if
                int max_interpolation_level=static_cast<uint> (log2(maxStep))+1;
                if (max_interpolation_level<=interpolation_level){ 
                    interpolation_level=max_interpolation_level;
                }
            }
            dimension_offsets[N - 1] = 1;
            for (int i = N - 2; i >= 0; i--) {
                dimension_offsets[i] = dimension_offsets[i + 1] * global_dimensions[i + 1];
            }
            dimension_sequences = std::vector<std::array<size_t, N>>();
            auto sequence = std::array<size_t, N>();
            for (size_t i = 0; i < N; i++) {
                sequence[i] = i;
            }
            do {
                dimension_sequences.push_back(sequence);
            } while (std::next_permutation(sequence.begin(), sequence.end()));  
            
            
        }
       
        void build_grid(const Config &conf, T *data,size_t maxStep,int tuning=0){
            
            assert(maxStep>0);

           
            if(tuning>1)
                return;

            if (N==1){
                for (size_t x=maxStep*(tuning==1);x<conf.dims[0];x+=maxStep){

                    quantizer.insert_unpred(*(data+x));
                    quant_inds[quant_index++] = 0;
                       
                    
                }
            }

            else if (N==2){
                for (size_t x=maxStep*(tuning==1);x<conf.dims[0];x+=maxStep){
                    for (size_t y=maxStep*(tuning==1);y<conf.dims[1];y+=maxStep){

                        quantizer.insert_unpred(*(data+x*conf.dims[1]+y));
                        quant_inds[quant_index++] = 0;

                    }
                }
            }
            else if(N==3){

                std::array<size_t,3>anchor_strides={maxStep,maxStep,maxStep};

                int fd=conf.frozen_dim;
                if(fd>=0)
                    anchor_strides[fd]=1;
                

                
                for (size_t x=anchor_strides[0]*(tuning==1);x<conf.dims[0];x+=anchor_strides[0]){
                    for (size_t y=anchor_strides[1]*(tuning==1);y<conf.dims[1];y+=anchor_strides[1]){
                        for(size_t z=anchor_strides[2]*(tuning==1);z<conf.dims[2];z+=anchor_strides[2]){
                            quantizer.insert_unpred(*(data+x*dimension_offsets[0]+y*dimension_offsets[1]+z) );
                            quant_inds[quant_index++] = 0;
                        }           
                    }
                }
            }
            else if(N==4){

                std::array<size_t,4>anchor_strides={maxStep,maxStep,maxStep,maxStep};

                int fd=conf.frozen_dim;
                if(fd>=0)
                    anchor_strides[fd]=1;
                

                
                for (size_t x=anchor_strides[0]*(tuning==1);x<conf.dims[0];x+=anchor_strides[0]){
                    for (size_t y=anchor_strides[1]*(tuning==1);y<conf.dims[1];y+=anchor_strides[1]){
                        for(size_t z=anchor_strides[2]*(tuning==1);z<conf.dims[2];z+=anchor_strides[2]){
                            for(size_t w=anchor_strides[3]*(tuning==1);w<conf.dims[3];w+=anchor_strides[3]){
                                quantizer.insert_unpred(*(data+x*dimension_offsets[0]+y*dimension_offsets[1]+z*dimension_offsets[2]+w) );
                                quant_inds[quant_index++] = 0;
                            }
                        }           
                    }
                }
            }

        }
 
        void recover_grid(T *decData,const std::array<size_t,N>& global_dimensions,size_t maxStep,size_t frozen_dim=-1){
            assert(maxStep>0);
            if (N==1){
                for (size_t x=0;x<global_dimensions[0];x+=maxStep){
                    decData[x]=quantizer.recover_unpred();
                    quant_index++;
                }
            }
            else if (N==2){
                for (size_t x=0;x<global_dimensions[0];x+=maxStep){
                    for (size_t y=0;y<global_dimensions[1];y+=maxStep){
                        decData[x*dimension_offsets[0]+y]=quantizer.recover_unpred();
                        quant_index++;
                    }
                }
            }
            else if(N==3){
                std::array<size_t,3>anchor_strides={maxStep,maxStep,maxStep};
                if(frozen_dim>=0)
                    anchor_strides[frozen_dim]=1;
                for (size_t x=0;x<global_dimensions[0];x+=anchor_strides[0]){
                    for (size_t y=0;y<global_dimensions[1];y+=anchor_strides[1]){
                        for(size_t z=0;z<global_dimensions[2];z+=anchor_strides[2]){
                            decData[x*dimension_offsets[0]+y*dimension_offsets[1]+z]=quantizer.recover_unpred();
                            quant_index++;
                        }    
                    }
                }

            }
            else if(N==4){
                std::array<size_t,4>anchor_strides={maxStep,maxStep,maxStep,maxStep};
                if(frozen_dim>=0)
                    anchor_strides[frozen_dim]=1;
                for (size_t x=0;x<global_dimensions[0];x+=anchor_strides[0]){
                    for (size_t y=0;y<global_dimensions[1];y+=anchor_strides[1]){
                        for(size_t z=0;z<global_dimensions[2];z+=anchor_strides[2]){
                            for(size_t w=0;w<global_dimensions[3];w+=anchor_strides[3]){
                                decData[x*dimension_offsets[0]+y*dimension_offsets[1]+z*dimension_offsets[2]+w]=quantizer.recover_unpred();
                                quant_index++;
                            }
                        }    
                    }
                }

            }
        }
        inline void quantize(size_t idx, T &d, T pred) {

                quant_inds[quant_index++] = (quantizer.quantize_and_overwrite(d, pred));
        }
        /*
        inline double quantize_tuning(size_t idx, T &d, T pred, int mode=1) {

            if (mode==1){
                T orig=d;
                quant_inds[quant_index++] = (quantizer.quantize_and_overwrite(d, pred));
                return (d-orig)*(d-orig);

            }
            else{
                double pred_error=fabs(d-pred);
                quant_inds[quant_index++] = quantizer.quantize_and_overwrite(d, pred,false);
                return pred_error;
            }
       
        } */

        inline void recover(size_t idx, T &d, T pred) {
            d = quantizer.recover(pred, quant_inds[quant_index++]);
        }

        inline double quantize_integrated(size_t idx, T &d, T pred, int mode=0){

            double pred_error=0;
            if(mode==-1){//recover
                d = quantizer.recover(pred, quant_inds[quant_index++]);
                return 0;
            }
            else if(mode==0){
                 quant_inds[quant_index++] = (quantizer.quantize_and_overwrite(d, pred));
                 return 0;
            }
            else if(mode==1){
                T orig=d;
                quant_inds[quant_index++] = (quantizer.quantize_and_overwrite(d, pred));
                return (d-orig)*(d-orig);
            }
            else{// if (mode==2){
                pred_error=fabs(d-pred);
                quantizer.quantize_and_overwrite(d, pred);//todo: check whether to set save_unpred = false
                return pred_error;
            }
        }
        


        double block_interpolation_1d_crossblock(T *data, const std::array<size_t,N> &begin_idx, const std::array<size_t,N> &end_idx,const size_t &direction,const size_t &math_stride, const std::string &interp_func, const PredictorBehavior pb,const Interp_Meta &meta,int cross_block=1,int tuning=0) {//cross block: 0: no cross 1: only front-cross 2: all cross
            size_t math_begin_idx=begin_idx[direction],math_end_idx=end_idx[direction];
            size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;

            if (n <= 1) {
                return 0;
            }
            double predict_error = 0;
            bool cross_back=cross_block>0;
            bool cross_front=cross_block>0;
            if(cross_front){
                for(size_t i=0;i<N;i++){
                    if(i!=direction and begin_idx[i]%(2*math_stride)!=0){
                        cross_front=false;
                        break;
                    }
                }
            }

            size_t begin=0,global_end_idx=global_dimensions[direction];
            for(size_t i=0;i<N;i++)
                begin+=dimension_offsets[i]*begin_idx[i];

            size_t stride=math_stride*dimension_offsets[direction];
            size_t stride2x = 2 * stride;
            

            int mode=(pb == PB_predict_overwrite)?tuning:-1;
          //  size_t quant_idx=quant_index;

            if (interp_func == "linear") {
                for (size_t i = 1; i + 1 < n; i += 2) {
                    T *d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if(cross_front and math_end_idx+math_stride<global_end_idx)
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                    else if (n < 3) {                              
                        predict_error+=quantize_integrated(d - data, *d, *(d - stride),mode);
                        } else {
                        predict_error+=quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                    }
                }

            } 

            else {

                size_t stride3x = 3 * stride;
                size_t math_stride2x=2*math_stride;
                size_t math_stride3x=3*math_stride;
                T *d;
                size_t i;
                if(!meta.adjInterp){
                    size_t i_start= (cross_back and math_begin_idx>=math_stride2x)?1:3;

                    for (i = i_start; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        predict_error+=quantize_integrated(d - data, *d,
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                        
                    }


                    std::vector<size_t> boundary;
                    if(i_start==3 or n<=4)
                        boundary.push_back(1);

                    if(n%2==1){
                        if(n>3)
                            boundary.push_back(n-2);
                    }
                    else{
                        if(n>4)
                            boundary.push_back(n-3);
                        if (n>2)
                            boundary.push_back(n-1);
                    }

                    for(auto i:boundary){
                        d = data + begin + i*stride;
                        size_t math_cur_idx=math_begin_idx+i*math_stride;
                        if( i>=3 or (cross_back and math_cur_idx>=math_stride3x) ){
                            if(i+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) ){
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                
                            }
                            else if(i+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx )){
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                                
                            }
                            else {

                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_linear1(*(d - stride3x), *(d - stride)),mode);
                                
                            }
                        }
                        else{
                            if(i+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) ){
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                
                            }
                            else if(i+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) ){
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_linear( *(d - stride), *(d + stride)),mode);
                                
                            }
                            else {
                                predict_error+=quantize_integrated(d - data, *d,
                                        *(d - stride),mode);
                                
                            }
                        }
                        
                    }
                }
                else{
                    size_t i_start= (cross_back and math_begin_idx>=math_stride2x)?1:5;
                    for (i = i_start; i + 3 < n; i += 4) {
                        
                        d = data + begin + i * stride;
                     
                        predict_error+=quantize_integrated(d - data, *d,
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                    }


                    std::vector<size_t> boundary;
                    if(n<=4){
                        boundary.push_back(1);
                    }
                    else{
                        if(i_start==5)
                            boundary.push_back(1);
                        int temp=n%4;
                        if(temp==0)
                            temp=4;
                        if(temp!=1){
                            boundary.push_back(n+1-temp);
                        }
                    }

                    for(auto i:boundary){
                        d = data + begin + i*stride;
                        size_t math_cur_idx=math_begin_idx+i*math_stride;
                        if(i>3 or (cross_back and math_cur_idx>=math_stride3x)){
                            if(i+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) )
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                            else if(i+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                            else{
                                //std::cout<<"dwa"<<i<<std::endl;
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_linear1(*(d - stride3x), *(d - stride)),mode);
                            }
                        }
                        else{
                            if(i+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) )
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                            else if(i+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_linear( *(d - stride), *(d + stride)),mode);
                            else 
                                predict_error+=quantize_integrated(d - data, *d,
                                        *(d - stride),mode);
                        }
                    }

                    for (i = 3; i + 3 < n; i += 4) {
                        d = data + begin + i * stride;

                        predict_error+=quantize_integrated(d - data, *d,
                                   interp_cubic_adj(meta.cubicSplineType,*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x),*(d + stride3x)),mode);
                        //predict_error+=quantize_integrated(d - data, *d,
                         //           interp_cubic_3(*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x)),mode);
                    }

                    size_t temp=n%4;
                    if(temp!=3 and n>temp+1){

                        i=n-1-temp;
                        
                       
                        d = data + begin + i*stride;
                        size_t math_cur_idx=math_begin_idx+i*math_stride;
 
                        if(i>3 or (cross_back and math_cur_idx>=math_stride3x)){
                            if(i+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx)   )
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_cubic_adj2(meta.cubicSplineType,*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d + stride3x)),mode);//to determine,another choice is noadj cubic.
                            else if(i+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx )  )
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_quad_2_adj(*(d - stride2x), *(d - stride), *(d + stride)),mode);
                            else {
                                predict_error+=quantize_integrated(d - data, *d,
                                        lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                            }
                        }
                        else{
                            if(i+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx)  )
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                            else if(i+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                predict_error+=quantize_integrated(d - data, *d,
                                        interp_linear( *(d - stride), *(d + stride)),mode);
                            else 
                                predict_error+=quantize_integrated(d - data, *d,
                                        *(d - stride),mode);
                        }
                    }
                    
                }
       
            }
           // quant_index=quant_idx;
            return predict_error;
        }
      
        double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride,const std::string &interp_func,const PredictorBehavior pb,const Interp_Meta &meta,int tuning=0) {
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }
            double predict_error = 0;
            size_t stride2x=2*stride;
            size_t stride3x = 3 * stride;
            size_t stride5x = 5 * stride;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
           // size_t quant_idx=quant_index;
            if (interp_func == "linear" || n < 5) {
                for (size_t i = 1; i + 1 < n; i += 2) {
                    T *d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if (n < 3) {                              
                        predict_error+=quantize_integrated(d - data, *d, *(d - stride),mode);
                        } else {
                        predict_error+=quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                    }
                }

            } 
            
            else {
                T *d;
                size_t i;
                if(!meta.adjInterp){
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        predict_error+=quantize_integrated(d - data, *d,
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                    }

                    d = data + begin + stride;

                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),mode);
                    d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        predict_error+=
                        quantize_integrated(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)),mode);
                    }
                }
                else{
                    //i=1
                    d = data + begin + stride;
                    
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),mode);
                    //predict_error+=quantize_integrated(d - data, *d, interp_cubic_front_adj(*(d -stride),*(d + stride), *(d+stride2x), *(d + stride3x)),mode);
                    for (i = 5; i + 3 < n; i += 4) {
                        
                        d = data + begin + i * stride;
                     
                        predict_error+=quantize_integrated(d - data, *d,
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                    }

                    //i=n-3 or n-2

                    if(i<n-1){
                        d = data + begin + i * stride;
                       
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                        //predict_error+=quantize_integrated(d - data, *d, interp_cubic_back_adj(*(d -stride3x),*(d - stride2x), *(d-stride), *(d + stride)),mode);
                    }
                    //i=n-1
                    else if(i<n){
                         d = data + begin + (n - 1) * stride;
             
                        predict_error+=
                        //quantize_integrated(d - data, *d, interp_quad_3_adj(*(d - stride3x), *(d - stride2x), *(d - stride)),mode);
                        quantize_integrated(d - data, *d, *(d - stride),mode);//to determine
                        //quantize_integrated(d - data, *d, *(d - stride),mode);

                    }


                    for (i = 3; i + 3 < n; i += 4) {
                        d = data + begin + i * stride;

                        predict_error+=quantize_integrated(d - data, *d,
                                   interp_cubic_adj(meta.cubicSplineType,*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x),*(d + stride3x)),mode);
                        //predict_error+=quantize_integrated(d - data, *d,
                         //           interp_cubic_3(*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x)),mode);
                    }
                    //i=n-3 or n-2
                    if(i<n-1){
                        d = data + begin + i * stride;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x), *(d - stride), *(d + stride)),mode);
                        //predict_error+=quantize_integrated(d - data, *d, interp_cubic_back_adj(*(d -stride3x),*(d - stride2x), *(d-stride), *(d + stride)),mode);
                    }
                    //i=n-1
                    else if(i<n){
                         d = data + begin + (n - 1) * stride;

                        predict_error+=
                        //quantize_integrated(d - data, *d, interp_quad_3_adj(*(d - stride3x), *(d - stride2x), *(d - stride)),mode);
                        quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x),*(d - stride)) ,mode);//to determine
                        //quantize_integrated(d - data, *d, *(d - stride),mode);
                    }
                }
                
            }
           // quant_index=quant_idx;
            return predict_error;
        } 
        double block_interpolation_1d_crossblock_2d(T *data, const std::array<size_t,N> &begin_idx, const std::array<size_t,N> &end_idx,const size_t &direction, std::array<size_t,N> &steps,const size_t &math_stride, const std::string &interp_func, const PredictorBehavior pb,const Interp_Meta &meta,int cross_block=1,int tuning=0) {//cross block: 0: no cross 1: only front-cross 2: all cross
            
            for(size_t i=0;i<N;i++){
                if(end_idx[i]<begin_idx[i])
                    return 0;
            }

            size_t math_begin_idx=begin_idx[direction],math_end_idx=end_idx[direction];
            size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;
            if (n <= 1) {
                return 0;
            }
          //  size_t quant_idx=quant_index;

            double predict_error = 0;
            bool cross_back=cross_block>0;
            bool global_cross_front=cross_block>0;
            
            
            
            size_t begin=0,global_end_idx=global_dimensions[direction];
            for(size_t i=0;i<N;i++)
                begin+=dimension_offsets[i]*begin_idx[i];

            size_t stride=math_stride*dimension_offsets[direction];
            std::array<size_t,N>begins,ends,strides;
            for(size_t i=0;i<N;i++){
                begins[i]=0;
                ends[i]=end_idx[i]-begin_idx[i]+1;
                strides[i]=dimension_offsets[i];
            }
            strides[direction]=stride;
            

            size_t stride2x = 2 * stride;
            

            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            

            if (interp_func == "linear") {
                begins[direction]=1;
                ends[direction]=n-1;
                steps[direction]=2;
                for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                    for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                        T *d = data + begin + i * strides[0]+j*strides[1];
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);


                    }

                }
                if (n % 2 == 0) {
                    begins[direction]=n-1;
                    ends[direction]=n;
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            bool cross_front=global_cross_front;
                            if(cross_front){
                                std::array<size_t,N>idxs{begin_idx[0]+i,begin_idx[1]+j};
                                for(size_t t=0;t<N;t++){
                                    if(t!=direction and idxs[t]%(2*math_stride)!=0){
                                        cross_front=false;
                                        break;
                                    }
                                }
                            }
                            T *d = data + begin + i * strides[0]+j*strides[1];
                            
                            if(cross_front and math_end_idx+math_stride<global_end_idx)
                                predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                            else if (n < 3) {                              
                                predict_error+=quantize_integrated(d - data, *d, *(d - stride),mode);
                                } else {
                                predict_error+=quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                            }
                        
                        }
                    }
                    
                }

            } 
            /*
            else if (interp_func == "quad"){
                T *d= data + begin +  stride;
                predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                for (size_t i = 3; i + 1 < n; i += 2) {
                    T *d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x),*(d - stride), *(d + stride)),mode);
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    predict_error+=quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                }


            }*/
            else {
                size_t stride3x = 3 * stride;
                size_t math_stride2x=2*math_stride;
                size_t math_stride3x=3*math_stride;
                T *d;
                //size_t i;

                if(!meta.adjInterp){
                    size_t i_start= (cross_back and math_begin_idx>=math_stride2x)?1:3;

                    begins[direction]=i_start;
                    ends[direction]=(n>=3)?(n-3):0;
                    steps[direction]=2;
                  
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            d = data + begin + i * strides[0]+j*strides[1];
                            predict_error+=quantize_integrated(d - data, *d,
                                        interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                            
                        }
                        
                    }

                  
                    std::vector<size_t> boundary;
                    if(i_start==3 or n<=4)
                        boundary.push_back(1);

                    if(n%2==1){
                        if(n>3)
                            boundary.push_back(n-2);
                    }
                    else{
                        if(n>4)
                            boundary.push_back(n-3);
                        if (n>2)
                            boundary.push_back(n-1);
                    }

                    for(auto ii:boundary){
                       // std::cout<<ii<<std::endl;

                        begins[direction]=ii;
                        ends[direction]=ii+1;

                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                
                                bool cross_front=global_cross_front;
                                if(cross_front){
                                    std::array<size_t,N>idxs{begin_idx[0]+i,begin_idx[1]+j};
                                    for(size_t t=0;t<N;t++){
                                        if(t!=direction and idxs[t]%(2*math_stride)!=0){
                                            cross_front=false;
                                            break;
                                        }
                                    }
                                }
                                d = data + begin + i * strides[0]+j*strides[1];
                                size_t main_idx=ii;
                                size_t math_cur_idx=math_begin_idx+main_idx*math_stride;
                                if( main_idx>=3 or (cross_back and math_cur_idx>=math_stride3x) ){
                                    if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) ){
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                        
                                    }
                                    else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx )){
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                                        
                                    }
                                    else {
                                        //if(mode==0)
                                        //std::cout<<"n-1 "<<main_idx<<std::endl;
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_linear1(*(d - stride3x), *(d - stride)),mode);
                                        
                                    }
                                }
                                else{
                                    if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) ){
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                        
                                    }
                                    else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) ){
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_linear( *(d - stride), *(d + stride)),mode);
                                        
                                    }
                                    else {
                                        predict_error+=quantize_integrated(d - data, *d,
                                                *(d - stride),mode);
                                        
                                    }
                                }
                                
                            }
                            
                        }
                    }
                }
                else{
                    //i=1

                    size_t i_start= (cross_back and math_begin_idx>=math_stride2x)?1:5;
                    begins[direction]=i_start;
                    ends[direction]=(n>=3)?(n-3):0;
                    steps[direction]=4;
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            
                    
                            d = data + begin + i * strides[0]+j*strides[1];
                         
                            predict_error+=quantize_integrated(d - data, *d,
                                        interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                            
                        }
                    }


                    std::vector<size_t> boundary;
                    if(n<=4){
                        boundary.push_back(1);
                    }
                    else{
                        if(i_start==5)
                            boundary.push_back(1);
                        int temp=n%4;
                        if(temp==0)
                            temp=4;
                        if(temp!=1){
                            boundary.push_back(n+1-temp);
                        }
                    }


                    for(auto ii:boundary){

                        begins[direction]=ii;
                        ends[direction]=ii+1;

                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                
                                bool cross_front=global_cross_front;
                                if(cross_front){
                                    std::array<size_t,N>idxs{begin_idx[0]+i,begin_idx[1]+j};
                                    for(size_t t=0;t<N;t++){
                                        if(t!=direction and idxs[t]%(2*math_stride)!=0){
                                            cross_front=false;
                                            break;
                                        }
                                    }
                                }
                                d = data + begin + i * strides[0]+j*strides[1];
                                size_t main_idx=ii;
                   
                                size_t math_cur_idx=math_begin_idx+ main_idx*math_stride;
                                if(main_idx>3 or (cross_back and math_cur_idx>=math_stride3x)){
                                    if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) )
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                    else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                                    else{
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_linear1(*(d - stride3x), *(d - stride)),mode);
                                    }
                                }
                                else{
                                    if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) )
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                    else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_linear( *(d - stride), *(d + stride)),mode);
                                    else 
                                        predict_error+=quantize_integrated(d - data, *d,
                                                *(d - stride),mode);
                                }
                                
                            }
                        }
                    }
                    begins[direction]=3;
                    ends[direction]=(n>=3)?(n-3):0;
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                          
                    
                            d = data + begin + i * strides[0]+j*strides[1];

                            predict_error+=quantize_integrated(d - data, *d,
                                       interp_cubic_adj(meta.cubicSplineType,*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x),*(d + stride3x)),mode);
                            //predict_error+=quantize_integrated(d - data, *d,
                             //           interp_cubic_3(*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x)),mode);
                            
                        }
                    }
                    //if(end_idx[0]==503 )
                     //   std::cout<<"r4"<<std::endl;

                    size_t temp=n%4;
                    if(temp!=3 and n>temp+1){

                        size_t ii =n-1-temp;
                        begins[direction]=n-1-temp;
                        ends[direction]=begins[direction]+1;
                        
                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                
                                bool cross_front=global_cross_front;
                                if(cross_front){
                                    std::array<size_t,N>idxs{begin_idx[0]+i,begin_idx[1]+j};
                                    for(size_t t=0;t<N;t++){
                                        if(t!=direction and idxs[t]%(2*math_stride)!=0){
                                            cross_front=false;
                                            break;
                                        }
                                    }
                                }
                
                                d = data + begin + i * strides[0]+j*strides[1];
                                size_t main_idx=ii;
                                size_t math_cur_idx=math_begin_idx+main_idx*math_stride;
                                if(main_idx>3 or (cross_back and math_cur_idx>=math_stride3x)){
                                    if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx)   )
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_cubic_adj2(meta.cubicSplineType,*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d + stride3x)),mode);//to determine,another choice is noadj cubic.
                                    else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx )  )
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_quad_2_adj(*(d - stride2x), *(d - stride), *(d + stride)),mode);
                                    else {
                                        predict_error+=quantize_integrated(d - data, *d,
                                                lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                                    }
                                }
                                else{
                                    if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx)  )
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                    else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                        predict_error+=quantize_integrated(d - data, *d,
                                                interp_linear( *(d - stride), *(d + stride)),mode);
                                    else 
                                        predict_error+=quantize_integrated(d - data, *d,
                                                *(d - stride),mode);
                                }
                                
                            }
                        }      
                    
                    }
                }
            }

          //  quant_index=quant_idx;
            return predict_error;
        }
        double block_interpolation_1d_crossblock_3d(T *data, const std::array<size_t,N> &begin_idx, const std::array<size_t,N> &end_idx,const size_t &direction, std::array<size_t,N> &steps,const size_t &math_stride, const std::string &interp_func, const PredictorBehavior pb,const Interp_Meta &meta,int cross_block=1,int tuning=0) {//cross block: 0: no cross 1: only front-cross 2: all cross
            
            for(size_t i=0;i<N;i++){
                if(end_idx[i]<begin_idx[i])
                    return 0;
            }

            size_t math_begin_idx=begin_idx[direction],math_end_idx=end_idx[direction];
            size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;

            if (n <= 1) {
                return 0;
            }
            //size_t quant_idx=quant_index;

            double predict_error = 0;
            bool cross_back=cross_block>0;
            bool global_cross_front=cross_block>0;
            
            
            
            size_t begin=0,global_end_idx=global_dimensions[direction];
            for(size_t i=0;i<N;i++)
                begin+=dimension_offsets[i]*begin_idx[i];

            size_t stride=math_stride*dimension_offsets[direction];
            std::array<size_t,N>begins,ends,strides;
            for(size_t i=0;i<N;i++){
                begins[i]=0;
                ends[i]=end_idx[i]-begin_idx[i]+1;
                strides[i]=dimension_offsets[i];
            }
            strides[direction]=stride;
            

            size_t stride2x = 2 * stride;
            

            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            

            if (interp_func == "linear") {
                begins[direction]=1;
                ends[direction]=n-1;
                steps[direction]=2;
                for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                    for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                        for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                            T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);



                        }
                    }

                }
                if (n % 2 == 0) {
                    begins[direction]=n-1;
                    ends[direction]=n;
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                bool cross_front=global_cross_front;
                                if(cross_front){
                                    std::array<size_t,N>idxs{begin_idx[0]+i,begin_idx[1]+j,begin_idx[2]+k};
                                    for(size_t t=0;t<N;t++){
                                        if(t!=direction and idxs[t]%(2*math_stride)!=0){
                                            cross_front=false;
                                            break;
                                        }
                                    }
                                }
                                T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                                
                                if(cross_front and math_end_idx+math_stride<global_end_idx)
                                    predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                                else if (n < 3) {                              
                                    predict_error+=quantize_integrated(d - data, *d, *(d - stride),mode);
                                    } else {
                                    predict_error+=quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                                }
                            }
                        }
                    }
                    
                }

            } 
            /*
            else if (interp_func == "quad"){
                T *d= data + begin +  stride;
                predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                for (size_t i = 3; i + 1 < n; i += 2) {
                    T *d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x),*(d - stride), *(d + stride)),mode);
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    predict_error+=quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                }


            }*/
            else {
                
                size_t stride3x = 3 * stride;
                size_t math_stride2x=2*math_stride;
                size_t math_stride3x=3*math_stride;
                T *d;
                //size_t i;

                if(!meta.adjInterp){
                    size_t i_start= (cross_back and math_begin_idx>=math_stride2x)?1:3;

                    begins[direction]=i_start;
                    ends[direction]=(n>=3)?(n-3):0;
                    steps[direction]=2;
                  
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                                predict_error+=quantize_integrated(d - data, *d,
                                            interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                            }
                        }
                        
                    }

                  
                    std::vector<size_t> boundary;
                    if(i_start==3 or n<=4)
                        boundary.push_back(1);

                    if(n%2==1){
                        if(n>3)
                            boundary.push_back(n-2);
                    }
                    else{
                        if(n>4)
                            boundary.push_back(n-3);
                        if (n>2)
                            boundary.push_back(n-1);
                    }

                    for(auto ii:boundary){
                        begins[direction]=ii;
                        ends[direction]=ii+1;

                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    bool cross_front=global_cross_front;
                                    if(cross_front){
                                        std::array<size_t,N>idxs{begin_idx[0]+i,begin_idx[1]+j,begin_idx[2]+k};
                                        for(size_t t=0;t<N;t++){
                                            if(t!=direction and idxs[t]%(2*math_stride)!=0){
                                                cross_front=false;
                                                break;
                                            }
                                        }
                                    }
                                    d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                                    size_t main_idx=ii;
                                   size_t math_cur_idx=math_begin_idx+main_idx*math_stride;
                                   if( main_idx>=3 or (cross_back and math_cur_idx>=math_stride3x) ){
                                        if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) ){
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                            
                                        }
                                        else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx )){
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                                            
                                        }
                                        else {
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_linear1(*(d - stride3x), *(d - stride)),mode);
                                            
                                        }
                                    }
                                    else{
                                        if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) ){
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                            
                                        }
                                        else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) ){
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_linear( *(d - stride), *(d + stride)),mode);
                                            
                                        }
                                        else {
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    *(d - stride),mode);
                                            
                                        }
                                    }
                                }
                            }
                            
                        }
                    }
                }
                else{
                    //i=1
                    size_t i_start= (cross_back and math_begin_idx>=math_stride2x)?1:5;
                    begins[direction]=i_start;
                    ends[direction]=(n>=3)?(n-3):0;
                    steps[direction]=4;
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                    
                                d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                             
                                predict_error+=quantize_integrated(d - data, *d,
                                            interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                            }
                        }
                    }


                    std::vector<size_t> boundary;
                    if(n<=4){
                        boundary.push_back(1);
                    }
                    else{
                        if(i_start==5)
                            boundary.push_back(1);
                        int temp=n%4;
                        if(temp==0)
                            temp=4;
                        if(temp!=1){
                            boundary.push_back(n+1-temp);
                        }
                    }

                    for(auto ii:boundary){

                        begins[direction]=ii;
                        ends[direction]=ii+1;

                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    bool cross_front=global_cross_front;
                                    if(cross_front){
                                        std::array<size_t,N>idxs{begin_idx[0]+i,begin_idx[1]+j,begin_idx[2]+k};
                                        for(size_t t=0;t<N;t++){
                                            if(t!=direction and idxs[t]%(2*math_stride)!=0){
                                                cross_front=false;
                                                break;
                                            }
                                        }
                                    }
                                    d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                                    size_t main_idx=ii;
                       
                                    size_t math_cur_idx=math_begin_idx+ main_idx*math_stride;
                                    if(main_idx>3 or (cross_back and math_cur_idx>=math_stride3x)){
                                        if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) )
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                        else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                                        else{
                                            //std::cout<<"dwa"<<i<<std::endl;
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_linear1(*(d - stride3x), *(d - stride)),mode);
                                        }
                                    }
                                    else{
                                        if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx) )
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                        else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_linear( *(d - stride), *(d + stride)),mode);
                                        else 
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    *(d - stride),mode);
                                    }
                                }
                            }
                        }
                    }
                    begins[direction]=3;
                    ends[direction]=(n>=3)?(n-3):0;
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                    
                                d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];

                                predict_error+=quantize_integrated(d - data, *d,
                                           interp_cubic_adj(meta.cubicSplineType,*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x),*(d + stride3x)),mode);
                                //predict_error+=quantize_integrated(d - data, *d,
                                 //           interp_cubic_3(*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x)),mode);
                            }
                        }
                    }

                    size_t temp=n%4;
                    if(temp!=3 and n>temp+1){

                        size_t ii =n-1-temp;
                        begins[direction]=n-1-temp;
                        ends[direction]=begins[direction]+1;
                        
                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    bool cross_front=global_cross_front;
                                    if(cross_front){
                                        std::array<size_t,N>idxs{begin_idx[0]+i,begin_idx[1]+j,begin_idx[2]+k};
                                        for(size_t t=0;t<N;t++){
                                            if(t!=direction and idxs[t]%(2*math_stride)!=0){
                                                cross_front=false;
                                                break;
                                            }
                                        }
                                    }
                    
                                    d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                                    size_t main_idx=ii;
                                    size_t math_cur_idx=math_begin_idx+main_idx*math_stride;
                                    if(main_idx>3 or (cross_back and math_cur_idx>=math_stride3x)){
                                        if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx)   )
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_cubic_adj2(meta.cubicSplineType,*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d + stride3x)),mode);//to determine,another choice is noadj cubic.
                                        else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx )  )
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_quad_2_adj(*(d - stride2x), *(d - stride), *(d + stride)),mode);
                                        else {
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                                        }
                                    }
                                    else{
                                        if(main_idx+3<n or (cross_front and math_cur_idx+math_stride3x<global_end_idx)  )
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_quad_1( *(d - stride), *(d + stride), *(d + stride3x)),mode);
                                        else if(main_idx+1<n or (cross_front and math_cur_idx+math_stride<global_end_idx ) )
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    interp_linear( *(d - stride), *(d + stride)),mode);
                                        else 
                                            predict_error+=quantize_integrated(d - data, *d,
                                                    *(d - stride),mode);
                                    }
                                }
                            }
                        }      
                    
                    }
                }
               
            }

          //  quant_index=quant_idx;
            return predict_error;
        }
        double block_interpolation_1d_regressive(T *data, size_t begin, size_t end, size_t stride,const std::string &interp_func,const PredictorBehavior pb,const Interp_Meta &meta,const std::vector<float>& coeff,int tuning=0) {
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }
            double predict_error = 0;
            size_t stride2x=2*stride;
            size_t stride3x = 3 * stride;
            size_t stride5x = 5 * stride;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
           // size_t quant_idx=quant_index;
           
            if (interp_func == "linear" || n < 5) {
                for (size_t i = 1; i + 1 < n; i += 2) {
                    T *d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if (n < 4) {                              
                        predict_error+=quantize_integrated(d - data, *d, *(d - stride),mode);
                        } else {
                        predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)),mode);
                    }
                }

            } 

            else {
                T *d;
                size_t i;
                if(!meta.adjInterp){

                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        //predict_error+=quantize_integrated(d - data, *d,
                        //            interp_cubic(meta.cubicSplineType,*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                        predict_error+=quantize_integrated(d - data, *d,
                                    coeff[0]*(*(d - stride3x))+coeff[1]*(*(d - stride))+coeff[2]*( *(d + stride))+coeff[3]*( *(d + stride3x)),mode);
                    }

                    d = data + begin + stride;

                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),mode);
                    d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        predict_error+=
                        quantize_integrated(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)),mode);
                    }
                }
                else{
                    //i=1
                    d = data + begin + stride;
                    
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),mode);
                    //predict_error+=quantize_integrated(d - data, *d, interp_cubic_front_adj(*(d -stride),*(d + stride), *(d+stride2x), *(d + stride3x)),mode);
                    for (i = 5; i + 3 < n; i += 4) {
                        
                        d = data + begin + i * stride;
                     
                        predict_error+=quantize_integrated(d - data, *d,
                                     coeff[0]*(*(d - stride3x))+coeff[1]*(*(d - stride))+coeff[2]*( *(d + stride))+coeff[3]*( *(d + stride3x)),mode);
                    }

                    //i=n-3 or n-2

                    if(i<n-1){
                        d = data + begin + i * stride;
                       
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                        //predict_error+=quantize_integrated(d - data, *d, interp_cubic_back_adj(*(d -stride3x),*(d - stride2x), *(d-stride), *(d + stride)),mode);
                    }
                    //i=n-1
                    else if(i<n){
                         d = data + begin + (n - 1) * stride;
             
                        predict_error+=
                        //quantize_integrated(d - data, *d, interp_quad_3_adj(*(d - stride3x), *(d - stride2x), *(d - stride)),mode);
                        quantize_integrated(d - data, *d, *(d - stride),mode);//to determine
                        //quantize_integrated(d - data, *d, *(d - stride),mode);

                    }


                    for (i = 3; i + 3 < n; i += 4) {
                        d = data + begin + i * stride;

                        predict_error+=quantize_integrated(d - data, *d,
                                   interp_cubic_adj(meta.cubicSplineType,*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x),*(d + stride3x)),mode);
                        //predict_error+=quantize_integrated(d - data, *d,
                         //           interp_cubic_3(*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x)),mode);
                    }
                    //i=n-3 or n-2
                    if(i<n-1){
                        d = data + begin + i * stride;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x), *(d - stride), *(d + stride)),mode);
                        //predict_error+=quantize_integrated(d - data, *d, interp_cubic_back_adj(*(d -stride3x),*(d - stride2x), *(d-stride), *(d + stride)),mode);
                    }
                    //i=n-1
                    else if(i<n){
                         d = data + begin + (n - 1) * stride;

                        predict_error+=
                        //quantize_integrated(d - data, *d, interp_quad_3_adj(*(d - stride3x), *(d - stride2x), *(d - stride)),mode);
                        quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x),*(d - stride)) ,mode);//to determine
                        //quantize_integrated(d - data, *d, *(d - stride),mode);
                    }
                }

            }

           // quant_index=quant_idx;
            return predict_error;
        } 


        double block_interpolation_2d(T *data, size_t begin1, size_t end1, size_t begin2, size_t end2, size_t stride1,size_t stride2,const std::string &interp_func,const PredictorBehavior pb,const std::array<float,2> &dim_coeffs,const Interp_Meta &meta,int tuning=0) {
            size_t n = (end1 - begin1) / stride1 + 1;
            if (n <= 1) {
                return 0;
            }
            size_t m = (end2 - begin2) / stride2 + 1;
            if (m <= 1) {
                return 0;
            }


            double predict_error = 0;
            
            float coeff_x=(dim_coeffs[0])/((dim_coeffs[0])+(dim_coeffs[1])),coeff_y=1-coeff_x;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            size_t begin=begin1+begin2,end=end1+end2;
           // size_t quant_idx=quant_index;
            
            if (interp_func == "linear"||(n<5 &m<5)) {
               
        
                for (size_t i = 1; i + 1 < n; i += 2) {
                    for(size_t j=1;j+1<m;j+=2){
                        T *d = data + begin + i* stride1+j*stride2;
                        //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1), *(d + stride1),*(d - stride2), *(d + stride2)),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_linear(*(d - stride1), *(d + stride1))+coeff_y*interp_linear(*(d - stride2), *(d + stride2)),mode);

                    }
                    if(m%2 ==0){
                        T *d = data + begin + i * stride1+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1), *(d + stride1)),mode);//to determine whether 2d or 1d 
                    }
                }
                if (n % 2 == 0) {
                    for(size_t j=1;j+1<m;j+=2){
                        T *d = data + begin + (n-1) * stride1+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride2), *(d + stride2)),mode);//to determine whether 2d or 1d 
                    }
                    if(m%2 ==0){
                        T *d = data + begin + (n-1) * stride1+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d - stride1-stride2), *(d - stride1), *(d - stride2)),mode);//to determine whether use lorenzo or not
                    }          
                }
                    
            }
                    
            else{//cubic

                if(n<5){//m>=5
                    begin=begin1+begin2+stride1,end=begin+(m-1)*stride2;
                    for(size_t i=1;i<n;i+=2){
                        
                        predict_error+=block_interpolation_1d(data,  begin, end,  stride2,interp_func,pb,meta,tuning);
                        begin+=2*stride1;
                        end+=2*stride1;
                    }
                    return predict_error;
                }
                else if(m<5){//n>=5
                    begin=begin1+begin2+stride2,end=begin+(n-1)*stride1;
                    for(size_t j=1;j<m;j+=2){
                        
                        predict_error+=block_interpolation_1d(data,  begin, end,  stride1,interp_func,pb,meta,tuning);
                        begin+=2*stride2;
                        end+=2*stride2;
                    }
                    return predict_error;

                }
                size_t stride3x1=3*stride1,stride3x2=3*stride2,stride5x1=5*stride1,stride5x2=5*stride2,stride2x1=2*stride1,stride2x2=2*stride2;
                //adaptive todo
              
                   
                size_t i,j;
                T *d;
                if(!meta.adjInterp){
                    for (i = 3; i + 3 < n; i += 2) {
                       
                        for(j=3;j+3<m;j+=2){
                            d = data + begin + i* stride1+j*stride2;


                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                        ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //j=1
                        d = data + begin + i* stride1+stride2;
                        //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                        //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);//to determine
                        //j=m-3 or m-2
                        d = data +begin + i* stride1+j*stride2;
                        //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                        //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        //j=m-1
                        if(m%2 ==0){
                            d = data + begin + i * stride1+(m-1)*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }
                    }
                    //i=1
                    for(j=3;j+3<m;j+=2){
                        d = data + begin + stride1+j*stride2;
                        
                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=1
                    d = data + begin + stride1+stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                            
                    //j=m-3 or m-2
                    d = data +begin + stride1+j*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                    //j=m-1
                    if(m%2 ==0){
                        d = data + begin + stride1+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        d = data + begin + i*stride1+j*stride2;
                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);


                    }
                    //j=1
                    d = data + begin + i*stride1+stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                    
                    //j=m-3 or m-2
                    d = data +begin + i*stride1+j*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                    //j=m-1
                    if(m%2 ==0){
                        d = data + begin + i * stride1+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            d = data + begin + (n-1)*stride1+j*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=1
                        d = data + begin + (n-1)*stride1+stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        //j=m-3 or m-2
                        d = data +begin + (n-1)*stride1+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        //j=m-1
                        if(m%2 ==0){
                            d = data + begin + (n-1) * stride1+(m-1)*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)), interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)) ,mode);
                        } 
                    }
                }
                else{
                    size_t j_start;
                    //first half (non-adj)
                    for (i = 3; i + 3 < n; i += 2) {
                        j_start= (i%4==1)?5:3;
                        for(j=j_start;j+3<m;j+=4){

 
                            d = data + begin + i* stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                        ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d , coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=1
                        if(j_start==5){
                            
                            d = data + begin + i* stride1+stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);//to determine
                        }
                        
                        //j=m-3 or m-2 or j=m-1
                        if(j<m){
                            d = data +begin + i* stride1+j*stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                            //predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                ,mode);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);

                        }
                        
                            /*                          
                        //j=m-1
                        else if(j<m){
                            d = data + begin1 + i * stride1+begin2+j*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }

                        */
                        
                    }
                    //i=1
                    for(j=5;j+3<m;j+=4){
                        d = data + begin + stride1+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=1
                    d = data + begin + stride1+stride2;

                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                    //j=m-3 or m-2
                    if(j<m-1){
                        d = data +begin + stride1+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                    }
                    else if(j<m){//j=m-1

                        d = data + begin + stride1+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    }


                    //i=n-3 or n-2
                    // std::cout<<"f3"<<std::endl;
                    j_start= (i%4==1)?5:3;
                    for(j=j_start;j+3<m;j+=4){
   
                        d = data + begin + i*stride1+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);

                    }
                    
                    //j=1
                    if(j_start==5){
    
                        d = data + begin + i*stride1+stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                    }
                    //j=m-3 or m-2
                    if(j<m-1){
  
                        d = data +begin + i*stride1+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                    }                    
                    //j=m-1
                    else if(j<m){
   
                        d = data + begin + i * stride1+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                    }


                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        j_start= ((n-1)%4==1)?5:3;
                        for(j=j_start;j+3<m;j+=4){
 
                            d = data + begin + (n-1)*stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }

                        //j=1
                        if(j_start==5){
 
                            d = data + begin + (n-1)*stride1+stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
 
                            d = data +begin+ (n-1)*stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        else if(j<m){
 
                            d = data + begin+ (n-1) * stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)), interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)), mode);
                        } 
                    }

                    //second half (adj)
                    for (i = 3; i + 3 < n; i += 2) {
                        j_start= (i%4==1)?3:5;
                        for(j=j_start;j+3<m;j+=4){
                            
                            d = data + begin + i* stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1))
                            //                                        ,interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  );,mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1))
                                                                            +coeff_y*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2)) ,mode);
                        }
                        //j=1
                        if(j_start==5){
                     
                            d = data + begin + i* stride1+stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1)),mode);//to determine
                        }
                        //j=m-3 or m-2 or m-1
                        if(j<m){
         
                            d = data +begin + i* stride1+j*stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1)),mode);//to determine
                        }
                        /*
                        //j=m-1
                        else if(j<m){
                            d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }
                        */

                    }

                    //i=1
                    for(j=3;j+3<m;j+=4){
                      
                        d = data + begin + stride1+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2)) ,mode);
                    }
                    /*
                    //j=1
                    d = data + begin1 + stride1+ begin2+stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                    */    
                    //j=m-3 or m-2
                    if(j<m-1){
                       
                        d = data +begin + stride1+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)), interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    //j=m-1
                    else if(j<m){
                        
                        d = data + begin + stride1+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),mode);//to determine
                    }

                    //i= n-3 or n-2
                    j_start= (i%4==1)?3:5;
                    for(j=j_start;j+3<m;j+=4){
                        
                        d = data + begin + i* stride1+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  ,mode);
                    }
                    //j=1
                    if(j_start==5){
                       
                        d = data + begin + i*stride1+stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                        //                                                            , interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)),mode);
                    }
                    //j=m-3 or m-2
                    if(j<m-1){
                       
                        d = data +begin + i*stride1+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)), interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    //j=m-1
                    else if(j<m){
                       
                        d = data + begin + i * stride1+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)) ,mode);//to determine
                    }
                    
                    //i==n-1
                    if (n % 2 == 0) {
                        j_start= ((n-1)%4==1)?3:5;
                        for(j=j_start;j+3<m;j+=4){
                            
                            d = data + begin + (n-1)* stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  ,mode);
                        }
                        //j=1
                        if(j_start==5){
                            
                            d = data + begin + (n-1)*stride1+stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            d = data +begin + (n-1)*stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        else if(j<m){
                            d = data + begin + (n-1) * stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)),mode);
                        } 
                    }
                }  
            }      

           // quant_index=quant_idx;
            return predict_error;
        }
        double block_interpolation_2d_crossblock(T *data, const std::array<size_t,N> &begin_idx, const std::array<size_t,N> &end_idx,const std::array<size_t,2> &directions,const size_t &math_stride, const std::string &interp_func, const PredictorBehavior pb,const std::array<float,2> &dim_coeffs,const Interp_Meta &meta,int cross_block=1,int tuning=0) {
            size_t direction1=directions[0],direction2=directions[1];
            size_t math_begin_idx1=begin_idx[direction1],math_end_idx1=end_idx[direction1],math_begin_idx2=begin_idx[direction2],math_end_idx2=end_idx[direction2];
            size_t n = (math_end_idx1 - math_begin_idx1) / math_stride + 1, m = (math_end_idx2 - math_begin_idx2) / math_stride + 1;
            bool cross_back=cross_block>0;

            if (n <= 1||m<=1) {
                return 0;
            }
            size_t real_n=cross_back?(math_end_idx1 / math_stride + 1):n,real_m=cross_back?(math_end_idx2 / math_stride + 1):m;
            double predict_error = 0;
            
            float coeff_x=(dim_coeffs[0])/((dim_coeffs[0])+(dim_coeffs[1])),coeff_y=1-coeff_x;

            size_t begin=0;//,global_end_idx1=global_dimensions[direction1],global_end_idx2=global_dimensions[direction2];
            for(size_t i=0;i<N;i++)
                begin+=dimension_offsets[i]*begin_idx[i];

            size_t stride1=math_stride*dimension_offsets[direction1],stride2=math_stride*dimension_offsets[direction2];

            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            //size_t quant_idx=quant_index;
         


            if (interp_func == "linear"||(real_n<5 and real_m<5)) {
               
        
                for (size_t i = 1; i + 1 < n; i += 2) {
                    for(size_t j=1;j+1<m;j+=2){
                        T *d = data + begin + i* stride1+j*stride2;
                        //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1), *(d + stride1),*(d - stride2), *(d + stride2)),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_linear(*(d - stride1), *(d + stride1))+coeff_y*interp_linear(*(d - stride2), *(d + stride2)),mode);

                    }
                    if(m%2 ==0){
                        T *d = data + begin+i * stride1+(m-1)*stride2;
                        
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1), *(d + stride1)),mode);//to determine whether 2d or 1d 
                    }
                }
                if (n % 2 == 0) {
                    
                    
                    for(size_t j=1;j+1<m;j+=2){
                        T *d = data + begin + (n-1) * stride1+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride2), *(d + stride2)),mode);//to determine whether 2d or 1d 
                    }

                    if(m%2 ==0){
                        T *d = data + begin + (n-1) * stride1+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d - stride1-stride2), *(d - stride1), *(d - stride2)),mode);//to determine whether use lorenzo or not
                    }
                            
                }
                    
            }
                    
            else{//cubic

                 if(real_n<5){//real_m>=5
                    std::array<size_t,N> new_begin_idx=begin_idx,new_end_idx=end_idx;
                    for(size_t i=1;i<n;i+=2){
                        new_begin_idx[direction1]=math_begin_idx1+i*math_stride;
                        new_end_idx[direction1]=math_begin_idx1+i*math_stride;
                        predict_error+=block_interpolation_1d_crossblock(data, new_begin_idx,new_end_idx,  direction2,math_stride,interp_func,pb,meta,cross_block,tuning);
                    }
                    return predict_error;
                }
                else if(real_m<5){//real_n>=5
                    std::array<size_t,N> new_begin_idx=begin_idx,new_end_idx=end_idx;
                    for(size_t j=1;j<m;j+=2){
                        new_begin_idx[direction2]=math_begin_idx2+j*math_stride;
                        new_end_idx[direction2]=math_begin_idx2+j*math_stride;
                        predict_error+=block_interpolation_1d_crossblock(data, new_begin_idx,new_end_idx,  direction1,math_stride,interp_func,pb,meta,cross_block,tuning);
                    }
                    return predict_error;

                }
                size_t stride2x1 = 2 * stride1,stride2x2 = 2 * stride2;
                size_t stride3x1 = 3 * stride1,stride3x2 = 3 * stride2;
                size_t math_stride2x=2*math_stride;
                //size_t math_stride3x=3*math_stride;
                //adaptive todo
              
                size_t i,j;
                T *d;
                size_t i_start=(cross_back and math_begin_idx1>=math_stride2x)?1:3;
                size_t j_start=(cross_back and math_begin_idx2>=math_stride2x)?1:3;
                if(!meta.adjInterp){

                    bool i1_b=(i_start==3) and n>4;
                    //bool in_b= (n%2==0) and n>2;
                    bool j1_b=(j_start==3) and m>4;
                    //bool jm_b= (m%2==0) and m>2;

                    for (i = i_start; i + 3 < n; i += 2) {
                       
                        for(j=j_start;j+3<m;j+=2){
                            d = data + begin +i* stride1+j*stride2;


                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                        ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //j=1
                        if(j1_b){
                            d = data + begin+i* stride1+stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);//to determine
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            d = data +begin + i* stride1+j*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                        //j=m-1
                        if(m%2==0){
                            d = data + begin + i * stride1+(m-1)*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }
                    }
                    //i=1
                    if(i1_b){
                        for(j=j_start;j+3<m;j+=2){
                            d = data + begin + stride1+j*stride2;
                            
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //j=1
                        if(j1_b){
                            d = data + begin + stride1+stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                        }   
                        //j=m-3 or m-2
                        if(j<m-1){
                            d = data +begin + stride1+j*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        if(m%2==0){
                            d = data + begin + stride1+(m-1)*stride2;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                    }
                    //i= n-3 or n-2
                    if(i<n-1){
                        for(j=j_start;j+3<m;j+=2){
                            d = data + begin + i*stride1+j*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);


                        }
                        //j=1
                        if(j1_b){
                            d = data + begin + i*stride1+stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            d = data +begin + i*stride1+j*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        if(m%2==0){
                            d = data + begin + i * stride1+(m-1)*stride2;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }
                    //i=n-1 (odd)
                    if (n%2==0) {
                        for(j=j_start;j+3<m;j+=2){
                            d = data + begin + (n-1)*stride1+j*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=1
                        if(j1_b){
                            d = data + begin + (n-1)*stride1+stride2;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            d = data +begin + (n-1)*stride1+j*stride2;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        if(m%2==0){
                            d = data + begin + (n-1) * stride1+(m-1)*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)), interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)) ,mode);
                        } 
                    }
                }
                else{
                    size_t j_start_temp=(j_start==1)?1:5;
                    //first half (non-adj)

                    for (i = i_start; i + 3 < n; i += 2) {
                        j_start= (i%4==1)?j_start_temp:3;
                        for(j=j_start;j+3<m;j+=4){
                            

 
                            d = data + begin + i* stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                        ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d , coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=1

                        if(j_start==5 and m>4){
                            
                            d = data + begin + i* stride1+stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);//to determine
                        }
                        
                        //j=m-3 or m-2 or j=m-1
                        if(j<m){
                            d = data +begin + i* stride1+j*stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                            //predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                ,mode);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);

                        }
                        
                            /*                          
                        //j=m-1
                        else if(j<m){
                            d = data + begin1 + i * stride1+begin2+j*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }

                        */
                        
                    }
                    //i=1
                    if(i_start==3 and n>4){
                        for(j=j_start_temp;j+3<m;j+=4){
                            d = data + begin+ stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //j=1
                        if(j_start_temp==5 and m>4){
                            d = data + begin + stride1+stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            d = data +begin + stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        else if(j<m){//j=m-1

                            d = data + begin + stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                    }

                    //i=n-3 or n-2
                    if(i<n-1){
                        j_start= (i%4==1)?j_start_temp:3;
                        for(j=j_start;j+3<m;j+=4){
       
                            d = data + begin + i*stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);

                        }
                        
                        //j=1
                        if(j_start==5 and m>4){
        
                            d = data + begin + i*stride1+stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
      
                            d = data +begin + i*stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        }                    
                        //j=m-1
                        else if(j<m){
       
                            d = data + begin + i * stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }


                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        j_start= ((n-1)%4==1)?j_start_temp:3;
                        for(j=j_start;j+3<m;j+=4){
 
                            d = data + begin + (n-1)*stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }

                        //j=1
                        if(j_start==5 and m>4){
 
                            d = data + begin + (n-1)*stride1+stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
 
                            d = data +begin + (n-1)*stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        else if(j<m){
 
                            d = data + begin + (n-1) * stride1+j*stride2;

                            //redict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)), interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)), mode);
                        } 
                    }

                    //second half (adj)
                    for (i = i_start; i + 3 < n; i += 2) {
                        j_start= (i%4==1)?3:j_start_temp;
                        for(j=j_start;j+3<m;j+=4){
                            
                            d = data + begin + i* stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1))
                            //                                        ,interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  );,mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1))
                                                                            +coeff_y*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2)) ,mode);
                        }
                        //j=1

                        if(j_start==5 and m>4){
                     
                            d = data + begin + i* stride1+stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1)),mode);//to determine
                        }

                        //j=m-3 or m-2 or m-1
                        if(j<m){
         
                            d = data +begin+ i* stride1+j*stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1)),mode);//to determine
                        }
                        /*
                        //j=m-1
                        else if(j<m){
                            d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }
                        */

                    }

                    //i=1
                    if(i_start==3 and n>4){
                        for(j=3;j+3<m;j+=4){
                          
                            d = data + begin + stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2)) ,mode);
                        }

                        /*
                        //j=1
                        d = data + begin1 + stride1+ begin2+stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                        */    
                        //j=m-3 or m-2
                        if(j<m-1){
                           
                            d = data +begin + stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)), interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                            +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        //j=m-1
                        else if(j<m){
                            d = data + begin + stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),mode);//to determine
                        }
                    }
                    //i= n-3 or n-2
                    if(i<n-1){
                        j_start= (i%4==1)?3:j_start_temp;
                        for(j=j_start;j+3<m;j+=4){
                            
                            d = data + begin + i* stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  ,mode);
                        }
                        //j=1
                        if(j_start==5 and m>4){
                           
                            d = data + begin + i*stride1+stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                            //                                                            , interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)),mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                           
                            d = data +begin + i*stride1+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)), interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        //j=m-1
                        else if(j<m){
                           
                            d = data + begin + i * stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)) ,mode);//to determine
                        }
                    }
                    
                    //i==n-1
                    if (n % 2 == 0) {
                        j_start= ((n-1)%4==1)?3:j_start_temp;
                        for(j=j_start;j+3<m;j+=4){
                            
                            d = data + begin + (n-1)* stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  ,mode);
                        }
                        //j=1
                        if(j_start==5 and m>4){
                            
                            d = data + begin + (n-1)*stride1+stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            d = data +begin + (n-1)*stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        else if(j<m){
                            d = data + begin + (n-1) * stride1+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)),mode);
                        } 
                    }
                }
    
            }    

          //  quant_index=quant_idx;  
            return predict_error;
        }
        double block_interpolation_2d_crossblock_3d(T *data, const std::array<size_t,N> &begin_idx, const std::array<size_t,N> &end_idx,const std::array<size_t,2> &directions, std::array<size_t,N> &steps,const size_t &math_stride, const std::string &interp_func, const PredictorBehavior pb,const std::array<float,2> &dim_coeffs,const Interp_Meta &meta,int cross_block=1,int tuning=0) {
            for(size_t i=0;i<N;i++){
                if(end_idx[i]<begin_idx[i])
                    return 0;
            }
            
            size_t direction1=directions[0],direction2=directions[1];
            size_t math_begin_idx1=begin_idx[direction1],math_end_idx1=end_idx[direction1],math_begin_idx2=begin_idx[direction2],math_end_idx2=end_idx[direction2];
            size_t n = (math_end_idx1 - math_begin_idx1) / math_stride + 1, m = (math_end_idx2 - math_begin_idx2) / math_stride + 1;
            bool cross_back=cross_block>0;

            if (n <= 1||m<=1) {
                return 0;
            }
            size_t real_n=cross_back?(math_end_idx1 / math_stride + 1):n,real_m=cross_back?(math_end_idx2 / math_stride + 1):m;
            
            double predict_error = 0;


            
            float coeff_x=(dim_coeffs[0])/((dim_coeffs[0])+(dim_coeffs[1])),coeff_y=1-coeff_x;

            size_t begin=0;//,global_end_idx1=global_dimensions[direction1],global_end_idx2=global_dimensions[direction2];
            for(size_t i=0;i<N;i++)
                begin+=dimension_offsets[i]*begin_idx[i];
            size_t stride1=math_stride*dimension_offsets[direction1],stride2=math_stride*dimension_offsets[direction2];
            std::array<size_t,N>begins,ends,strides;
            for(size_t i=0;i<N;i++){
                begins[i]=0;
                ends[i]=end_idx[i]-begin_idx[i]+1;
                strides[i]=dimension_offsets[i];
            }
            strides[direction1]=stride1;
            strides[direction2]=stride2;

            int mode=(pb == PB_predict_overwrite)?tuning:-1;
           // size_t quant_idx=quant_index;
           


            if (interp_func == "linear"||(real_n<5 and real_m<5)) {
               
                begins[direction1]=1;
                ends[direction1]=n-1;
                begins[direction2]=1;
                ends[direction2]=m-1;
                steps[direction1]=2;
                steps[direction2]=2;
                for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                    for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                        for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                            T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                            //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1), *(d + stride1),*(d - stride2), *(d + stride2)),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_linear(*(d - stride1), *(d + stride1))+coeff_y*interp_linear(*(d - stride2), *(d + stride2)),mode);
                        }

                    }
                }
                if(m%2 ==0){
                    begins[direction2]=m-1;
                    ends[direction2]=m;
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1), *(d + stride1)),mode);//to determine whether 2d or 1d 
                            }
                        }
                    }
                }
                if (n % 2 == 0) {
                    begins[direction1]=n-1;
                    ends[direction1]=n;
                    begins[direction2]=1;
                    ends[direction2]=m-1;
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                                predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride2), *(d + stride2)),mode);//to determine whether 2d or 1d 
                            }
                        }
                    }

                    if(m%2 ==0){
                        begins[direction2]=m-1;
                        ends[direction2]=m;
                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    //std::cout<<"q4"<<std::endl;
                                    T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                                    predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d - stride1-stride2), *(d - stride1), *(d - stride2)),mode);//to determine whether use lorenzo or not
                                }
                            }
                        }

                    }
                            
                }
                    
            }
                    
            else{//cubic
                 if(real_n<5){//real_m>=5
                    std::array<size_t,N> new_begin_idx=begin_idx,new_end_idx=end_idx,new_steps=steps;
                    new_begin_idx[direction1]=math_begin_idx1+math_stride;
                    new_end_idx[direction1]=math_begin_idx1+(n-1)*math_stride;
                    new_steps[direction1]=2*math_stride;
                    predict_error+=block_interpolation_1d_crossblock_3d(data, new_begin_idx,new_end_idx,  direction2,new_steps,math_stride,interp_func,pb,meta,cross_block,tuning);
                return predict_error;
                }
                else if(real_m<5){//real_n>=5
                    std::array<size_t,N> new_begin_idx=begin_idx,new_end_idx=end_idx,new_steps=steps;
                    new_begin_idx[direction2]=math_begin_idx2+math_stride;
                    new_end_idx[direction2]=math_begin_idx2+(m-1)*math_stride;
                    new_steps[direction2]=2*math_stride;
                    predict_error+=block_interpolation_1d_crossblock_3d(data, new_begin_idx,new_end_idx,  direction1,new_steps,math_stride,interp_func,pb,meta,cross_block,tuning);
                    return predict_error;

                }


                size_t stride2x1 = 2 * stride1,stride2x2 = 2 * stride2;
                size_t stride3x1 = 3 * stride1,stride3x2 = 3 * stride2;
                size_t math_stride2x=2*math_stride;
                //size_t math_stride3x=3*math_stride;
                //adaptive todo
              
                size_t i,j;
                T *d;
                size_t i_start=(cross_back and math_begin_idx1>=math_stride2x)?1:3;
                size_t j_start=(cross_back and math_begin_idx2>=math_stride2x)?1:3;
                if(!meta.adjInterp){
                    bool i1_b=(i_start==3) and n>4;
                    bool j1_b=(j_start==3) and m>4;
                    steps[direction1]=2;
                    steps[direction2]=2;
                
                    //i=1
                    if(i1_b){

                        begins[direction1]=1;
                        ends[direction1]=2;
                        //j=1
                        if(j1_b){
                            begins[direction2]=1;
                            ends[direction2]=2;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];
                          
                                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                                    }

                                }
                            }
                        }
                        begins[direction2]=j_start;
                        ends[direction2]=(m>=3)?(m-3):0;

                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];                  
                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                                }
                            }
                        }
                           
                        //j=m-3 or m-2
                        if(m>2){
                            begins[direction2]=m+m%2-3;
                            ends[direction2]=begins[direction2]+1;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                        //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                                    }
                                }
                            }
                        }

                        //j=m-1
                        if(m%2==0){
                            begins[direction2]=m-1;
                            ends[direction2]=m;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                                    }
                                }
                            }
                        }
                    }
                    begins[direction1]=i_start;
                    ends[direction1]=(n>=3)?(n-3):0;

                    
                        //j=1
                    if(j1_b){
                        begins[direction2]=1;
                        ends[direction2]=2;
                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];               
                                    //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);//to determine
                                }
                            }
                        }
                    }
                    begins[direction2]=j_start;
                    ends[direction2]=(m>=3)?(m-3):0;
    
                    for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                        for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                            for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                

                                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                //                                        ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                    }

                    if(m>2){
                        begins[direction2]=m+m%2-3;
                        ends[direction2]=begins[direction2]+1;
                         for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                             for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  

                                    //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                                }
                            }
                        }
                    }
                    //j=m-1
                    if(m%2==0){
                        begins[direction2]=m-1;
                        ends[direction2]=m;
                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                    //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                 ,mode);
                                }
                            }
                        }
                    }

                    //i= n-3 or n-2
                    if(n>2){
                        begins[direction1]=n+n%2-3;
                        ends[direction1]=begins[direction1]+1;
                        //j=1
                        if(j1_b){
                            begins[direction2]=1;
                            ends[direction2]=2;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                                    }
                                }
                            }
                        }
                        begins[direction2]=j_start;
                        ends[direction2]=(m>=3)?(m-3):0;
                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                                }

                            }
                        }
                        
                        //j=m-3 or m-2
                        if(m>2){
                            begins[direction2]=m+m%2-3;
                            ends[direction2]=begins[direction2]+1;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                                    }
                                }
                            }
                        }
                        //j=m-1
                        if(m%2==0){
                            begins[direction2]=m-1;
                            ends[direction2]=m;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                                    }
                                }
                            }
                        }
                    }
                    //i=n-1 (odd)
                    if (n%2==0) {
                        begins[direction1]=n-1;
                        ends[direction1]=n;
                        //j=1
                        if(j1_b){
                            begins[direction2]=1;
                            ends[direction2]=2;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                                    }
                                }
                            }
                        }
                        begins[direction2]=j_start;
                        ends[direction2]=(m>=3)?(m-3):0;
                        for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                            for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                    T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                                }
                            }
                        }
                        
                        //j=m-3 or m-2
                        if(m>2){
                            begins[direction2]=m+m%2-3;
                            ends[direction2]=begins[direction2]+1;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                        predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                                    }
                                }
                            }
                        }
                        //j=m-1
                        if(m%2==0){
                            begins[direction2]=m-1;
                            ends[direction2]=m;
                            for(size_t i=begins[0];i<ends[0];i+=steps[0]){
                                for(size_t j=begins[1];j<ends[1];j+=steps[1]){
                                    for(size_t k=begins[2];k<ends[2];k+=steps[2]){
                                        T *d = data + begin + i * strides[0]+j*strides[1]+k*strides[2];  
                                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)), interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),mode);
                                        predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)) ,mode);
                                    }
                                }
                            }
                        } 
                    }
                }
                else{
                    if(direction1!=2 and direction2!=2){//temp. Too hard to generalize....
                        size_t j_start_temp=(j_start==1)?1:5;
                        size_t k_start=begins[2],k_end=ends[2],k_step=steps[2],k_stride=strides[2];
                        //first half (non-adj)
                        for (i = i_start; i + 3 < n; i += 2) {
                            j_start= (i%4==1)?j_start_temp:3;
                            for(j=j_start;j+3<m;j+=4){
                                for(size_t k=k_start;k<k_end;k+=k_step){

     
                                    d = data + begin + i* stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                                        ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                                    predict_error+=quantize_integrated(d - data, *d , coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                                }
                            }
                            //j=1
                            if(j_start==5 and m>4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                
                                    d = data + begin + i* stride1+stride2+k*k_stride;

                                    //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);//to determine
                                }
                            }
                            
                            //j=m-3 or m-2 or j=m-1
                            if(j<m){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data +begin + i* stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                                    //predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                                ,mode);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                                }

                            }
                            
                                /*                          
                            //j=m-1
                            else if(j<m){
                                d = data + begin1 + i * stride1+begin2+j*stride2;
                                //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                             ,mode);
                            }

                            */
                            
                        }
                        //i=1
                        if(i_start==3 and n>4){
                            for(j=j_start_temp;j+3<m;j+=4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data + begin+ stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                                }
                            }
                            //j=1
                            if(j_start_temp==5 and m>4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data + begin + stride1+stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                                }                           
                            }
                            //j=m-3 or m-2
                            if(j<m-1){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data +begin + stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                                }
                            }
                            else if(j<m){//j=m-1
                                for(size_t k=k_start;k<k_end;k+=k_step){

                                    d = data + begin + stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                                }
                            }
                        }

                        //i=n-3 or n-2
                        if(i<n-1){
                            j_start= (i%4==1)?j_start_temp:3;
                            for(j=j_start;j+3<m;j+=4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
           
                                    d = data + begin + i*stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                                }

                            }
                            
                            //j=1
                            if(j_start==5 and m>4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
            
                                    d = data + begin + i*stride1+stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                                }
                            }
                            //j=m-3 or m-2
                            if(j<m-1){
                                for(size_t k=k_start;k<k_end;k+=k_step){
          
                                    d = data +begin + i*stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                                }
                            }                    
                            //j=m-1
                            else if(j<m){
                                for(size_t k=k_start;k<k_end;k+=k_step){
           
                                    d = data + begin + i * stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                                }
                            }
                        }


                        //i=n-1 (odd)
                        if (n % 2 == 0) {
                            j_start= ((n-1)%4==1)?j_start_temp:3;
                            for(j=j_start;j+3<m;j+=4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
     
                                    d = data + begin + (n-1)*stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)) ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                                }
                            }

                            //j=1
                            if(j_start==5 and m>4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
     
                                    d = data + begin + (n-1)*stride1+stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                                }
                            }
                            //j=m-3 or m-2
                            if(j<m-1){
                                for(size_t k=k_start;k<k_end;k+=k_step){
     
                                    d = data +begin + (n-1)*stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                                }
                            }
                            //j=m-1
                            else if(j<m){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data + begin + (n-1) * stride1+j*stride2+k*k_stride;

                                    //redict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)), interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),mode);
                                    predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)), mode);
                                }
                            } 
                        }

                        //second half (adj)
                        for (i = i_start; i + 3 < n; i += 2) {
                            j_start= (i%4==1)?3:j_start_temp;
                            for(j=j_start;j+3<m;j+=4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                
                                    d = data + begin + i* stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1))
                                    //                                        ,interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  );,mode);
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2)) ,mode);
                                }
                            }
                            //j=1
                            if(j_start==5 and m>4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                         
                                    d = data + begin + i* stride1+stride2+k*k_stride;

                                    //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1)),mode);//to determine
                                }
                            }
                            //j=m-3 or m-2 or m-1
                            if(j<m){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data +begin+ i* stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                    //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1)),mode);//to determine
                                }
                            }
                            /*
                            //j=m-1
                            else if(j<m){
                                d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                                //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                             ,mode);
                            }
                            */

                        }

                        //i=1
                        if(i_start==3 and n>4){
                            for(j=3;j+3<m;j+=4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                              
                                    d = data + begin + stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2)) ,mode);
                                }
                            }

                            /*
                            //j=1
                            d = data + begin1 + stride1+ begin2+stride2;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                            */    
                            //j=m-3 or m-2
                            if(j<m-1){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                               
                                    d = data +begin + stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)), interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                                }
                            }
                            //j=m-1
                            else if(j<m){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data + begin + stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),mode);//to determine
                                }
                            }
                        }
                        //i= n-3 or n-2
                        if(i<n-1){
                            j_start= (i%4==1)?3:j_start_temp;
                            for(j=j_start;j+3<m;j+=4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                
                                    d = data + begin + i* stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  ,mode);
                                }
                            }
                            //j=1
                            if(j_start==5 and m>4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                               
                                    d = data + begin + i*stride1+stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                    //                                                            , interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ),mode);
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)),mode);
                                }
                            }
                            //j=m-3 or m-2
                            if(j<m-1){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data +begin + i*stride1+j*stride2+k*k_stride;

                                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)), interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                    +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                                }
                            }
                            //j=m-1
                            else if(j<m){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data + begin + i * stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)) ,mode);//to determine
                                }
                            }
                        }
                        
                        //i=n-1               
                        if (n % 2 == 0) {
                            j_start= ((n-1)%4==1)?3:j_start_temp;
                            for(j=j_start;j+3<m;j+=4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                
                                    d = data + begin + (n-1)* stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  ,mode);
                                }
                            }
                            //j=1
                            if(j_start==5 and m>4){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                
                                    d = data + begin + (n-1)*stride1+stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                                }
                            }
                            //j=m-3 or m-2
                            if(j<m-1){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data +begin + (n-1)*stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ,mode);
                                }
                            }
                            //j=m-1
                            else if(j<m){
                                for(size_t k=k_start;k<k_end;k+=k_step){
                                    d = data + begin + (n-1) * stride1+j*stride2+k*k_stride;

                                    predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)),mode);
                                }
                            } 
                        }

                    }
                    else{
                        size_t sub_direction=3-direction1-direction2;
                        size_t sub_start=begin_idx[sub_direction],sub_end=end_idx[sub_direction],sub_step=steps[sub_direction];
                        std::array<size_t,N>temp_start=begin_idx,temp_end=end_idx;
                        for(size_t sub=sub_start;sub<=sub_end;sub+=sub_step){
                          
                            temp_start[sub_direction]=temp_end[sub_direction]=sub;
                            predict_error+=block_interpolation_2d_crossblock(data, temp_start, temp_end,directions,math_stride, interp_func, pb,dim_coeffs,meta,cross_block,tuning);
                        }
                        return predict_error;



                    }
                }
    
            }      
           //  quant_index=quant_idx;

           
            return predict_error;
        }
        double block_interpolation_3d(T *data, size_t begin1, size_t end1, size_t begin2, size_t end2, size_t begin3, size_t end3, size_t stride1,size_t stride2,size_t stride3, const std::string &interp_func, const PredictorBehavior pb,const std::array<float,3> &dim_coeffs,const Interp_Meta &meta,int tuning=0) {
            size_t n = (end1 - begin1) / stride1 + 1;
            if (n <= 1) {
                return 0;
            }
            size_t m = (end2 - begin2) / stride2 + 1;
            if (m <= 1) {
                return 0;
            }
            size_t p = (end3 - begin3) / stride3 + 1;
            if (p <= 1) {
                return 0;
            }
            size_t begin=begin1+begin2+begin3,end=end1+end2+end3;
            double predict_error = 0;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;

            float coeff_x=dim_coeffs[0]/(dim_coeffs[0]+dim_coeffs[1]+dim_coeffs[2]);
            float coeff_y=dim_coeffs[1]/(dim_coeffs[0]+dim_coeffs[1]+dim_coeffs[2]);
            float coeff_z=1-coeff_x-coeff_y;

            float coeff_x_xy=(coeff_x)/(coeff_x+coeff_y),coeff_y_xy=1-coeff_x_xy;
            float coeff_x_xz=(coeff_x)/(coeff_x+coeff_z),coeff_z_xz=1-coeff_x_xz;
            float coeff_y_yz=(coeff_y)/(coeff_y+coeff_z),coeff_z_yz=1-coeff_y_yz;
           // size_t quant_idx=quant_index;
         


            if (interp_func == "linear" || (n<5 and m<5 and p<5) ){//nmpcond temp added
                
                for (size_t i = 1; i + 1 < n; i += 2) {
                    for(size_t j=1;j+1<m;j+=2){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin + i* stride1+j*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_linear(*(d - stride1), *(d + stride1))
                                                                            +coeff_y*interp_linear(*(d - stride2), *(d + stride2))
                                                                            +coeff_z*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin + i* stride1+j*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_linear(*(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_linear(*(d - stride2), *(d + stride2)),mode);
                        }

                    }
                    if(m%2 ==0){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin + i* stride1+(m-1)*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_linear(*(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin + i* stride1+(m-1)*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1), *(d + stride1)),mode);
                        }
                    }      
                }
                if (n % 2 == 0) {
                    for(size_t j=1;j+1<m;j+=2){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin + (n-1)* stride1+j*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_linear(*(d - stride2), *(d + stride2))
                                                                                +coeff_z_yz*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin + (n-1)* stride1+j*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride2), *(d + stride2)),mode);
                        }
                    }
                    if(m%2 ==0){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin + (n-1)* stride1+(m-1)*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin + (n-1)* stride1+(m-1)*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                        }
                    }           
                }
            }
            else{//cubic

                if(n<5){
                    if(m<5){//p>=5
                        begin=begin1+begin2+begin3,end=begin+(p-1)*stride3;
                        for(size_t i=1;i<n;i+=2){
                            for(size_t j=1;j<m;j+=2){
                                predict_error+=block_interpolation_1d(data,  begin+i*stride1+j*stride2, end+i*stride1+j*stride2,  stride3,interp_func,pb,meta,tuning);
                                
                            }
                        
                        }
                        return predict_error;
                    }
                    else if(p<5){//m>=5
                        begin=begin1+begin2+begin3,end=begin+(m-1)*stride2;
                        for(size_t i=1;i<n;i+=2){
                            for(size_t k=1;k<p;k+=2){
                            
                                predict_error+=block_interpolation_1d(data,  begin+i*stride1+k*stride3, end+i*stride1+k*stride3,  stride2,interp_func,pb,meta,tuning);
                                
                            }
                       
                        }
                        return predict_error;

                    }
                    else{//mp>=5
                        begin2+=begin1+stride1,end2+=begin1+stride1;
                        for(size_t i=1;i<n;i+=2){
                            predict_error+=block_interpolation_2d(data,  begin2, end2,begin3,end3,  stride2,stride3,interp_func,pb,std::array<float,2>{coeff_y_yz,coeff_z_yz},meta,tuning);
                            begin2+=2*stride1;
                            end2+=2*stride1;
                        }
                        return predict_error;
                    }
                    
                }
                else if(m<5){//n>=5

                    if(p<5){
                         begin=begin1+begin2+begin3,end=begin+(n-1)*stride1;
                        for(size_t j=1;j<m;j+=2){
                            for(size_t k=1;k<p;k+=2){
                            
                                predict_error+=block_interpolation_1d(data,  begin+j*stride2+k*stride3, end+j*stride2+k*stride3,  stride1,interp_func,pb,meta,tuning);
                                
                            }
                        
                        }
                        return predict_error;

                    }
                    else{//np>=5
                        begin1+=begin2+stride2,end1+=begin2+stride2;
                        for(size_t j=1;j<m;j+=2){
                            predict_error+=block_interpolation_2d(data,  begin1, end1,begin3,end3,  stride1,stride3,interp_func,pb,std::array<float,2>{coeff_x_xz,coeff_z_xz},meta,tuning);
                            begin1+=2*stride2;
                            end1+=2*stride2;
                        }
                        return predict_error;

                    }
                    

                }
                else if(p<5){//mn>=5
                    begin2+=begin3+stride3,end2+=begin3+stride3;
                    for(size_t k=1;k<p;k+=2){
                        predict_error+=block_interpolation_2d(data,  begin1, end1,begin2,end2,  stride1,stride2,interp_func,pb,std::array<float,2>{coeff_x_xy,coeff_x_xy},meta,tuning);
                        begin2+=2*stride3;
                        end2+=2*stride3;
                    }
                    return predict_error;

                }

                size_t stride3x1=3*stride1,stride3x2=3*stride2,stride5x1=5*stride1,stride5x2=5*stride2,stride3x3=3*stride3,stride5x3=5*stride3,stride2x1=2*stride1,stride2x2=2*stride2,stride2x3=2*stride3;
                //adaptive todo
              
                   
                size_t i,j,k;
                T *d;
                if(!meta.adjInterp){
                    for (i = 3; i + 3 < n; i += 2) {
                        for(j=3;j+3<m;j+=2){
                            for(k=3;k+3<p;k+=2){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=1
                            d = data + begin + i* stride1+j*stride2+stride3;
                            /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);//should or how we ave for different interps?
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        
                            //k=p-3 or p-2
                            d = data + begin + i* stride1+j*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + i* stride1+j*stride2+(p-1)*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                                */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                                 
                            }
                        }
                        //j=1
                        for(k=3;k+3<p;k+=2){
                            d = data + begin + i* stride1+stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                         }
                        //k=1
                        d = data + begin + i* stride1+stride2+stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);

                        //k=p-3 or p-2
                        d = data + begin + i* stride1+stride2+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin + i* stride1+stride2+(p-1)*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                        //j=m-3 or m-2
                        for(k=3;k+3<p;k+=2){
                            d = data + begin + i* stride1+j*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                         }
                        //k=1
                        d = data + begin + i* stride1+j*stride2+stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        //k=p-3 or p-2
                        d = data + begin + i* stride1+j*stride2+k*stride3;
                        /*
                        predict_error+=quantize_tuning(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin + i* stride1+j*stride2+(p-1)*stride3;
                            /*
                            predict_error+=quantize_tuning(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            for(k=3;k+3<p;k+=2){
                                d = data + begin + i* stride1+(m-1)*stride2+k*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                    interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                                */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1

                            d = data + begin + i* stride1+(m-1)*stride2+stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                    interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            //k=p-3 or p-2
                            d = data + begin + i* stride1+(m-1)*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                    interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            //k=p-1
                            if(p%2==0){
                                d = data + begin+ i* stride1+(m-1)*stride2+(p-1)*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                    interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                                */
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                        }
                    }
                    //i=1
                    for(j=3;j+3<m;j+=2){
                        for(k=3;k+3<p;k+=2){
                            d = data + begin +  stride1+j*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin +  stride1+j*stride2+stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        //k=p-3 or p-2
                        d = data + begin + stride1+j*stride2+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin +  stride1+j*stride2+(p-1)*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    for(k=3;k+3<p;k+=2){
                        d = data + begin +  stride1+stride2+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                            interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                            interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    d = data + begin + stride1+stride2+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                        +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                    //k=p-3 or p-2
                    d = data + begin +  stride1+stride2+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                        +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    //k=p-1
                    if(p%2==0){
                        d = data + begin +  stride1+stride2+(p-1)*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=m-3 or m-2
                    for(k=3;k+3<p;k+=2){
                        d = data + begin +  stride1+j*stride2+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), 
                            interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),
                            interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    d = data + begin + stride1+j*stride2+stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                    +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    //k=p-3 or p-2
                    d = data + begin +  stride1+j*stride2+k*stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                    +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                    //k=p-1
                    if(p%2==0){
                        d = data + begin + stride1+j*stride2+(p-1)*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        for(k=3;k+3<p;k+=2){
                            d = data + begin +  stride1+(m-1)*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin +  stride1+(m-1)*stride2+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        //k=p-3 or p-2
                        d = data + begin + stride1+(m-1)*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin + stride1+(m-1)*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        for(k=3;k+3<p;k+=2){
                            d = data + begin + i* stride1+j*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=1
                        d = data + begin +  i*stride1+j*stride2+stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        //k=p-3 or p-2
                        d = data + begin +i* stride1+j*stride2+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin +  i*stride1+j*stride2+(p-1)*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    for(k=3;k+3<p;k+=2){
                        d = data + begin +  i*stride1+stride2+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                            interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                            interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    d = data + begin + i*stride1+stride2+stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                    +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                    +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    //k=p-3 or p-2
                    d = data + begin +  i*stride1+stride2+k*stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                    +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                    +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    //k=p-1
                    if(p%2==0){
                        d = data + begin +  i*stride1+stride2+(p-1)*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=m-3 or m-2
                    for(k=3;k+3<p;k+=2){
                        d = data + begin +  i*stride1+j*stride2+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), 
                            interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),
                            interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    d = data + begin + i*stride1+j*stride2+stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                    +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                    +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                    //k=p-3 or p-2
                    d = data + begin +  i*stride1+j*stride2+k*stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                    +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                    +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    //k=p-1
                    if(p%2==0){
                        d = data + begin + i*stride1+j*stride2+(p-1)*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        for(k=3;k+3<p;k+=2){
                            d = data + begin +  i*stride1+(m-1)*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), 
                                interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin +  i*stride1+(m-1)*stride2+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        //k=p-3 or p-2
                        d = data + begin + i*stride1+(m-1)*stride2+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin + i*stride1+(m-1)*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            for(k=3;k+3<p;k+=2){
                                d = data + begin + (n-1)* stride1+j*stride2+k*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                                */
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            d = data + begin +  (n-1)*stride1+j*stride2+stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            //k=p-3 or p-2
                            d = data + begin +(n-1)* stride1+j*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            //k=p-1
                            if(p%2==0){
                                d = data + begin +  (n-1)*stride1+j*stride2+(p-1)*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                    interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                                */
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=1
                        for(k=3;k+3<p;k+=2){
                            d = data + begin + (n-1)*stride1+stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin + (n-1)*stride1+stride2+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                        +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        //k=p-3 or p-2
                        d = data + begin +  (n-1)*stride1+stride2+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                        +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin + (n-1)*stride1+stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //j=m-3 or m-2
                        for(k=3;k+3<p;k+=2){
                            d = data + begin + (n-1)*stride1+j*stride2+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin + (n-1)*stride1+j*stride2+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))
                                                                        +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        //k=p-3 or p-2
                        d = data + begin +  (n-1)*stride1+j*stride2+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                        +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin + (n-1)*stride1+j*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            for(k=3;k+3<p;k+=2){
                                d = data + begin +  (n-1)*stride1+(m-1)*stride2+k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            d = data + begin +  (n-1)*stride1+(m-1)*stride2+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);

                            //k=p-3 or p-2
                            d = data + begin + (n-1)*stride1+(m-1)*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + (n-1)*stride1+(m-1)*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }
                }
                else{

                    size_t k_start;
                    //first half (non-adj) 


                    size_t ks1=3,ks2=5;
                 
                    for (i = 3; i + 3 < n; i += 2) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i* stride1 +j*stride2 +stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                            }
                        
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                            }

                        }
                        //j=1
                        k_start=(i+1)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin + i* stride1  +stride2 +k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin + i* stride1 +stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }

                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin + i* stride1 +stride2 +k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }

                        //j=m-3 or m-2 or m-1
                        while(j<m){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1  +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin + i* stride1 +j*stride2 +stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                            j+=2;
                        }
                    }
                    //i=1
                    
                    for(j=3;j+3<m;j+=2){
                        k_start=(1+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin +  stride1 +j*stride2 +k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin +  stride1 +j*stride2 +stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if (k<p){
                            d = data + begin + stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
    
                    }
                    //j=1 (i+j=2)
                    k_start=ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin +  stride1+  +stride2 +k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                            interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                            interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin + stride1 +stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                            +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin +  stride1 +stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                            +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                    }

                    //k=p-1
                    else if(k<p){
                        d = data + begin +  stride1 +stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                    }
                    //j=m-3 or m-2
                    k_start=(1+j)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin +  stride1+  +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin + stride1 +j*stride2 +stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin +  stride1 +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin + stride1 +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    //j=m-1 (i=1)
                    if(m%2 ==0){
                        k_start=(m%4==0)?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin +  stride1  +(m-1)*stride2 +k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin +  stride1 +(m-1)*stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin + stride1 +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin + stride1 +(m-1)*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        k_start=(i+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin + i* stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin +  i*stride1 +j*stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin +i* stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    k_start=(i+1)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin +  i*stride1  +stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin + i*stride1 +stride2 +stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                        +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin +  i*stride1 +stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                        +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin +  i*stride1 +stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=m-3 or m-2
                    k_start=(i+j)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin +  i*stride1  +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin + i*stride1 +j*stride2 +stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                        +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin +  i*stride1 +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin + i*stride1 +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        k_start=(i+m-1)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin +  i*stride1 +(m-1)*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin +  i*stride1 +(m-1)*stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin + i*stride1 +(m-1)*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin + i*stride1 +(m-1)*stride2 +(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(n-1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + (n-1)* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  (n-1)*stride1 +j*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin +(n-1)* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            
                        }
                        //j=1
                        k_start=(n%4==0)?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin + (n-1)*stride1  +stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin + (n-1)*stride1 +stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                            +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin +  (n-1)*stride1 +stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                            +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin + (n-1)*stride1 +stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        k_start=(n-1+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin + (n-1)*stride1  +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin + (n-1)*stride1 +j*stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                            +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin +  (n-1)*stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                            +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin + (n-1)*stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(n+m-2)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  (n-1)*stride1  +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  (n-1)*stride1 +(m-1)*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + (n-1)*stride1 +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + (n-1)*stride1 +(m-1)*stride2 +(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }


                    //}
                    
                    //second half (adj)
                    ks1=5;
                    ks2=3;

                    for (i = 3; i + 3 < n; i += 2) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i* stride1 +j*stride2 +stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ,mode);
                            }
                        
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ,mode);                            }

                        }
                        //j=1
                        k_start=(i+1)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin + i* stride1+  stride2 +k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ,mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin + i* stride1 +stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                        }

                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin + i* stride1 +stride2 +k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                        }

                        //j=m-3 or m-2 or m-1
                        while(j<m){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+  j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin + i* stride1 +j*stride2 +stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                            j+=2;
                        }
                    }
                    //i=1
                    for(j=3;j+3<m;j+=2){
                        k_start=(1+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin +  stride1 +j*stride2 +k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin +  stride1 +j*stride2 +stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if (k<p){
                            d = data + begin + stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
    
                    }
                    //j=1 (i+j=2)
                    k_start=ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin +  stride1+  +stride2 +k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                            interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                            interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin + stride1 +stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                            +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                            +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin +  stride1 +stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                            +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                            +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                    }

                    //k=p-1
                    else if(k<p){
                        d = data + begin +  stride1 +stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y_xy*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                    }
                    //j=m-3 or m-2
                    k_start=(1+j)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin +  stride1+  +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin + stride1 +j*stride2 +stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin +  stride1 +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin + stride1 +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y_xy*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    //j=m-1 (i=1)
                    if(m%2 ==0){
                        k_start=(m%4==0)?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin +  stride1  +(m-1)*stride2 +k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin +  stride1 +(m-1)*stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                            +coeff_z_xz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin + stride1 +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                +coeff_z_xz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin + stride1 +(m-1)*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),mode);
                        }
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        k_start=(i+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin + i* stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin +  i*stride1 +j*stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin +i* stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    k_start=(i+1)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin +  i*stride1+  +stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin + i*stride1 +stride2 +stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                        +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin +  i*stride1 +stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                        +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin +  i*stride1 +stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)),mode);
                    }
                    //j=m-3 or m-2
                    k_start=(i+j)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin +  i*stride1+ j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin + i*stride1 +j*stride2 +stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                        +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin +  i*stride1 +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin + i*stride1 +j*stride2 +k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        k_start=(i+m-1)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin +  i*stride1+ (m-1)*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin +  i*stride1 +(m-1)*stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_z_xz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin + i*stride1 +(m-1)*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_z_xz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin + i*stride1 +(m-1)*stride2 +(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(n-1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + (n-1)* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) 
                                                                                +coeff_z_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  (n-1)*stride1 +j*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin +(n-1)* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            
                        }
                        //j=1
                        k_start=(n%4==0)?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin + (n-1)*stride1  +stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin + (n-1)*stride1 +stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                            +coeff_z_yz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin +  (n-1)*stride1 +stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                            +coeff_z_yz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin + (n-1)*stride1 +stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        k_start=(n-1+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin + (n-1)*stride1+  +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin + (n-1)*stride1 +j*stride2 +stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                            +coeff_z_yz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin +  (n-1)*stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                            +coeff_z_yz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin + (n-1)*stride1 +j*stride2 +k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(n+m-2)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  (n-1)*stride1  +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  (n-1)*stride1 +(m-1)*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + (n-1)*stride1 +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + (n-1)*stride1 +(m-1)*stride2 +(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }

                }
            }

          //  quant_index=quant_idx;
            return predict_error;
        }
        double block_interpolation_3d_crossblock(T *data, const std::array<size_t,N> &begin_idx, const std::array<size_t,N> &end_idx,const std::array<size_t,3> &directions,const size_t &math_stride, const std::string &interp_func, const PredictorBehavior pb,const std::array<float,3> &dim_coeffs,const Interp_Meta &meta,int cross_block=1,int tuning=0) {
            size_t direction1=directions[0],direction2=directions[1],direction3=directions[2];
            size_t math_begin_idx1=begin_idx[direction1],math_end_idx1=end_idx[direction1],math_begin_idx2=begin_idx[direction2],math_end_idx2=end_idx[direction2],math_begin_idx3=begin_idx[direction3],math_end_idx3=end_idx[direction3];
            size_t n = (math_end_idx1 - math_begin_idx1) / math_stride + 1, m = (math_end_idx2 - math_begin_idx2) / math_stride + 1, p = (math_end_idx3 - math_begin_idx3) / math_stride + 1;

            bool cross_back=cross_block>0;

            if (n <= 1||m<=1||p<=1) {
                return 0;
            }
            size_t real_n=cross_back?(math_end_idx1 / math_stride + 1):n,real_m=cross_back?(math_end_idx2 / math_stride + 1):m,real_p=cross_back?(math_end_idx3 / math_stride + 1):p;
            double predict_error = 0;
            
            float coeff_x=dim_coeffs[0]/(dim_coeffs[0]+dim_coeffs[1]+dim_coeffs[2]);
            float coeff_y=dim_coeffs[1]/(dim_coeffs[0]+dim_coeffs[1]+dim_coeffs[2]);
            float coeff_z=1-coeff_x-coeff_y;

            float coeff_x_xy=(coeff_x)/(coeff_x+coeff_y),coeff_y_xy=1-coeff_x_xy;
            float coeff_x_xz=(coeff_x)/(coeff_x+coeff_z),coeff_z_xz=1-coeff_x_xz;
            float coeff_y_yz=(coeff_y)/(coeff_y+coeff_z),coeff_z_yz=1-coeff_y_yz;

            size_t begin=0;//,global_end_idx1=global_dimensions[direction1],global_end_idx2=global_dimensions[direction2],global_end_idx3=global_dimensions[direction3];
            for(size_t i=0;i<N;i++)
                begin+=dimension_offsets[i]*begin_idx[i];

            size_t stride1=math_stride*dimension_offsets[direction1],stride2=math_stride*dimension_offsets[direction2],stride3=math_stride*dimension_offsets[direction3];
            
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
           // size_t quant_idx=quant_index;
            
            if (interp_func == "linear" || (real_n<5 and real_m<5 and real_p<5) ){//nmpcond temp added
                
                for (size_t i = 1; i + 1 < n; i += 2) {
                    for(size_t j=1;j+1<m;j+=2){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin + i* stride1+j*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_linear(*(d - stride1), *(d + stride1))
                                                                            +coeff_y*interp_linear(*(d - stride2), *(d + stride2))
                                                                            +coeff_z*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin + i* stride1+j*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_linear(*(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_linear(*(d - stride2), *(d + stride2)),mode);
                        }

                    }
                    if(m%2 ==0){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin + i* stride1+(m-1)*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_linear(*(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin + i* stride1+(m-1)*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1), *(d + stride1)),mode);
                        }
                    }      
                }
                if (n % 2 == 0) {
                    for(size_t j=1;j+1<m;j+=2){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin + (n-1)* stride1+j*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_linear(*(d - stride2), *(d + stride2))
                                                                                +coeff_z_yz*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin + (n-1)* stride1+j*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride2), *(d + stride2)),mode);
                        }
                    }
                    if(m%2 ==0){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin + (n-1)* stride1+(m-1)*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin + (n-1)* stride1+(m-1)*stride2+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                        }
                    }           
                }
            }
            else{//cubic
                if(real_n<5){
                    if(real_m<5){//real_p>=5
                        std::array<size_t,N>new_begin_idx=begin_idx,new_end_idx=end_idx;
                        for(size_t i=1;i<n;i+=2){
                            for(size_t j=1;j<m;j+=2){
                                new_end_idx[direction1]=new_begin_idx[direction1]=math_begin_idx1+i*math_stride;
                                new_end_idx[direction2]=new_begin_idx[direction2]=math_begin_idx2+j*math_stride;
                                predict_error+=block_interpolation_1d_crossblock(data,  new_begin_idx, new_end_idx, direction3,math_stride,interp_func,pb,meta,cross_block,tuning);
                                
                            }
                        
                        }
                        return predict_error;
                    }
                    else if(real_p<5){//m>=5
                        std::array<size_t,N>new_begin_idx=begin_idx,new_end_idx=end_idx;
                        for(size_t i=1;i<n;i+=2){
                            for(size_t k=1;k<p;k+=2){
                                new_end_idx[direction1]=new_begin_idx[direction1]=math_begin_idx1+i*math_stride;
                                new_end_idx[direction3]=new_begin_idx[direction3]=math_begin_idx3+k*math_stride;
                                predict_error+=block_interpolation_1d_crossblock(data,  new_begin_idx, new_end_idx, direction2,math_stride,interp_func,pb,meta,cross_block,tuning);
                                
                            }
                        
                        }
                        return predict_error;

                    }
                    else{//mp>=5
                        std::array<size_t,N>new_begin_idx=begin_idx,new_end_idx=end_idx;
                        for(size_t i=1;i<n;i+=2){
                            new_end_idx[direction1]=new_begin_idx[direction1]=math_begin_idx1+i*math_stride;
                            predict_error+=block_interpolation_2d_crossblock(data,  new_begin_idx, new_end_idx,std::array<size_t,2>{direction2,direction3}
                                                                            ,math_stride,interp_func,pb,std::array<float,2>{coeff_y_yz,coeff_z_yz},meta,cross_block,tuning);
                        }
                        return predict_error;
                    }
                    
                }
                else if(real_m<5){//real_n>=5

                    if(real_p<5){
                        std::array<size_t,N>new_begin_idx=begin_idx,new_end_idx=end_idx;
                        for(size_t j=1;j<m;j+=2){
                            for(size_t k=1;k<p;k+=2){
                                new_end_idx[direction2]=new_begin_idx[direction2]=math_begin_idx2+j*math_stride;
                                new_end_idx[direction3]=new_begin_idx[direction3]=math_begin_idx3+k*math_stride;
                                predict_error+=block_interpolation_1d_crossblock(data,  new_begin_idx, new_end_idx, direction1,math_stride,interp_func,pb,meta,cross_block,tuning);
                                
                            }
                        
                        }
                        return predict_error;

                    }
                    else{//np>=5
                        std::array<size_t,N>new_begin_idx=begin_idx,new_end_idx=end_idx;
                        for(size_t j=1;j<m;j+=2){
                            new_end_idx[direction2]=new_begin_idx[direction2]=math_begin_idx2+j*math_stride;
                            predict_error+=block_interpolation_2d_crossblock(data,  new_begin_idx, new_end_idx,std::array<size_t,2>{direction1,direction3}
                                                                            ,math_stride,interp_func,pb,std::array<float,2>{coeff_x_xz,coeff_z_xz},meta,cross_block,tuning);
                        }
                        return predict_error;
                    }
                    

                }
                else if(real_p<5){//mn>=5
                    std::array<size_t,N>new_begin_idx=begin_idx,new_end_idx=end_idx;
                    for(size_t k=1;k<p;k+=2){
                        new_end_idx[direction3]=new_begin_idx[direction3]=math_begin_idx3+k*math_stride;
                        predict_error+=block_interpolation_2d_crossblock(data,  new_begin_idx, new_end_idx,std::array<size_t,2>{direction1,direction2}
                                                                        ,math_stride,interp_func,pb,std::array<float,2>{coeff_x_xy,coeff_y_xy},meta,cross_block,tuning);
                    }
                    return predict_error;

                }
                size_t stride3x1=3*stride1,stride3x2=3*stride2,stride3x3=3*stride3,stride2x1=2*stride1,stride2x2=2*stride2,stride2x3=2*stride3;
                size_t math_stride2x=2*math_stride;
                //size_t math_stride3x=3*math_stride;
                //adaptive todo
              
                   
                size_t i,j,k;
                T *d;
                size_t i_start=(cross_back and math_begin_idx1>=math_stride2x)?1:3;
                size_t j_start=(cross_back and math_begin_idx2>=math_stride2x)?1:3;
                size_t k_start=(cross_back and math_begin_idx3>=math_stride2x)?1:3;
                bool i1_b=(i_start==3) and n>4;
                bool j1_b=(j_start==3) and m>4;
                bool k1_b=(k_start==3) and p>4;
                


                if(!meta.adjInterp){
                    

                    
                    for (i = i_start; i + 3 < n; i += 2) {
                        for(j=j_start;j+3<m;j+=2){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin + i* stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }         
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + i* stride1+j*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);  
                            }
                        }
                        //j=1
                        if(j1_b){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin + i* stride1+stride2+k*stride3;

                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                             }
                            //k=1
                            if(k1_b){
                                d = data + begin + i* stride1+stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }

                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + i* stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + i* stride1+stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                             }
                            //k=1
                            if(k1_b){
                                d = data + begin + i* stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + i* stride1+j*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                        }
                        if(m%2 ==0){//j=m-1
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin + i* stride1+(m-1)*stride2+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin + i* stride1+(m-1)*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + i* stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + i* stride1+(m-1)*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                        }
                    }
                    //i=1
                    if(i1_b){
                        for(j=j_start;j+3<m;j+=2){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin +  stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                                +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin +  stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin +  stride1+j*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=1
                        if(j1_b){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin +  stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin + stride1+stride2+stride3;
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                    +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                    +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin +  stride1+stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin +  stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin + stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + stride1+j*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        if(m%2 ==0){//j=m-1
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin +  stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin +  stride1+(m-1)*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + stride1+(m-1)*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                        }
                    }
                    //i= n-3 or n-2
                    if(i<n-1){
                        for(j=j_start;j+3<m;j+=2){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                                +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin +  i*stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin +  i*stride1+j*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=1
                        if(j1_b){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin +  i*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin + i*stride1+stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  i*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                                +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin +  i*stride1+stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin +  i*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin + i*stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  i*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                                +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + i*stride1+j*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        if(m%2 ==0){//j=m-1
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin +  i*stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin +  i*stride1+(m-1)*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + i*stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + i*stride1+(m-1)*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                            }
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=j_start;j+3<m;j+=2){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin + (n-1)* stride1+j*stride2+k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin +  (n-1)*stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +(n-1)* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin +  (n-1)*stride1+j*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=1
                        if(j1_b){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin + (n-1)*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin + (n-1)*stride1+stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  (n-1)*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                                +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + (n-1)*stride1+stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin + (n-1)*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin + (n-1)*stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))
                                                                                +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  (n-1)*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                                +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + (n-1)*stride1+j*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        if(m%2 ==0){//j=m-1
                            for(k=k_start;k+3<p;k+=2){
                                d = data + begin +  (n-1)*stride1+(m-1)*stride2+k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k1_b){
                                d = data + begin +  (n-1)*stride1+(m-1)*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }

                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + (n-1)*stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            if(p%2==0){
                                d = data + begin + (n-1)*stride1+(m-1)*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }
                }
                else{


                    //first half (non-adj) 

                    size_t temp_k_start=k_start;
                    size_t ks1=3,ks2=(temp_k_start==1)?1:5;

                    

                    for (i = i_start; i + 3 < n; i += 2) {
                        for(j=j_start;j+3<m;j+=2){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i* stride1+j*stride2+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                            }
                        
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                            }

                        }
                        //j=1
                        if(j1_b){
                            k_start=(i+1)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+stride2+k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i* stride1+stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }

                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1+stride2+k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                        }

                        //j=m-3 or m-2 or m-1
                        if(j<m-1){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin + i* stride1+j*stride2+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                           
                        }
                        if(m%2==0){
                            k_start=(i+m-1)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin + i* stride1+(m-1)*stride2+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }

                        }
                    }
                    //i=1
                    if(i1_b){
                    
                        for(j=j_start;j+3<m;j+=2){
                            k_start=(1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  stride1+j*stride2+k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                                +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  stride1+j*stride2+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if (k<p){
                                d = data + begin + stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
        
                        }
                        //j=1 (i+j=2)
                        if(j1_b){
                            k_start=ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  stride1+stride2+k*stride3;
                              
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + stride1+stride2+stride3;
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                    +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                    +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  stride1+stride2+k*stride3;
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                    +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                    +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }

                            //k=p-1
                            else if(k<p){
                                d = data + begin +  stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            k_start=(1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        //j=m-1 (i=1)
                        if(m%2 ==0){
                            k_start=(m%4==0)?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  stride1+(m-1)*stride2+k*stride3;
                              
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  stride1+(m-1)*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + stride1+(m-1)*stride2+k*stride3;
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                    +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                        }
                    }
                    //i= n-3 or n-2
                    if(i<n-1){
                        for(j=j_start;j+3<m;j+=2){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                                +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  i*stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin +i* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=1
                        if(j1_b){
                            k_start=(i+1)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  i*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i*stride1+stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  i*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin +  i*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  i*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i*stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                                +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  i*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + i*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(i+m-1)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  i*stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  i*stride1+(m-1)*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + i*stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + i*stride1+(m-1)*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                            }
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=j_start;j+3<m;j+=2){
                            k_start=(n-1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + (n-1)* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z_yz*interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  (n-1)*stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin +(n-1)* stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            
                        }
                        //j=1
                        if(j1_b){
                            k_start=(n%4==0)?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + (n-1)*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + (n-1)*stride1+stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  (n-1)*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + (n-1)*stride1+stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            k_start=(n-1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + (n-1)*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + (n-1)*stride1+j*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                                +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  (n-1)*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                                +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + (n-1)*stride1+j*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(n+m-2)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  (n-1)*stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  (n-1)*stride1+(m-1)*stride2+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + (n-1)*stride1+(m-1)*stride2+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + (n-1)*stride1+(m-1)*stride2+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }


                    //}
                    
                    //second half (adj)
                    ks1=(temp_k_start==1)?1:5;
                    ks2=3;

                    for (i = i_start; i + 3 < n; i += 2) {
                        for(j=j_start;j+3<m;j+=2){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i* stride1 +j*stride2 +stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ,mode);
                            }
                        
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ,mode);                            }

                        }
                        //j=1
                        if(j1_b){
                            k_start=(i+1)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+  stride2 +k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ,mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i* stride1 +stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }

                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1 +stride2 +k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                        }

                        //j=m-3 or m-2 
                        if(j<m-1){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+  j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin + i* stride1 +j*stride2 +stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                            
                        }
                        //j=m-1
                        if(m%2==0){
                            k_start=(i+m-1)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1+  (m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin + i* stride1 +(m-1)*stride2 +stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin + i* stride1 +(m-1)*stride2 +k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }

                        }
                    }
                    //i=1
                    if(i1_b){
                        for(j=j_start;j+3<m;j+=2){
                            k_start=(1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  stride1 +j*stride2 +k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2))
                                                                                +coeff_z_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  stride1 +j*stride2 +stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if (k<p){
                                d = data + begin + stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
        
                        }
                        //j=1 (i+j=2)
                        if(j1_b){
                            k_start=ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  stride1  +stride2 +k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + stride1 +stride2 +stride3;
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                    +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                                    +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  stride1 +stride2 +k*stride3;
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                    +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                                    +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                            }

                            //k=p-1
                            else if(k<p){
                                d = data + begin +  stride1 +stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                +coeff_y_xy*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            k_start=(1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  stride1+ j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + stride1 +j*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                +coeff_y_xy*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        //j=m-1 (i=1)
                        if(m%2 ==0){
                            k_start=(m%4==0)?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  stride1+ (m-1)*stride2 +k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  stride1 +(m-1)*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                +coeff_z_xz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + stride1 +(m-1)*stride2 +k*stride3;
                                    predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                    +coeff_z_xz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + stride1 +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),mode);
                            }
                        }
                    }
                    //i= n-3 or n-2
                    if(i<n-1){
                        for(j=j_start;j+3<m;j+=2){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + i* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2))
                                                                                +coeff_z_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  i*stride1 +j*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin +i* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=1
                        if(j1_b){
                            k_start=(i+1)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  i*stride1+  stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i*stride1 +stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                                +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  i*stride1 +stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                                +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin +  i*stride1 +stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)),mode);
                            }
                        }
                        if(j<m-1){
                            //j=m-3 or m-2
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  i*stride1  +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + i*stride1 +j*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                                +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  i*stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                                +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + i*stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(i+m-1)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  i*stride1 +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  i*stride1 +(m-1)*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + i*stride1 +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + i*stride1 +(m-1)*stride2 +(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),mode);
                            }
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=j_start;j+3<m;j+=2){
                            k_start=(n-1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + (n-1)* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) 
                                                                                +coeff_z_yz*interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  (n-1)*stride1 +j*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin +(n-1)* stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            
                        }
                        //j=1
                        if(j1_b){
                            k_start=(n%4==0)?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + (n-1)*stride1+ stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + (n-1)*stride1 +stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                                +coeff_z_yz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  (n-1)*stride1 +stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                                +coeff_z_yz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + (n-1)*stride1 +stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                            }
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            k_start=(n-1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin + (n-1)*stride1+  j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin + (n-1)*stride1 +j*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                                +coeff_z_yz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin +  (n-1)*stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                                +coeff_z_yz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + (n-1)*stride1 +j*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                            }
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(n+m-2)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin +  (n-1)*stride1+  (m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_cubic_adj(meta.cubicSplineType,*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin +  (n-1)*stride1 +(m-1)*stride2 +stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin + (n-1)*stride1 +(m-1)*stride2 +k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin + (n-1)*stride1 +(m-1)*stride2 +(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }
                }
            }

           // quant_index=quant_idx;
            return predict_error;


        }


        double block_interpolation_2d_cross(T *data, size_t begin1, size_t end1, size_t begin2, size_t end2, size_t stride1,size_t stride2, const std::string &interp_func, const PredictorBehavior pb,const Interp_Meta &meta,int tuning=0) {
            size_t n = (end1 - begin1) / stride1 + 1;
            if (n <= 1) {
                return 0;
            }
            size_t m = (end2 - begin2) / stride2 + 1;
            if (m <= 1) {
                return 0;
            }

            double predict_error = 0;
            size_t stride3x1=3*stride1,stride3x2=3*stride2,stride5x1=5*stride1,stride5x2=5*stride2;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
           // size_t quant_idx=quant_index;
           
            if (interp_func == "linear"|| n<5 || m<5 ) {//nmcond temp added
                
                for (size_t i = 1; i + 1 < n; i += 2) {
                    for(size_t j=1;j+1<m;j+=2){
                        T *d = data + begin1 + i* stride1+begin2+j*stride2;
                       
                        predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);

                    }
                    if(m%2 ==0){
                        T *d = data + begin1 + i * stride1+begin2+(m-1)*stride2;

                        if(i<3 or i+3>=n or m<4)
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1 - stride2), *(d + stride1 - stride2)),mode);//this is important. Not sure whether it is good.
                        else
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_linear1(*(d - stride3x1 - stride3x2),*(d - stride1 - stride2))
                                                            , interp_linear1(*(d + stride3x1 - stride3x2),*(d + stride1 - stride2))),mode);//this is important. Not sure whether it is good.
                    }
                }
                if (n % 2 == 0) {
                    for(size_t j=1;j+1<m;j+=2){

                        T *d = data + begin1 + (n-1) * stride1+begin2+j*stride2;
                        if(n<4 or j<3 or j+3>=m)
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1 - stride2), *(d - stride1 + stride2)),mode);//this is important. Not sure whether it is good.
                        else
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_linear1(*(d - stride3x1 - stride3x2),*(d - stride1 - stride2))
                                                            , interp_linear1(*(d - stride3x1 + stride3x2),*(d - stride1 + stride2))),mode);//this is important. Not sure whether it is good.
                    }
                    if(m%2 ==0){
                        T *d = data + begin1 + (n-1) * stride1+begin2+(m-1)*stride2;
                        if(n<4 or m<4)
                            predict_error+=quantize_integrated(d - data, *d, *(d - stride1 - stride2),mode);//this is important. Not sure whether it is good.
                        else
                            predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d - stride3x1 - stride3x2),*(d - stride1 - stride2) ),mode);//this is important. Not sure whether it is good.
                    }          
                }
                    
            }
            else{//cubic
                //adaptive todo
                size_t i,j;
                T *d;
                for (i = 3; i + 3 < n; i += 2) {
                    for(j=3;j+3<m;j+=2){
                        d = data + begin1 + i* stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1-stride3x2), *(d - stride1-stride2), *(d + stride1+stride2), *(d + stride3x1+stride3x2))
                                                        ,interp_cubic(meta.cubicSplineType,*(d +stride3x1- stride3x2), *(d +stride1- stride2), *(d -stride1+ stride2), *(d -stride3x1+ stride3x2)) ),mode);
                    }
                    //j=1
                    d = data + begin1 + i* stride1+ begin2+stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1( *(d - stride1-stride2),*(d + stride1+stride2),*(d + stride3x1+stride3x2) )
                                                    ,interp_quad_1( *(d + stride1-stride2),*(d - stride1+stride2),*(d - stride3x1+stride3x2) ) ),mode);
                                       
                    //j=m-3 or m-2
                    d = data +begin1 + i* stride1+ begin2+j*stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2( *(d - stride3x1-stride3x2),*(d - stride1-stride2),*(d + stride1+stride2) )
                                                    ,interp_quad_2( *(d + stride3x1-stride3x2),*(d + stride1-stride2),*(d - stride1+stride2) ) ),mode);
                    
                    //j=m-1
                    if(m%2 ==0){
                        d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                        if(i>=5 and i+5<n)
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3( *(d - stride5x1-stride5x2),*(d - stride3x1-stride3x2),*(d - stride1-stride2) )
                                                   ,interp_quad_2( *(d + stride5x1-stride5x2),*(d + stride3x1-stride3x2),*(d + stride1-stride2) ) ),mode);
                        else
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_linear1( *(d - stride3x1-stride3x2),*(d - stride1-stride2) )
                                                    ,interp_linear1(*(d + stride3x1-stride3x2),*(d + stride1-stride2) ) ),mode);
                    }
                }
               
                //i=1
                for(j=3;j+3<m;j+=2){
                    d = data + begin1 + stride1+begin2+j*stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1( *(d - stride1-stride2),*(d + stride1+stride2),*(d + stride3x1+stride3x2) )
                                                    ,interp_quad_1( *(d - stride1+stride2),*(d + stride1-stride2),*(d + stride3x1-stride3x2) ) ),mode);
                }
                //j=1

                d = data + begin1 + stride1+ begin2+stride2;
                //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);//2d linear
                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad1( *(d - stride1-stride2),*(d + stride1+stride2),*(d + stride3x1+stride3x2) )
                //                                    ,interp_linear( *(d + stride1-stride2),*(d - stride1+stride2) ) ),mode);//2d linear+quad
                predict_error+=quantize_integrated(d - data, *d, interp_quad_1( *(d - stride1-stride2),*(d + stride1+stride2),*(d + stride3x1+stride3x2) ),mode);//1d quad

                //j=m-3 or m-2
                d = data +begin1 + stride1+ begin2+j*stride2;

                //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);//2d linear
                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad2( *(d + stride3x1-stride3x2),*(d + stride1-stride2),*(d - stride1+stride2) )
                //                                    ,interp_linear( *(d - stride1-stride2),*(d + stride1+stride2) ) ),mode);//2d linear+quad
                predict_error+=quantize_integrated(d - data, *d, interp_quad_2( *(d + stride3x1-stride3x2),*(d + stride1-stride2),*(d - stride1+stride2) ),mode);//1d quad
                
                //j=m-1
                if(m%2 ==0){
                    d = data + begin1 + stride1+begin2+(m-1)*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1-stride2), *(d + stride1-stride2)),mode);//1d linear
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d + stride3x1-stride3x2), *(d + stride1-stride2)),mode);//1d cross linear
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1-stride2), *(d + stride1-stride2),*(d + stride3x1-stride2)),mode);//1d quad

                }

                //i= n-3 or n-2
                for(j=3;j+3<m;j+=2){
                   
                    d = data + begin1 + i*stride1+begin2+j*stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2( *(d - stride3x1-stride3x2),*(d - stride1-stride2),*(d + stride1+stride2) )
                                                    ,interp_quad_2( *(d - stride3x1+stride3x2),*(d - stride1+stride2),*(d + stride1-stride2) ) ),mode);

                }
                //j=1
                d = data + begin1 + i*stride1+ begin2+stride2;
                //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);//2d linear
                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad1( *(d + stride1-stride2),*(d - stride1+stride2),*(d - stride3x1+stride3x2) )
                //                                    ,interp_linear( *(d - stride1-stride2),*(d + stride1+stride2) ) ),mode);//2d linear+quad
                predict_error+=quantize_integrated(d - data, *d, interp_quad_1( *(d + stride1-stride2),*(d - stride1+stride2),*(d - stride3x1+stride3x2) ),mode);//1d quad
                
                //j=m-3 or m-2
                d = data +begin1 + i*stride1+ begin2+j*stride2;
           
                //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);//2d linear
                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad2( *(d - stride3x1-stride3x2),*(d - stride1-stride2),*(d + stride1+stride2) )
                //                                    ,interp_linear( *(d + stride1-stride2),*(d - stride1+stride2) ) ),mode);//2d linear+quad
                predict_error+=quantize_integrated(d - data, *d, interp_quad_2( *(d - stride3x1-stride3x2),*(d - stride1-stride2),*(d + stride1+stride2) ),mode);//1d quad
                
                //j=m-1
                if(m%2 ==0){
                    d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1-stride2), *(d + stride1-stride2)),mode);//1d linear
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d + stride3x1-stride3x2), *(d + stride1-stride2)),mode);//1d cross linear
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d + stride1-stride2), *(d - stride1-stride2),*(d - stride3x1-stride2)),mode);//1d quad
                }

                //i=n-1 (odd)
                if (n % 2 == 0) {
                    for(j=3;j+3<m;j+=2){
                        d = data + begin1 + (n-1)*stride1+begin2+j*stride2;
                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_linear1(*(d - stride3x1-stride3x2), *(d - stride1-stride2)) ,interp_linear1(*(d - stride3x1+stride3x2), *(d - stride1+stride2)) ),mode);//2d cross
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride1-stride3x2), *(d -  stride1-stride2), *(d - stride1+ stride2), *(d - stride1+ stride3x2)),mode);//1d cubic


                    }
                    //j=1
                    d = data + begin1 + (n-1)*stride1+ begin2+stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear1( *(d - stride3x1+stride3x2), *(d - stride1+stride2)),mode);//1d linear
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1-stride2), *(d - stride1+stride2),*(d - stride1+stride3x2)),mode);//1d quad

                    //j=m-3 or m-2
                    d = data +begin1 + (n-1)*stride1+ begin2+j*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear1( *(d - stride3x1-stride3x2), *(d - stride1-stride2)),mode);//1d linear
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride1-stride3x2), *(d - stride1-stride2),*(d - stride1+stride2)),mode);//1d quad
                    //j=m-1
                    if(m%2 ==0){
                        d = data + begin1 + (n-1) * stride1+begin2+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_3(*(d - stride5x1-stride5x2), *(d - stride3x1-stride3x2), *(d - stride1-stride2)),mode);
                    }
                }         
            }

           // quant_index=quant_idx;
            return predict_error;
        }

        double block_interpolation_2d_aftercross(T *data, size_t begin1, size_t end1, size_t begin2, size_t end2, size_t stride1,size_t stride2, const std::string &interp_func, const PredictorBehavior pb,const Interp_Meta &meta,int tuning=0) {
            size_t n = (end1 - begin1) / stride1 + 1;
            
            size_t m = (end2 - begin2) / stride2 + 1;
            if (n<=1&& m <= 1) {
                return 0;
            }
            double predict_error = 0;
            size_t stride3x1=3*stride1,stride3x2=3*stride2,stride5x1=5*stride1,stride5x2=5*stride2;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
           // size_t quant_idx=quant_index;
           
            if (interp_func == "linear"|| n<5 || m<5 ) {//nmcond temp added
                size_t i,j;
                for (i = 1; i + 1 < n; i += 1) {
                    for(j=1+(i%2);j+1<m;j+=2){
                        T *d = data + begin1 + i* stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 ), *(d + stride1 ),*(d  + stride2), *(d  - stride2)),mode);

                    }

                    //j=0
                    if(i%2==1 and begin2==0){
                        T *d = data + begin1 + i* stride1+begin2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1 ), *(d + stride1 )),mode);
                    }
                    //j=m-1, j wont be 0
                    if(j==m-1){
                        T *d = data + begin1 + i* stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1 ), *(d + stride1 )),mode);
                    }
                }
                //i=0
                if(begin1==0){
                    for(j=1;j+1<m;j+=2){
                        T *d = data + begin1 +begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d  + stride2), *(d  - stride2)),mode);
                    }
                
                //j=m-1, j wont be 0
                    if(j==m-1){
                        T *d = data + begin1 +begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, *(d-stride2),mode);//for simplicity,may extend to 2d.
                    }
                }
                //i=n-1
                if(n>1){
                    for(j=1+(n-1)%2;j+1<m;j+=2){
                        T *d = data + begin1 +(n-1)*stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d  + stride2), *(d  - stride2)),mode);
                    }
                    //j=0
                    if((n-1)%2==1 and begin2==0){
                        
                        T *d = data + begin1 +(n-1)*stride1+begin2;

                        predict_error+=quantize_integrated(d - data, *d, *(d-stride1),mode);//for simplicity,may extend to 2d.
                    }
                    //j=m-1, j wont be 0
                    if( j==m-1){
                        
                        T *d = data + begin1 +(n-1)*stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, *(d-stride1),mode);//for simplicity,may extend to 2d.
                    }
                }
                    
            }
            else{//cubic
                //adaptive todo
                size_t i,j;
                T *d;
                for (i = 3; i + 3 < n; i += 1) {
                    for(j=3+(i%2);j+3<m;j+=2){
                        d = data + begin1 + i* stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                        ,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d+ stride2), *(d + stride3x2)) ),mode);
                    }
                    //j=0
                    if(i%2==1 and begin2==0){
                        d = data + begin1 + i* stride1+begin2;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    }
                    //j=1 or 2 
                    d = data + begin1 + i* stride1+begin2+(1+(i%2))*stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    
                    //j=m-3 or m-2, j wont be 2.
                    d = data + begin1 + i* stride1+begin2+j*stride2;
                    predict_error+=quantize_integrated(d - data, *d,  interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                   
                    //j=m-1
                    if(j+2==m-1){
                        d = data + begin1 + i* stride1+begin2+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d,  interp_cubic(meta.cubicSplineType,*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    }
                   
                }
                std::vector<size_t> boundary_is=(n>5)?std::vector<size_t>{0,1,2,n-3,n-2,n-1}:std::vector<size_t>{0,1,2,n-2,n-1};
                
                for(auto ii:boundary_is){
                    if(ii==0 and begin1!=0)
                        continue;
                    for(j=3+(ii%2);j+3<m;j+=2){
                        d = data + begin1+ii*stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d,interp_cubic(meta.cubicSplineType,*(d - stride3x2), *(d - stride2), *(d+ stride2), *(d + stride3x2) ),mode);
                    }

                    std::vector<size_t> boundary_js=(ii%2)?std::vector<size_t>{0,2,j}:std::vector<size_t>{1,j};
                    if(j+2==m-1)
                        boundary_js.push_back(m-1);

                    for(auto jj:boundary_js){
                        if(begin2!=0 and jj==0)
                            continue;
                        d = data + begin1 + ii* stride1+begin2+jj*stride2;
                        T v1;
                        if(ii==0){
                            if(n>5)
                                v1=interp_quad_3(*(d + stride5x1), *(d+ stride3x1), *(d + stride1) );
                            else
                                v1=interp_linear1( *(d+ stride3x1), *(d + stride1) );
                        }
                        else if(ii==1)
                            v1=interp_quad_1(*(d - stride1), *(d+ stride1), *(d + stride3x1) );
                        else if (ii==n-2)
                            v1=interp_quad_2(*(d - stride3x1), *(d- stride1), *(d + stride1) );
                        else if (ii==n-1){
                            if(n>5)
                                v1=interp_quad_3(*(d - stride5x1), *(d- stride3x1), *(d - stride1) );
                            else
                                v1=interp_linear1( *(d- stride3x1), *(d - stride1) );
                        }
                        else{//i==2 or n-3
                            if(n==5)
                                v1=interp_linear(*(d - stride1), *(d+ stride1));
                            else if (ii==2)
                                v1=interp_quad_1(*(d - stride1), *(d+ stride1), *(d + stride3x1) );
                            else
                                v1=interp_quad_2(*(d - stride3x1), *(d- stride1), *(d + stride1) );
                        }

                        T v2;
                        if(jj==0){
                            if(m>5)
                                v2=interp_quad_3( *(d + stride5x2), *(d+ stride3x2), *(d + stride2) );
                            else
                                v2=interp_linear1( *(d+ stride3x2), *(d + stride2) );
                        }
                        else if (jj==1 or jj==2)
                            v2=interp_quad_1( *(d - stride2), *(d+ stride2), *(d + stride3x2) );
                        else if(jj==m-1){
                            if(m>5)
                                v2=interp_quad_3( *(d - stride5x2), *(d- stride3x2), *(d - stride2) );
                            else
                                v2=interp_linear1( *(d- stride3x2), *(d - stride2) );
                        }
                        else
                            v2=interp_quad_2( *(d - stride3x2), *(d- stride2), *(d + stride2) );

                        predict_error+=quantize_integrated(d - data, *d,interp_linear( v1,v2 ),mode);

                    }

                }
            }

           // quant_index=quant_idx;
            return predict_error;
        }
        
        template<uint NN = N>
        typename std::enable_if<NN == 1, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const Interp_Meta & meta, size_t stride = 1,int tuning=0,int cross_block=0) {//regressive to reduce into meta.
            if(!cross_block)
                return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb,meta,tuning);
            else
                return block_interpolation_1d_crossblock(data, begin, end, 0,stride, interp_func, pb,meta,1,tuning);

        }


        template<uint NN = N>
        typename std::enable_if<NN == 2, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const Interp_Meta & meta, size_t stride = 1,int tuning=0,int cross_block=0) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            //bool full_adjacent_interp=false;
            uint8_t paradigm=meta.interpParadigm;
            uint8_t direction=meta.interpDirection;
            assert(direction<2);
            if(paradigm==0){
                const std::array<size_t, N> dims = dimension_sequences[direction];
                std::array<size_t, N>steps;
                std::array<size_t, N> begin_idx=begin,end_idx=end;
                steps[dims[0]]=1;
                begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
                steps[dims[1]]=stride2x;
                
                
                predict_error += block_interpolation_1d_crossblock_2d(data, begin_idx,
                                                                    end_idx,dims[0],steps,
                                                                    stride , interp_func, pb,meta,cross_block,tuning);
                

                begin_idx[dims[1]]=begin[dims[1]];
                begin_idx[dims[0]]=(begin[dims[0]] ? begin[dims[0]] + stride : 0);
                steps[dims[0]]=stride;
               
                predict_error += block_interpolation_1d_crossblock_2d(data, begin_idx,
                                                                    end_idx,dims[1],steps,
                                                                    stride , interp_func, pb,meta,cross_block,tuning);
            }
            
            else{// if(paradigm<3){//md or hd
                const std::array<size_t, N> dims = dimension_sequences[0];
                std::array<float,2>dim_coeffs={meta.dimCoeffs[0],meta.dimCoeffs[1]};

                std::array<size_t, N>steps;
                std::array<size_t, N> begin_idx=begin,end_idx=end;
                steps[dims[0]]=1;
                begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
                steps[dims[1]]=stride2x;

                predict_error += block_interpolation_1d_crossblock_2d(data, begin_idx,
                                                                    end_idx,dims[0],steps,
                                                                    stride , interp_func, pb,meta,cross_block,tuning);
                
            
                begin_idx[dims[1]]=begin[dims[1]];

                begin_idx[dims[0]]=(begin[dims[0]] ? begin[dims[0]] + stride2x : 0);
           
                steps[dims[0]]=stride2x;
             

                predict_error += block_interpolation_1d_crossblock_2d(data, begin_idx,
                                                                    end_idx,dims[1],steps,
                                                                    stride , interp_func, pb,meta,cross_block,tuning);
              
                begin_idx=begin,end_idx=end;
                
                predict_error += block_interpolation_2d_crossblock(data, begin_idx,
                                                            end_idx,dims,
                                                            stride , interp_func, pb,dim_coeffs,meta,cross_block,tuning);

                

            }
            return predict_error;
        }




        template<uint NN = N>
        typename std::enable_if<NN == 3, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const Interp_Meta & meta, size_t stride = 1,int tuning=0,int cross_block=0) {//cross block: 0 or conf.num

            double predict_error = 0;
            size_t stride2x = stride * 2;
            uint8_t paradigm=meta.interpParadigm;
            uint8_t direction=meta.interpDirection;
            bool fallback_2d=direction>=6;
            if (fallback_2d){
                direction-=6;
            }
            assert(direction<6);
            if(paradigm==0){
                const std::array<size_t, N> dims = dimension_sequences[direction];
                if(!fallback_2d){
                        std::array<size_t, N>steps;
                        std::array<size_t, N> begin_idx=begin,end_idx=end;
                        steps[dims[0]]=1;
                        begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
                        begin_idx[dims[2]]=(begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
                        steps[dims[1]]=stride2x;
                        steps[dims[2]]=stride2x;
                        predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                            end_idx,dims[0],steps,
                                                                            stride , interp_func, pb,meta,cross_block,tuning);
                        
                        begin_idx[dims[1]]=begin[dims[1]];

                        begin_idx[dims[0]]=(begin[dims[0]] ? begin[dims[0]] + stride : 0);
                        steps[dims[0]]=stride;

                        predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                            end_idx,dims[1],steps,
                                                                            stride , interp_func, pb,meta,cross_block,tuning);

                        
                        begin_idx[dims[2]]=begin[dims[2]];

                        begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride : 0);
                        steps[dims[1]]=stride;
                        predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                            end_idx,dims[2],steps,
                                                                            stride , interp_func, pb,meta,cross_block,tuning);
                }
                else{
                    std::array<size_t, N>steps;
                    std::array<size_t, N> begin_idx=begin,end_idx=end;
                    steps[dims[1]]=1;
                    begin_idx[dims[0]]=(begin[dims[0]] ? begin[dims[0]] + 1 : 0);
                    begin_idx[dims[2]]=(begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
                    steps[dims[0]]=1;
                    steps[dims[2]]=stride2x;


                    predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                        end_idx,dims[1],steps,
                                                                        stride , interp_func, pb,meta,cross_block,tuning);
                   
            

                    begin_idx[dims[2]]=begin[dims[2]];
                    begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride : 0);
                    steps[dims[1]]=stride;
                    predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                        end_idx,dims[2],steps,
                                                                        stride , interp_func, pb,meta,cross_block,tuning);
                }
            }
            
            else {//if (paradigm==1){
                std::array<float,3>dim_coeffs=meta.dimCoeffs;
                if(!fallback_2d){
                    const std::array<size_t, N> dims = dimension_sequences[0];
                    std::array<size_t, N>steps;
                    std::array<size_t, N> begin_idx=begin,end_idx=end;
                    steps[dims[0]]=1;
                    begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
                    begin_idx[dims[2]]=(begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
                    steps[dims[1]]=stride2x;
                    steps[dims[2]]=stride2x;
                    predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                        end_idx,dims[0],steps,
                                                                        stride , interp_func, pb,meta,cross_block,tuning);
                    
                    begin_idx[dims[1]]=begin[dims[1]];
                    begin_idx[dims[0]]=(begin[dims[0]] ? begin[dims[0]] + stride2x : 0);
                    steps[dims[0]]=stride2x;


                    predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                        end_idx,dims[1],steps,
                                                                        stride , interp_func, pb,meta,cross_block,tuning);
                    begin_idx[dims[2]]=begin[dims[2]];

                    begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
                    steps[dims[1]]=stride2x;


                    predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                        end_idx,dims[2],steps,
                                                                        stride , interp_func, pb,meta,cross_block,tuning);
                    begin_idx=begin,end_idx=end;
                    begin_idx[dims[2]]=(begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
                    steps[dims[2]]=stride2x;
                    predict_error += block_interpolation_2d_crossblock_3d(data, begin_idx,
                                                                end_idx,std::array<size_t,2>{dims[0],dims[1]},steps,
                                                                stride , interp_func, pb,std::array<float,2>{dim_coeffs[dims[0]],dim_coeffs[dims[1]]},meta,cross_block,tuning);
                    begin_idx[dims[2]]=begin[dims[2]];
                    begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
                    steps[dims[1]]=stride2x;  
                        predict_error += block_interpolation_2d_crossblock_3d(data, begin_idx,
                                                                end_idx,std::array<size_t,2>{dims[0],dims[2]},steps,
                                                                stride , interp_func, pb,std::array<float,2>{dim_coeffs[dims[0]],dim_coeffs[dims[2]]},meta,cross_block,tuning);
                    begin_idx[dims[1]]=begin[dims[1]];
                    begin_idx[dims[0]]=(begin[dims[0]] ? begin[dims[0]] + stride2x : 0);
                    steps[dims[0]]=stride2x;
                        predict_error += block_interpolation_2d_crossblock_3d(data, begin_idx,
                                                                end_idx,std::array<size_t,2>{dims[1],dims[2]},steps,
                                                                stride , interp_func, pb,std::array<float,2>{dim_coeffs[dims[1]],dim_coeffs[dims[2]]},meta,cross_block,tuning);
                    begin_idx=begin,end_idx=end;
                    predict_error += block_interpolation_3d_crossblock(data, begin_idx,
                                                                end_idx,dims,
                                                                stride , interp_func, pb,dim_coeffs,meta,cross_block,tuning);
                                                                    
                }
                else{
                    const std::array<size_t, N> dims = dimension_sequences[direction];
                    std::array<size_t, N>steps;
                    std::array<size_t, N> begin_idx=begin,end_idx=end;
                    steps[dims[1]]=1;
                    begin_idx[dims[0]]=(begin[dims[0]] ? begin[dims[0]] + 1 : 0);
                    begin_idx[dims[2]]=(begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
                    steps[dims[0]]=1;
                    steps[dims[2]]=stride2x;
                    predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                        end_idx,dims[1],steps,
                                                                        stride , interp_func, pb,meta,cross_block,tuning);
                    

                    begin_idx[dims[2]]=begin[dims[2]];
                    begin_idx[dims[1]]=(begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
                    steps[dims[1]]=stride2x;
                    predict_error += block_interpolation_1d_crossblock_3d(data, begin_idx,
                                                                        end_idx,dims[2],steps,
                                                                        stride , interp_func, pb,meta,cross_block,tuning);
                    begin_idx[dims[1]]=begin[dims[1]];
                    predict_error += block_interpolation_2d_crossblock_3d(data, begin_idx,
                                                                end_idx,std::array<size_t,2>{dims[1],dims[2]},steps,
                                                                stride , interp_func, pb,std::array<float,2>{dim_coeffs[dims[1]],dim_coeffs[dims[2]]},meta,cross_block,tuning);

                }

                
            }
            return predict_error;
            
        }


        template<uint NN = N>
        typename std::enable_if<NN == 4, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const Interp_Meta & meta, size_t stride = 1,int tuning=0,int cross_block=0) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            //uint8_t paradigm=meta.interpParadigm;
            uint8_t direction=meta.interpDirection;
            assert(direction<24);
            //max_error = 0;
            const std::array<size_t, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning);
                    }
                }
            }
//            printf("%.8f ", max_error);
           // max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                    }
                }
            }
//            printf("%.8f ", max_error);
            //max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              begin[dims[2]] * dimension_offsets[dims[2]] +
                                              t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning);
                    }
                }
            }

//            printf("%.8f ", max_error);
          //  max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0); k <= end[dims[2]]; k += stride) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                begin[dims[3]] * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[3]] - begin[dims[3]]) *
                                                                dimension_offsets[dims[3]],
                                                                stride * dimension_offsets[dims[3]], interp_func, pb,meta,tuning);
                    }
                }
            }
//            printf("%.8f \n", max_error);
            return predict_error;
        }


        bool anchor=false;
        int interpolation_level = -1;
        uint blocksize;
        Interp_Meta interp_meta;

        double eb_ratio = 0.5;
        double alpha;
        double beta;
        std::vector<std::string> interpolators = {"linear", "cubic","quad"};
        int *quant_inds;
        //std::vector<bool> mark;
        size_t quant_index = 0; 
        size_t maxStep=0;
        //double max_error;
        Quantizer quantizer;
        size_t num_elements;

        std::array<size_t, N> global_dimensions;
        std::array<size_t, N> dimension_offsets;
        std::vector<std::array<size_t, N>> dimension_sequences;


        int levelwise_predictor_levels = 0;
        bool blockwiseTuning = false;
        std::vector <uint8_t> interpAlgo_list;
        std::vector <uint8_t> interpDirection_list;
        //std::vector <uint8_t> cubicSplineType_list;
        


        std::vector <Interp_Meta> interpMeta_list;
        int fixBlockSize;
            //int trimToZero;

        int frozen_dim=-1;
            
        int cross_block=0;
           
        

        std::vector<float> prediction_errors;//for test, to delete in final version. The float time is to match the vector in config.
        //int peTracking=0;//for test, to delete in final version

        size_t cur_level; //temp for "adaptive anchor stride";
        //size_t min_anchor_level;//temp for "adaptive anchor stride";
       // double anchor_threshold=0.0;//temp for "adaptive anchor stride";



    };
template <class T, uint N, class Quantizer>
QoZInterpolationDecomposition<T, N, Quantizer> make_decomposition_interpolation(const Config &conf, Quantizer quantizer) {
    return QoZInterpolationDecomposition<T, N, Quantizer>(conf, quantizer);
}


}


}  // namespace SZ3

#endif
