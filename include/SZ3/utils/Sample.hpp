

#ifndef SZ_SAMPLE_HPP
#define SZ_SAMPLE_HPP

namespace SZ3 {

    
    template<class T, uint N>
    inline void
    sample_block_4d(T *data, std::vector<T> & sampling_data , std::vector<size_t> &dims, std::vector<size_t> &starts,size_t block_size) {
        assert(dims.size() == N);
        assert(starts.size() == N);
        
        
        size_t sample_num = block_size*block_size*block_size;
        //std::vector<T> sampling_data(sample_num, 0);

       
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
        //return sampling_data;
    }
    


    template<class T, uint N>
    inline void
    profiling_block(T *data, std::vector<size_t> &dims, std::vector< std::vector<size_t> > &starts,size_t block_size, double abseb,size_t stride = 4) {
        assert(dims.size() == N);
        if (stride == 0)
            stride = block_size;
        if constexpr (N==3){
            size_t dimx=dims[0],dimy=dims[1],dimz=dims[2],dimyz=dimy*dimz;
            
            for (size_t i = 0; i < dimx-block_size; i+=block_size) {
                for (size_t j = 0; j < dimy-block_size; j+=block_size) {
                    for (size_t k = 0; k < dimz-block_size; k+=block_size) {
                        //std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                        size_t start_idx=i*dimyz+j*dimz+k;
                        T min=data[start_idx];
                        T max=data[start_idx];
                        for (int ii=0;ii<=block_size;ii+=stride){
                            for(int jj=0;jj<=block_size;jj+=stride){
                                for (int kk=0;kk<=block_size;kk+=stride){
                                    size_t cur_idx=start_idx+ii*dimyz+jj*dimz+kk;
                                    T cur_value=data[cur_idx];
                                    if (cur_value<min)
                                        min=cur_value;
                                    else if (cur_value>max)
                                        max=cur_value;

                                }
                            }
                        }
                        if (max-min>abseb){
                            size_t a[3]={i,j,k};
                            starts.push_back(std::vector<size_t>(a,a+3));
                        }

                    }
                }
            }
        }
        else if constexpr (N==2){
             size_t dimx=dims[0],dimy=dims[1];
        
            for (size_t i = 0; i < dimx-block_size; i+=block_size) {
                for (size_t j = 0; j < dimy-block_size; j+=block_size) {
                    
                    size_t start_idx=i*dimy+j;
                    T min=data[start_idx];
                    T max=data[start_idx];
                    for (int ii=0;ii<=block_size;ii+=stride){
                        for(int jj=0;jj<=block_size;jj+=stride){
                               
                            size_t cur_idx=start_idx+ii*dimy+jj;
                            T cur_value=data[cur_idx];
                            if (cur_value<min)
                                min=cur_value;
                            else if (cur_value>max)
                                max=cur_value;

                        }
                    }
                        
                    if (max-min>abseb){
                         size_t a[2]={i,j};
                        starts.push_back(std::vector<size_t>(a,a+2));
                    }


                        
                        
                }
            }

        }
        //current has a problem. May return no blocks. Thinking how to better solve it.
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
       // return sampling_data;
    }
    


    template<class T, uint N>
    inline void
    sample_blocks(T *data, std::vector<T> & sampling_data, std::vector<size_t> &dims, std::vector<size_t> &starts,size_t block_size) {
        assert(dims.size() == N);
        assert(starts.size() == N);
        if constexpr (N==3){
        
            size_t sample_num = block_size*block_size*block_size;
            sampling_data.resize(sample_num, 0);

            size_t startx=starts[0],starty=starts[1],startz=starts[2],dimx=dims[0],dimy=dims[1],dimz=dims[2];
            size_t square_block_size=block_size*block_size,dimyz=dimy*dimz;
            for (size_t i = 0; i < block_size; i++) {
                for (size_t j = 0; j < block_size; j++) {
                    for (size_t k = 0; k < block_size; k++) {
                        size_t sample_idx=i*square_block_size+j*block_size+k;
                        size_t idx=(i+startx)*dimyz+(j+starty)*dimz+k+startz;
                        sampling_data[sample_idx]=data[idx];
                        
                    }
                }
            }
        }
        else if constexpr (N==2){
            size_t sample_num = block_size*block_size;
            sampling_data.resize(sample_num, 0);
            size_t startx=starts[0],starty=starts[1],dimx=dims[0],dimy=dims[1];
            
            for (size_t i = 0; i < block_size; i++) {
                for (size_t j = 0; j < block_size; j++) {
                    
                    size_t sample_idx=i*block_size+j;
                    size_t idx=(i+startx)*dimy+(j+starty);
                    sampling_data[sample_idx]=data[idx];
                        
                    
                }
            }

        }
        else if constexpr(N==1){
            size_t sample_num = block_size;
            sampling_data.resize(sample_num, 0);

            size_t startx=starts[0],dimx=dims[0];
            
            for (size_t i = 0; i < block_size; i++) {
                
                    
                size_t sample_idx=i;
                size_t idx=(i+startx);
                sampling_data[sample_idx]=data[idx];
                        
                    
                
            }

        }
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
       // return sampling_data;
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
            totalblock_num *= static_cast<int>((dims[i]-1)/sampleBlockSize);
        }               
        size_t idx=0;   
        if(profiling){
            size_t num_filtered_blocks=starts.size();    
            
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
            
            size_t sample_stride = static_cast<size_t>(1.0/sample_rate);
            if(sample_stride<=0)
                sample_stride=1;
            if constexpr (N==2){                        
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
            else if constexpr (N==3){                  
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
    }

   
    


};


#endif
