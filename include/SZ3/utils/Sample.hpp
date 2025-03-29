

#ifndef SZ_SAMPLE_HPP
#define SZ_SAMPLE_HPP

#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/CoeffRegression.hpp"
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

   
    


};


#endif
