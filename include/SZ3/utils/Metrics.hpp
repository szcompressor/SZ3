

#ifndef SZ_METRICS_HPP
#define SZ_METRICS_HPP
#include<cmath>
namespace SZ3{
    namespace QoZ {
        
        inline double calc_PSNR(const double & rng, const double & mse) {
            return 20*log10(rng)-10*log10(mse);
        }

        inline double calc_SSIM(const double &xr,const double &xm,const double &xv2,const double &ym,const double &yv2,const double &cov) {
            double c1=0.0001,c2=0.0009;
            if (xr>0){
                c1*=(xr*xr);
                c2*=xr*xr;
            }
            return (2*xm*ym+c1)*(2*cov+c2)/(  (xm*xm+ym*ym+c1) *(xv2+yv2+c2) );
           
        }
        template <class T>
        void blockwise_profiling(const T *data, const std::vector<size_t> &dims, const std::vector<size_t> &starts,const size_t &blocksize,double & mean,double & sigma2,double & range){
            size_t N=dims.size();

            if(N==2){
                size_t dimy=dims[1],startx=starts[0],starty=starts[1],element_num=blocksize*blocksize;
                double sum=0,max,min;
                sigma2=0;
                size_t start_idx=startx*dimy+starty;
                max=min=data[start_idx];
                for(size_t i=startx;i<startx+blocksize;i++){
                    for(size_t j=starty;j<starty+blocksize;j++){
                        size_t cur_idx=i*dimy+j;
                        T value=data[cur_idx];
                        sum+=value;
                        
                        max=value>max?value:max;
                        min=value<min?value:min;

                    }

                }
                mean=sum/element_num;

                for(size_t i=startx;i<startx+blocksize;i++){
                    for(size_t j=starty;j<starty+blocksize;j++){
                        size_t cur_idx=i*dimy+j;
                        T value=data[cur_idx];
                        sigma2+=(value-mean)*(value-mean);

                    }

                }
                
                range=max-min;
                sigma2/=element_num;
            }


            else if(N==3){
                size_t dimy=dims[1],dimz=dims[2],startx=starts[0],starty=starts[1],startz=starts[2],element_num=blocksize*blocksize*blocksize;
                size_t dimyz=dimy*dimz;
                double sum=0,max,min;
                sigma2=0;
                size_t start_idx=startx*dimyz+starty*dimz+startz;
                max=min=data[start_idx];
                for(size_t i=startx;i<startx+blocksize;i++){
                    for(size_t j=starty;j<starty+blocksize;j++){
                        for(size_t k=startz;k<startz+blocksize;k++){
                            size_t cur_idx=i*dimyz+j*dimz+k;
                            T value=data[cur_idx];
                            sum+=value;
                           
                            max=value>max?value:max;
                            min=value<min?value:min;
                        }

                    }

                }
                mean=sum/element_num;
                for(size_t i=startx;i<startx+blocksize;i++){
                    for(size_t j=starty;j<starty+blocksize;j++){
                        for(size_t k=startz;k<startz+blocksize;k++){
                            size_t cur_idx=i*dimyz+j*dimz+k;
                            T value=data[cur_idx];
                            sigma2+=(value-mean)*(value-mean);

                            
                        }

                    }

                }

               
                range=max-min;
                sigma2/=element_num;
            }
        }

        template <class T>
        double calc_blockwise_cov(const T *data,const T * data2,const std::vector<size_t> &dims, const std::vector<size_t> &starts,const size_t &blocksize,const double & mean=0,const double & mean2=0){
            size_t N=dims.size();

            if(N==2){
                size_t dimy=dims[1],startx=starts[0],starty=starts[1],element_num=blocksize*blocksize;
                double covsum=0;
                
                for(size_t i=startx;i<startx+blocksize;i++){
                    for(size_t j=starty;j<starty+blocksize;j++){
                        size_t cur_idx=i*dimy+j;
                        T value=data[cur_idx],value2=data2[cur_idx];
                        covsum+=(value-mean)*(value2-mean2);

                    }

                }
                return covsum/element_num;
            }

            

            else if(N==3){
                size_t dimy=dims[1],dimz=dims[2],startx=starts[0],starty=starts[1],startz=starts[2],element_num=blocksize*blocksize*blocksize;
                size_t dimyz=dimy*dimz;
                double covsum=0;
                for(size_t i=startx;i<startx+blocksize;i++){
                    for(size_t j=starty;j<starty+blocksize;j++){
                        for(size_t k=startz;k<startz+blocksize;k++){
                            size_t cur_idx=i*dimyz+j*dimz+k;
                            T value=data[cur_idx],value2=data2[cur_idx];
                            covsum+=value*value2;
                        }

                    }

                }
                return covsum/element_num-mean*mean2;
            }
            else{
                return 0;
            }
        }

        template <class T>
        double calc_blockwise_autocorrelation(const T *data,const T * data2,const std::vector<size_t> &dims, const std::vector<size_t> &starts,const size_t &blocksize){
             size_t N=dims.size();
             size_t element_num;
             std::vector<T>diffs;
             if(N==2){
                size_t dimx=dims[0],dimy=dims[1],startx=starts[0],starty=starts[1];
                element_num=blocksize*blocksize;
                diffs.resize(element_num,0);
                
                size_t idx=0;
                for(size_t i=startx;i<startx+blocksize;i++){
                    for(size_t j=starty;j<starty+blocksize;j++){
                        size_t cur_idx=i*dimy+j;
                        T value=data[cur_idx],value2=data2[cur_idx];
                        diffs[idx++]=value-value2;

                    }

                }
                
            }

            

            else if(N==3){
                size_t dimx=dims[0],dimy=dims[1],dimz=dims[2],startx=starts[0],starty=starts[1],startz=starts[2];
                size_t dimyz=dimy*dimz;
                element_num=blocksize*blocksize*blocksize;
                
                diffs.resize(element_num,0);
                size_t idx=0;
                for(size_t i=startx;i<startx+blocksize;i++){
                    for(size_t j=starty;j<starty+blocksize;j++){
                        for(size_t k=startz;k<startz+blocksize;k++){
                            size_t cur_idx=i*dimyz+j*dimz+k;
                            T value=data[cur_idx],value2=data2[cur_idx];
                            diffs[idx++]=value-value2;
                        }

                    }

                }
                
            }
            else{
                return 0;
            }
            double avg=0;
            for(size_t i=0;i<element_num;i++){
                avg+=diffs[i];
            }
            avg/=element_num;
            double cov=0;
            for(size_t i=0;i<element_num;i++){
                T value=diffs[i]-avg;
                cov+=value*value;
            }
            cov/=element_num;
            if (cov==0){
                return 1;
            }
            double sum=0;
            for (size_t i=0;i<element_num-1;i++){
                sum+=(diffs[i]-avg)*(diffs[i+1]-avg);
            }
            return sum/(element_num-1)/cov;





        }

        template <class T>
        double calc_autocorrelation(const T *data, const T * data2,const size_t &element_num){
             
             
            std::vector<T>diffs(element_num,0);

             
            double avg=0;
            for(size_t i=0;i<element_num;i++){
                T value=data[i]-data2[i];
                avg+=value;
                diffs[i]=value;
            }
            avg/=element_num;
            double cov=0;
            for(size_t i=0;i<element_num;i++){
                T value=diffs[i]-avg;
                cov+=value*value;
            }
            cov/=element_num;
            if (cov==0){
                return 1.0;
            }
            double sum=0;
            for (size_t i=0;i<element_num-1;i++){
                sum+=(diffs[i]-avg)*(diffs[i+1]-avg);
            }
            return sum/(element_num-1)/cov;





        }
   }


    
}


#endif //SZ_INTERPOLATORS_HPP
