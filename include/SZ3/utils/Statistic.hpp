//
// Created by Kai Zhao on 4/20/20.
//

#ifndef SZ_STATISTIC_HPP
#define SZ_STATISTIC_HPP

#include "Config.hpp"

namespace SZ {
    template<class T>
    T data_range(T *data, size_t num) {
        T max = data[0];
        T min = data[0];
        for (size_t i = 1; i < num; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }
        return max - min;
    }

    int factorial(int n) {
        return (n == 0) || (n == 1) ? 1 : n * factorial(n - 1);
    }

    double computeABSErrBoundFromPSNR(double psnr, double threshold, double value_range) {
        double v1 = psnr + 10 * log10(1 - 2.0 / 3.0 * threshold);
        double v2 = v1 / (-20);
        double v3 = pow(10, v2);
        return value_range * v3;
    }

    template<class T>
    void calAbsErrorBound(SZ::Config &conf, T *data) {
        if (conf.errorBoundMode != EB_ABS) {
            if (conf.errorBoundMode == EB_REL) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = conf.relErrorBound * SZ::data_range(data, conf.num);
            } else if (conf.errorBoundMode == EB_PSNR) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = computeABSErrBoundFromPSNR(conf.psnrErrorBound, 0.99, SZ::data_range(data, conf.num));
            } else if (conf.errorBoundMode == EB_L2NORM) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = sqrt(3.0 / conf.num) * conf.l2normErrorBound;
            } else if (conf.errorBoundMode == EB_ABS_AND_REL) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = std::min(conf.absErrorBound, conf.relErrorBound * SZ::data_range(data, conf.num));
            } else if (conf.errorBoundMode == EB_ABS_OR_REL) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = std::max(conf.absErrorBound, conf.relErrorBound * SZ::data_range(data, conf.num));
            } else {
                printf("Error, error bound mode not supported\n");
                exit(0);
            }
        }
    }

    template<typename Type>
    double autocorrelation1DLag1(const Type *data, size_t numOfElem, Type avg) {
        double cov = 0;
        for (size_t i = 0; i < numOfElem; i++) {
            cov += (data[i] - avg) * (data[i] - avg);
        }
        cov = cov / numOfElem;

        if (cov == 0) {
            return 0;
        } else {
            int delta = 1;
            double sum = 0;

            for (size_t i = 0; i < numOfElem - delta; i++) {
                sum += (data[i] - avg) * (data[i + delta] - avg);
            }
            return sum / (numOfElem - delta) / cov;
        }
    }

    // auxilliary functions for evaluating regional average
    template <class T>
    std::vector<T> compute_average(T const * data, uint32_t n1, uint32_t n2, uint32_t n3, int block_size){
        uint32_t dim0_offset = n2 * n3;
        uint32_t dim1_offset = n3;
        uint32_t num_block_1 = (n1 - 1) / block_size + 1;
        uint32_t num_block_2 = (n2 - 1) / block_size + 1;
        uint32_t num_block_3 = (n3 - 1) / block_size + 1;
        std::vector<T> aggregated = std::vector<T>();
        uint32_t index = 0;
        T const * data_x_pos = data;
        for(int i=0; i<num_block_1; i++){
            int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
            T const * data_y_pos = data_x_pos;
            for(int j=0; j<num_block_2; j++){
                int size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
                T const * data_z_pos = data_y_pos;
                for(int k=0; k<num_block_3; k++){
                    int size_3 = (k == num_block_3 - 1) ? n3 - k * block_size : block_size;
                    T const * cur_data_pos = data_z_pos;
                    int n_block_elements = size_1 * size_2 * size_3;
                    double sum = 0;
                    for(int ii=0; ii<size_1; ii++){
                        for(int jj=0; jj<size_2; jj++){
                            for(int kk=0; kk<size_3; kk++){
                                sum += *cur_data_pos;
                                cur_data_pos ++;
                            }
                            cur_data_pos += dim1_offset - size_3;
                        }
                        cur_data_pos += dim0_offset - size_2 * dim1_offset;
                    }
                    aggregated.push_back(sum / n_block_elements);
                    data_z_pos += size_3;
                }
                data_y_pos += dim1_offset * size_2;
            }
            data_x_pos += dim0_offset * size_1;
        }    
        return aggregated;
    }

    template<class T>
    T evaluate_L_inf(T const * data, T const * dec_data, uint32_t num_elements, bool normalized=true, bool verbose=false){
        T L_inf_error = 0;
        T L_inf_data = 0;
        for(int i=0; i<num_elements; i++){
            if(L_inf_data < fabs(data[i])) L_inf_data = fabs(data[i]);
            T error = data[i] - dec_data[i];
            if(L_inf_error < fabs(error)) L_inf_error = fabs(error);
        }
        if(verbose){
            std::cout << "Max absolute value = " << L_inf_data << std::endl;
            std::cout << "L_inf error = " << L_inf_error << std::endl;
            std::cout << "L_inf error (normalized) = " << L_inf_error / L_inf_data << std::endl;
        }
        return normalized ? L_inf_error / L_inf_data : L_inf_error;
    }

    template<class T>
    void evaluate_average(T const * data, T const * dec_data, double value_range, uint32_t n1, uint32_t n2, uint32_t n3, int block_size = 1){
        if(block_size == 0){
            double average = 0;
            double average_dec = 0;
            for(int i=0; i<n1 * n2 * n3; i++){
                average += data[i];
                average_dec += dec_data[i]; 
            }
            std::cout << "Error in global average = " << fabs(average - average_dec)/(n1*n2*n3) << std::endl;
        }
        else{
            auto average = compute_average(data, n1, n2, n3, block_size);
            auto average_dec = compute_average(dec_data, n1, n2, n3, block_size);
            auto error = evaluate_L_inf(average.data(), average_dec.data(), average.size(), false, false);
            std::cout << "L^infinity error of average with block size " << block_size << " = " << error << ", relative error = " << error * 1.0 / value_range << std::endl;
            // std::cout << "L^2 error of average with block size " << block_size << " = " << evaluate_L2(average.data(), average_dec.data(), average.size(), true, false) << std::endl;
        }
    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements, double &psnr, double &nrmse) {
        size_t i = 0;
        double Max = ori_data[0];
        double Min = ori_data[0];
        double diffMax = fabs(data[0] - ori_data[0]);
        double diff_sum = 0;
        double maxpw_relerr = 0;
        double sum1 = 0, sum2 = 0, l2sum = 0;
        for (i = 0; i < num_elements; i++) {
            sum1 += ori_data[i];
            sum2 += data[i];
            l2sum += data[i] * data[i];
        }
        double mean1 = sum1 / num_elements;
        double mean2 = sum2 / num_elements;

        double sum3 = 0, sum4 = 0;
        double sum = 0, prodSum = 0, relerr = 0;

        double *diff = (double *) malloc(num_elements * sizeof(double));

        for (i = 0; i < num_elements; i++) {
            diff[i] = data[i] - ori_data[i];
            diff_sum += data[i] - ori_data[i];
            if (Max < ori_data[i]) Max = ori_data[i];
            if (Min > ori_data[i]) Min = ori_data[i];
            double err = fabs(data[i] - ori_data[i]);
            if (ori_data[i] != 0) {
                relerr = err / fabs(ori_data[i]);
                if (maxpw_relerr < relerr)
                    maxpw_relerr = relerr;
            }

            if (diffMax < err)
                diffMax = err;
            prodSum += (ori_data[i] - mean1) * (data[i] - mean2);
            sum3 += (ori_data[i] - mean1) * (ori_data[i] - mean1);
            sum4 += (data[i] - mean2) * (data[i] - mean2);
            sum += err * err;
        }
        double std1 = sqrt(sum3 / num_elements);
        double std2 = sqrt(sum4 / num_elements);
        double ee = prodSum / num_elements;
        double acEff = ee / std1 / std2;

        double mse = sum / num_elements;
        double range = Max - Min;
        psnr = 20 * log10(range) - 10 * log10(mse);
        nrmse = sqrt(mse) / range;

        double normErr = sqrt(sum);
        double normErr_norm = normErr / sqrt(l2sum);

        printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf("Max absolute error = %.2G\n", diffMax);
        printf("Max relative error = %.2G\n", diffMax / (Max - Min));
        printf("Max pw relative error = %.2G\n", maxpw_relerr);
        printf("PSNR = %f, NRMSE= %.10G\n", psnr, nrmse);
        printf("normError = %f, normErr_norm = %f\n", normErr, normErr_norm);
        printf("acEff=%f\n", acEff);
//        printf("errAutoCorr=%.10f\n", autocorrelation1DLag1<double>(diff, num_elements, diff_sum / num_elements));
        free(diff);
    }

    template<typename Type>
    void verifyQoI(Type *ori_data, Type *data, std::vector<size_t> dims, int blockSize=1) {
        size_t num_elements = 1;
        for(const auto d:dims){
            num_elements *= d;
        }
        double psnr = 0;
        double nrmse = 0;
        verify(ori_data, data, num_elements, psnr, nrmse);
        Type max = data[0];
        Type min = data[0];
        for (size_t i = 1; i < num_elements; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }

        double max_abs_val = std::max(fabs(max), fabs(min));
        double max_abs_val_sq = max_abs_val * max_abs_val;
        double max_x_square_diff = 0;
        double max_log_diff = 0;
        for(int i=0; i<num_elements; i++){
            double x_square_diff = fabs(ori_data[i] * ori_data[i] - data[i] * data[i]);
            if(x_square_diff > max_x_square_diff) max_x_square_diff = x_square_diff;
            // if(x_square_diff / max_abs_val_sq > 1e-5){
            //     std::cout << i << ": ori = " << ori_data[i] << ", dec = " << data[i] << ", err = " << x_square_diff / max_abs_val_sq << std::endl;
            // }
            if(ori_data[i] == 0){
                // for log x only
                // if(data[i] != 0){
                //     std::cout << i << ": dec_data does not equal to 0\n";
                // }
            }
            else{
                if(ori_data[i] == data[i]) continue;
                double log2 = log(2);
                double log_diff = fabs(log(fabs(ori_data[i]))/log2 - log(fabs(data[i]))/log2);
                if(log_diff > max_log_diff) max_log_diff = log_diff;  
                // if(log_diff > 1){
                //     std::cout << i << ": " << ori_data[i] << " " << data[i] << std::endl;
                // }              
            }
        }

        printf("QoI error info:\n");
        printf("Max x^2 error = %.6G, relative x^2 error = %.6G\n", max_x_square_diff, max_x_square_diff / max_abs_val_sq);
        printf("Max log error = %.6G\n", max_log_diff);
        if(dims.size() == 3) evaluate_average(ori_data, data, max - min, dims[0], dims[1], dims[2], blockSize);

    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements) {
        double psnr, nrmse;
        verify(ori_data, data, num_elements, psnr, nrmse);
    }
};


#endif
