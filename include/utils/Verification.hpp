//
// Created by Kai Zhao on 4/20/20.
//

#ifndef SZ_VERIFICATION_HPP
#define SZ_VERIFICATION_HPP


namespace SZ {

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

    template<typename Type>
    void verify(const char *ori_datafile, const char *dec_datafile, size_t num_elements, double &psnr, double &nrmse) {
        std::ifstream fori(ori_datafile, std::ios::binary);
        std::ifstream fdec(dec_datafile, std::ios::binary);
        if (!fori || !fdec) {
            std::cout << " Error, Couldn't find the file" << "\n";
            exit(0);
        }

        size_t buffer_size = 1000 * 1000, start = 0;
        std::vector<Type> data(buffer_size), ori_data(buffer_size);

        double max = std::numeric_limits<Type>::min();
        double min = std::numeric_limits<Type>::max();
        double max_err = 0;
        double maxpw_relerr = 0;
        double l2_err = 0;

        while (start < num_elements) {
            size_t size = num_elements - start > buffer_size ? buffer_size : num_elements - start;
            fori.read(reinterpret_cast<char *>(&ori_data[0]), size * sizeof(Type));
            fdec.read(reinterpret_cast<char *>(&data[0]), size * sizeof(Type));

            for (size_t i = 0; i < size; i++) {
                if (max < ori_data[i]) max = ori_data[i];
                if (min > ori_data[i]) min = ori_data[i];
                double err = fabs(data[i] - ori_data[i]);
                if (ori_data[i] != 0) {
                    double relerr = err / fabs(ori_data[i]);
                    if (maxpw_relerr < relerr)
                        maxpw_relerr = relerr;
                }

                if (max_err < err)
                    max_err = err;
                l2_err += err * err;
            }
            start += size;
        }

        double mse = l2_err / num_elements;
        double range = max - min;
        psnr = 20 * log10(range) - 10 * log10(mse);
        nrmse = sqrt(mse) / range;

        printf("min=%.20G, max=%.20G, range=%.20G\n", min, max, range);
        printf("max absolute error = %.2G\n", max_err);
        printf("max relative error = %.2G\n", max_err / (max - min));
        printf("max pw relative error = %.2G\n", maxpw_relerr);
        printf("PSNR = %f, NRMSE= %.10G\n", psnr, nrmse);
    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements, double &psnr, double &nrmse) {
        size_t i = 0;
        double Max = ori_data[0];
        double Min = ori_data[0];
        double diffMax = fabs(data[0] - ori_data[0]);
        double diff_sum = 0;
        double maxpw_relerr = fabs((data[0] - ori_data[0]) / ori_data[0]);
        double sum1 = 0, sum2 = 0;
        for (i = 0; i < num_elements; i++) {
            sum1 += ori_data[i];
            sum2 += data[i];
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

        printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf("Max absolute error = %.2G\n", diffMax);
        printf("Max relative error = %.2G\n", diffMax / (Max - Min));
        printf("Max pw relative error = %.2G\n", maxpw_relerr);
        printf("PSNR = %f, NRMSE= %.10G\n", psnr, nrmse);
        printf("acEff=%f\n", acEff);
        printf("errAutoCorr=%.10f\n", autocorrelation1DLag1<double>(diff, num_elements, diff_sum / num_elements));
        free(diff);
    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements) {
        double psnr, nrmse;
        verify(ori_data, data, num_elements, psnr, nrmse);
    }
};


#endif //SZ_VERIFICATION_HPP
