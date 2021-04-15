#ifndef _meta_prediction_hpp
#define _meta_prediction_hpp

#include "meta_def.hpp"

namespace META {


// quantize with decompression data (Lorenzo)
    template<typename T>
    inline int
    quantize(float pred, T cur_data, double precision, double recip_precision, int capacity, int intv_radius,
             T *&unpredictable_data_pos, T *decompressed) {
        double diff = cur_data - pred;
        double quant_diff = fabs(diff) * recip_precision + 1;
        if (quant_diff < capacity) {
            quant_diff = (diff > 0) ? quant_diff : -quant_diff;
            int quant_index = (int) (quant_diff / 2) + intv_radius;
            T decompressed_data = pred + 2 * (quant_index - intv_radius) * precision;
            *decompressed = decompressed_data;
            if (fabs(decompressed_data - cur_data) <= precision) return quant_index;
        }
        *decompressed = cur_data;
        *(unpredictable_data_pos++) = cur_data;
        return 0;
    }

// return quantization index, no decompression data (regression)
    template<typename T>
    inline int
    quantize(float pred, T cur_data, double precision, double recip_precision, int capacity, int intv_radius,
             T *&unpredictable_data_pos) {
        double diff = cur_data - pred;
        double quant_diff = fabs(diff) * recip_precision + 1;
        if (quant_diff < capacity) {
            quant_diff = (diff > 0) ? quant_diff : -quant_diff;
            int quant_index = (int) (quant_diff / 2) + intv_radius;
            T decompressed_data = pred + 2 * (quant_index - intv_radius) * precision;
            if (fabs(decompressed_data - cur_data) <= precision) return quant_index;
        }
        *(unpredictable_data_pos++) = cur_data;
        return 0;
    }

// de-quantization, for regression
    template<typename T>
    inline T
    recover(float pred, double precision, int type_val, int intv_radius, const T *&unpredictable_data_pos) {
        if (type_val == 0) {
            return *(unpredictable_data_pos++);
        } else {
            return pred + 2 * (type_val - intv_radius) * precision;
        }
    }

// de-quantization, for use_lorenzo
    template<typename T>
    inline T
    recover(const meanInfo <T> &mean_info, float pred, double precision, int type_val, int intv_radius,
            const T *&unpredictable_data_pos) {
        if (type_val == 0) {
            return *(unpredictable_data_pos++);
        } else {
            if ((type_val == 1) && (mean_info.use_mean)) return mean_info.mean;
            return pred + 2 * (type_val - intv_radius) * precision;
        }
    }

    template<typename T>
    inline void
    compress_regression_coefficient_3d(const int coeff_num, const T *precisions, const T *recip_precisions, float *reg_params_pos,
                                       int *reg_params_type_pos, float *&reg_unpredictable_data_pos) {
        float *prev_reg_params = reg_params_pos - coeff_num;
        for (int i = 0; i < coeff_num; i++) {
            *(reg_params_type_pos++) = quantize(*prev_reg_params, *reg_params_pos, precisions[i], recip_precisions[i],
                                                RegCoeffCapacity, RegCoeffRadius, reg_unpredictable_data_pos, reg_params_pos);
            prev_reg_params++, reg_params_pos++;
        }
    }

}
#endif