#ifndef _meta_autotuning_3d_hpp
#define _meta_autotuning_3d_hpp

#include <cstddef>
#include "meta_prediction.hpp"
#include "meta_quantization.hpp"
#include "meta_encode.hpp"
#include "meta_lossless.hpp"
#include "meta_optimize_quant_intervals.hpp"
#include "meta_def.hpp"
#include "meta_utils.hpp"
#include <random>
#include <list>
#include <utils/Verification.hpp>

namespace META {
    using namespace std;

    template<typename T>
    inline void
    meta_block_error_estimation_3d(const T *data_pos, const float *reg_params_pos,
                                   const meanInfo<T> &mean_info, int x, int y, int z, size_t dim0_offset, size_t dim1_offset,
                                   T precision, double &err_lorenzo, double &err_lorenzo_2layer, double &err_reg,
                                   const int pred_dim,
                                   const bool use_lorenzo, const bool use_lorenzo_2layer,
                                   const bool use_regression) {
        T noise = 0;
        T noise_2layer = 0;
        const T *cur_data_pos = data_pos + x * dim0_offset + y * dim1_offset + z;
        T cur_data = *cur_data_pos;
        if (use_regression) {
            err_reg += fabs(cur_data - regression_predict_3d<T>(reg_params_pos, x, y, z));
        }
        double lorenzo_predict = 0;
        double lorenzo_2layer_predict = 0;
        if (pred_dim == 3) {
            if (use_lorenzo_2layer) {
                lorenzo_2layer_predict = lorenzo_predict_3d_2layer(cur_data_pos, dim0_offset, dim1_offset);
                noise_2layer = Lorenze2LayerNoise3d * precision;
            }
            if (use_lorenzo) {
                lorenzo_predict = lorenzo_predict_3d(cur_data_pos, dim0_offset, dim1_offset);
                noise = LorenzeNoise3d * precision;
            }
        } else if (pred_dim == 2) {
            if (use_lorenzo_2layer) {
                lorenzo_2layer_predict = lorenzo_predict_2d_2layer(cur_data_pos, dim0_offset, dim1_offset);
                noise_2layer = Lorenze2LayerNoise2d * precision;

            }
            if (use_lorenzo) {
                lorenzo_predict = lorenzo_predict_2d(cur_data_pos, dim0_offset, dim1_offset);
                noise = LorenzeNoise2d * precision;
            }
        } else {
            if (use_lorenzo_2layer) {
                lorenzo_2layer_predict = lorenzo_predict_1d_2layer(cur_data_pos, dim0_offset);
                noise_2layer = Lorenze2LayerNoise1d * precision;

            }
            if (use_lorenzo) {
                lorenzo_predict = lorenzo_predict_1d(cur_data_pos, dim0_offset);
                noise = LorenzeNoise1d * precision;
            }
        }
        err_lorenzo += mean_info.use_mean ? MIN(fabs(cur_data - mean_info.mean), fabs(cur_data - lorenzo_predict) + noise) :
                       fabs(cur_data - lorenzo_predict) + noise;
        err_lorenzo_2layer += mean_info.use_mean ? MIN(fabs(cur_data - mean_info.mean),
                                                       fabs(cur_data - lorenzo_2layer_predict) + noise_2layer) :
                              fabs(cur_data - lorenzo_2layer_predict) + noise_2layer;
    }


    template<typename T>
    inline int
    meta_blockwise_selection_3d(const T *data_pos, const meanInfo<T> &mean_info, size_t dim0_offset, size_t dim1_offset,
                                int min_size,
                                T precision, const float *reg_params_pos, const int pred_dim,
                                const bool use_lorenzo, const bool use_lorenzo_2layer, const bool use_regression) {
        double err_lorenzo = 0;
        double err_lorenzo_2layer = 0;
        double err_reg = 0;
        for (int i = 2; i < min_size - 1; i++) {
            int bmi = min_size - i;
            meta_block_error_estimation_3d(data_pos, reg_params_pos, mean_info, i, i, i, dim0_offset, dim1_offset,
                                           precision, err_lorenzo, err_lorenzo_2layer, err_reg, pred_dim, use_lorenzo,
                                           use_lorenzo_2layer, use_regression);
            meta_block_error_estimation_3d(data_pos, reg_params_pos, mean_info, i, i, bmi, dim0_offset,
                                           dim1_offset, precision, err_lorenzo, err_lorenzo_2layer, err_reg, pred_dim,
                                           use_lorenzo, use_lorenzo_2layer, use_regression);
            meta_block_error_estimation_3d(data_pos, reg_params_pos, mean_info, i, bmi, i, dim0_offset,
                                           dim1_offset, precision, err_lorenzo, err_lorenzo_2layer, err_reg, pred_dim,
                                           use_lorenzo, use_lorenzo_2layer, use_regression);
            meta_block_error_estimation_3d(data_pos, reg_params_pos, mean_info, i, bmi, bmi, dim0_offset,
                                           dim1_offset, precision, err_lorenzo, err_lorenzo_2layer, err_reg, pred_dim,
                                           use_lorenzo, use_lorenzo_2layer, use_regression);
        }
        if (min_size > 3) {
            meta_block_error_estimation_3d(data_pos, reg_params_pos, mean_info, min_size - 1, min_size - 1, min_size - 1,
                                           dim0_offset, dim1_offset,
                                           precision, err_lorenzo, err_lorenzo_2layer, err_reg, pred_dim, use_lorenzo,
                                           use_lorenzo_2layer, use_regression);
        }

        if (use_regression && (!use_lorenzo || err_reg <= err_lorenzo)
            && (!use_lorenzo_2layer || err_reg < err_lorenzo_2layer)) {
            return SELECTOR_REGRESSION;
        } else if (use_lorenzo_2layer && (!use_lorenzo || err_lorenzo_2layer <= err_lorenzo)
                   && (!use_regression || err_lorenzo_2layer <= err_reg)) {
            return SELECTOR_LORENZO_2LAYER;
        } else {
            return SELECTOR_LORENZO;
        }
    }


    template<typename T>
    unsigned char *
    meta_compress_3d(const T *data, size_t r1, size_t r2, size_t r3, double precision, size_t &compressed_size,
                     const meta_params &params, meta_compress_info &compress_info) {
//    ori_data_sp_float = (float *) data;
        DSize_3d size(r1, r2, r3, params.block_size);
        int capacity = 0; // num of quant intervals
        meanInfo<T> mean_info = optimize_quant_invl_3d(data, r1, r2, r3, precision, capacity);
        if (params.capacity > 0) {
            capacity = params.capacity;
        }
        int intv_radius = (capacity >> 1);
        int *type = (int *) malloc(size.num_elements * sizeof(int));
        int *indicator = (int *) malloc(size.num_blocks * sizeof(int));
        int *reg_params_type = (int *) malloc(RegCoeffNum3d * size.num_blocks * sizeof(int));
        float *reg_unpredictable_data = (float *) malloc(RegCoeffNum3d * size.num_blocks * sizeof(float));
        float *reg_unpredictable_data_pos = reg_unpredictable_data;

        // prepare unpred buffer for vectorization
        int est_unpred_count_per_index = size.num_blocks * size.block_size * 1;
        // if(!params.block_independant) est_unpred_count_per_index /= 20;
        T *unpred_data_buffer = (T *) malloc(size.block_size * size.block_size * est_unpred_count_per_index * sizeof(T));
        int *unpred_count_buffer = (int *) malloc(size.block_size * size.block_size * sizeof(int));
        memset(unpred_count_buffer, 0, size.block_size * size.block_size * sizeof(int));
//        T precision_t = (T) precision;
        size_t reg_count = 0;
        size_t lorenzo_count = 0;
        size_t lorenzo_2layer_count = 0;

        int *type_pos = type;
        int *indicator_pos = indicator;

        float *reg_params = (float *) malloc(RegCoeffNum3d * (size.num_blocks + 1) * sizeof(float));
        for (int i = 0; i < RegCoeffNum3d; i++) {
            reg_params[i] = 0;
        }
        float *reg_params_pos = reg_params + RegCoeffNum3d;
        int *reg_params_type_pos = reg_params_type;


        T reg_precisions[RegCoeffNum3d];
        T reg_recip_precisions[RegCoeffNum3d];
        for (int i = 0; i < RegCoeffNum3d - 1; i++) {
            reg_precisions[i] = params.regression_param_eb_linear;
            reg_recip_precisions[i] = 1.0 / reg_precisions[i];
        }
        reg_precisions[RegCoeffNum3d - 1] = params.regression_param_eb_independent;
        reg_recip_precisions[RegCoeffNum3d - 1] = 1.0 / reg_precisions[RegCoeffNum3d - 1];

        // maintain a buffer of (block_size+1)*(r2+1)*(r3+1)
        // 2-layer use_lorenzo
        size_t buffer_dim0_offset = (size.d2 + params.lorenzo_padding_layer) * (size.d3 + params.lorenzo_padding_layer);
        size_t buffer_dim1_offset = size.d3 + params.lorenzo_padding_layer;
        T *pred_buffer = (T *) malloc(
                (size.block_size + params.lorenzo_padding_layer) * (size.d2 + params.lorenzo_padding_layer) *
                (size.d3 + params.lorenzo_padding_layer) * sizeof(T));
        memset(pred_buffer, 0, (size.block_size + params.lorenzo_padding_layer) * (size.d2 + params.lorenzo_padding_layer) *
                               (size.d3 + params.lorenzo_padding_layer) * sizeof(T));
        int capacity_lorenzo = mean_info.use_mean ? capacity - 2 : capacity;
        auto *lorenzo_pred_and_quant = compress_lorenzo_3d_predict<T>;
        if (params.prediction_dim == 2) lorenzo_pred_and_quant = compress_lorenzo_3d_as2d_predict<T>;
        else if (params.prediction_dim == 1) lorenzo_pred_and_quant = compress_lorenzo_3d_as1d_predict<T>;
        T recip_precision = (T) 1.0 / precision;

        const T *x_data_pos = data;
        for (size_t i = 0; i < size.num_x; i++) {
            const T *y_data_pos = x_data_pos;
            T *pred_buffer_pos = pred_buffer;
            for (size_t j = 0; j < size.num_y; j++) {
                const T *z_data_pos = y_data_pos;
                for (size_t k = 0; k < size.num_z; k++) {
                    int size_x = ((i + 1) * size.block_size < size.d1) ? size.block_size : size.d1 - i * size.block_size;
                    int size_y = ((j + 1) * size.block_size < size.d2) ? size.block_size : size.d2 - j * size.block_size;
                    int size_z = ((k + 1) * size.block_size < size.d3) ? size.block_size : size.d3 - k * size.block_size;
                    int min_size = MIN(size_x, size_y);
                    min_size = MIN(min_size, size_z);

                    bool enable_regression = params.use_regression_linear && min_size >= 2;
//                bool enable_regression = params.use_regression_linear && min_size >= 1;

                    if (enable_regression) {
                        compute_regression_coeffcients_3d(z_data_pos, size_x, size_y, size_z, size.dim0_offset,
                                                          size.dim1_offset,
                                                          reg_params_pos);
                    }

                    int selection_result = meta_blockwise_selection_3d<T>(z_data_pos, mean_info, size.dim0_offset,
                                                                          size.dim1_offset,
                                                                          min_size, precision, reg_params_pos,
                                                                          params.prediction_dim,
                                                                          params.use_lorenzo, params.use_lorenzo_2layer,
                                                                          enable_regression);
                    *indicator_pos = selection_result;
                    if (selection_result == SELECTOR_REGRESSION) {
                        // regression
                        compress_regression_coefficient_3d(RegCoeffNum3d, reg_precisions, reg_recip_precisions,
                                                           reg_params_pos,
                                                           reg_params_type_pos,
                                                           reg_unpredictable_data_pos);
                        compress_regression_3d_predict<T>(z_data_pos, reg_params_pos, pred_buffer_pos, precision,
                                                          recip_precision, capacity, intv_radius,
                                                          size_x, size_y, size_z, buffer_dim0_offset,
                                                          buffer_dim1_offset, size.dim0_offset, size.dim1_offset,
                                                          type_pos, unpred_count_buffer, unpred_data_buffer,
                                                          est_unpred_count_per_index,
                                                          params.lorenzo_padding_layer);
                        reg_count++;
                        reg_params_pos += RegCoeffNum3d;
                        reg_params_type_pos += RegCoeffNum3d;
                    } else {
                        // Lorenzo
                        lorenzo_pred_and_quant(mean_info, z_data_pos, pred_buffer_pos, precision, recip_precision,
                                               capacity_lorenzo,
                                               intv_radius,
                                               size_x, size_y, size_z, buffer_dim0_offset, buffer_dim1_offset,
                                               size.dim0_offset,
                                               size.dim1_offset, type_pos, unpred_count_buffer, unpred_data_buffer,
                                               est_unpred_count_per_index,
                                               params.lorenzo_padding_layer, (selection_result == SELECTOR_LORENZO_2LAYER));
                        if (selection_result == SELECTOR_LORENZO_2LAYER) {
                            lorenzo_2layer_count++;
                        } else {
                            lorenzo_count++;
                        }
                    }
                    pred_buffer_pos += size.block_size;
                    indicator_pos++;
                    z_data_pos += size_z;
                }
                y_data_pos += size.block_size * size.dim1_offset;
                pred_buffer_pos += size.block_size * buffer_dim1_offset - size.block_size * size.num_z;
            }
            // copy bottom of buffer to top of buffer
            memcpy(pred_buffer, pred_buffer + size.block_size * buffer_dim0_offset,
                   params.lorenzo_padding_layer * buffer_dim0_offset * sizeof(T));
            x_data_pos += size.block_size * size.dim0_offset;
        }
        free(pred_buffer);
        free(reg_params);

//    printf("block %ld; lorenzo %ld, lorenzo_2layer %ld, regression %ld, poly regression %ld\n", size.num_blocks,
//        lorenzo_count, lorenzo_2layer_count, reg_count);
        compress_info.lorenzo_count = lorenzo_count;
        compress_info.lorenzo2_count = lorenzo_2layer_count;
        compress_info.regression_count = reg_count;
        compress_info.block_count = size.num_blocks;

        size_t est_size = size.num_elements * sizeof(T) * 2;

        unsigned char *compressed = (unsigned char *) malloc(est_size);
        unsigned char *compressed_pos = compressed;
        write_variable_to_dst(compressed_pos, params);
        write_variable_to_dst(compressed_pos, precision);
        write_variable_to_dst(compressed_pos, intv_radius);
        write_variable_to_dst(compressed_pos, mean_info);
        write_variable_to_dst(compressed_pos, reg_count);
        write_array_to_dst(compressed_pos, unpred_count_buffer, size.block_size * size.block_size);
        T *unpred_data_buffer_pos = unpred_data_buffer;
        for (int i = 0; i < size.block_size; i++) {
            for (int j = 0; j < size.block_size; j++) {
                write_array_to_dst(compressed_pos, unpred_data_buffer_pos, unpred_count_buffer[i * size.block_size + j]);
                unpred_data_buffer_pos += est_unpred_count_per_index;
            }
        }
        free(unpred_data_buffer);
        free(unpred_count_buffer);


        Huffman_encode_tree_and_data(SELECTOR_RADIUS, indicator, size.num_blocks, compressed_pos);
//	convertIntArray2ByteArray_fast_1b_to_result_sz(indicator, size.num_blocks, compressed_pos);

        if (reg_count) {
            encode_regression_coefficients(reg_params_type, reg_unpredictable_data, RegCoeffNum3d * reg_count,
                                           reg_unpredictable_data_pos - reg_unpredictable_data, compressed_pos);
        }

        Huffman_encode_tree_and_data(2 * capacity, type, size.num_elements, compressed_pos);

        compressed_size = compressed_pos - compressed;
        free(indicator);
        free(reg_params_type);
        free(reg_unpredictable_data);
        free(type);
        return compressed;
    }


    template<typename T>
    unsigned char *
    meta_compress_3d_sampling(const T *data, size_t r1, size_t r2, size_t r3, double precision, size_t &compressed_size,
                              const meta_params &params, meta_compress_info &compress_info) {
//    ori_data_sp_float = (float *) data;
        DSize_3d size(r1, r2, r3, params.block_size);
        int capacity = params.capacity;; // num of quant intervals
//    meanInfo<T> mean_info = optimize_quant_invl_3d(data, r1, r2, r3, precision, capacity);
        meanInfo<T> mean_info;
        {
            size_t r1 = size.d1 - params.lorenzo_padding_layer;
            size_t r2 = size.d2 - params.lorenzo_padding_layer;
            size_t r3 = size.d3 - params.lorenzo_padding_layer;
            auto block_size = size.block_size;
            size_t num_x = r1 / block_size;
            size_t num_y = r2 / block_size;
            size_t num_z = r3 / block_size;
            size.sample_distance = 0;
            double sample_ratio = 1;
            size.sample_nx = num_x, size.sample_ny = num_y, size.sample_nz = num_z;
            while (sample_ratio > params.sample_ratio) {
                size.sample_distance++;
                size.sample_nx = r1 / (size.sample_distance + block_size);
                size.sample_ny = r2 / (size.sample_distance + block_size);
                size.sample_nz = r3 / (size.sample_distance + block_size);
                sample_ratio = size.sample_nx * size.sample_ny * size.sample_nz * 1.0 / (num_x * num_y * num_z);
                if ((size.sample_nx == 1) && (size.sample_ny == 1) && (size.sample_nz == 1)) break;
            }
            size_t num_sample_blocks = size.sample_nx * size.sample_ny * size.sample_nz;
            sample_ratio = num_sample_blocks * 1.0 / (num_x * num_y * num_z);

            // modify size
            size.num_elements = num_sample_blocks * block_size * block_size * block_size;
            size.num_blocks = num_sample_blocks;
//        printf("sample_ratio = %.4f, num_blocks = %d %d %d\n", sample_ratio, size.sample_nx, size.sample_ny, size.sample_nz);
//        printf("num_elements = %ld, num_blocks = %ld\n", size.num_elements, size.num_blocks);

            compress_info.ori_bytes = size.num_elements * sizeof(float);
        }


        int intv_radius = (capacity >> 1);
        int *type = (int *) malloc(size.num_elements * sizeof(int));
        int *indicator = (int *) malloc(size.num_blocks * sizeof(int));
        int *reg_params_type = (int *) malloc(RegCoeffNum3d * size.num_blocks * sizeof(int));
        float *reg_unpredictable_data = (float *) malloc(RegCoeffNum3d * size.num_blocks * sizeof(float));
        float *reg_unpredictable_data_pos = reg_unpredictable_data;

        // prepare unpred buffer for vectorization
        int est_unpred_count_per_index = size.num_blocks * size.block_size * 1;
        // if(!params.block_independant) est_unpred_count_per_index /= 20;
        T *unpred_data_buffer = (T *) malloc(size.block_size * size.block_size * est_unpred_count_per_index * sizeof(T));
        int *unpred_count_buffer = (int *) malloc(size.block_size * size.block_size * sizeof(int));
        memset(unpred_count_buffer, 0, size.block_size * size.block_size * sizeof(int));
        // predict and quant on KNL
        T precision_t = (T) precision;
        size_t reg_count = 0;

        size_t lorenzo_count = 0;
        size_t lorenzo_2layer_count = 0;

        int *type_pos = type;
        int *indicator_pos = indicator;

        float *reg_params = (float *) malloc(RegCoeffNum3d * (size.num_blocks + 1) * sizeof(float));
        for (int i = 0; i < RegCoeffNum3d; i++) {
            reg_params[i] = 0;
        }
        float *reg_params_pos = reg_params + RegCoeffNum3d;
        int *reg_params_type_pos = reg_params_type;


        T reg_precisions[RegCoeffNum3d];
        T reg_recip_precisions[RegCoeffNum3d];
        for (int i = 0; i < RegCoeffNum3d - 1; i++) {
            reg_precisions[i] = params.regression_param_eb_linear;
            reg_recip_precisions[i] = 1.0 / reg_precisions[i];
        }
        reg_precisions[RegCoeffNum3d - 1] = params.regression_param_eb_independent;
        reg_recip_precisions[RegCoeffNum3d - 1] = 1.0 / reg_precisions[RegCoeffNum3d - 1];

        int capacity_lorenzo = mean_info.use_mean ? capacity - 2 : capacity;
        auto *lorenzo_pred_and_quant = compress_lorenzo_3d_predict<T>;
        if (params.prediction_dim == 2) lorenzo_pred_and_quant = compress_lorenzo_3d_as2d_predict<T>;
        else if (params.prediction_dim == 1) lorenzo_pred_and_quant = compress_lorenzo_3d_as1d_predict<T>;
        T recip_precision = (T) 1.0 / precision;

        // maintain a buffer of (block_size+1)*(r2+1)*(r3+1)
        // 2-layer lorenzo
        size_t buffer_dim0_offset =
                (size.block_size + params.lorenzo_padding_layer) * (size.block_size + params.lorenzo_padding_layer);
        size_t buffer_dim1_offset = size.block_size + params.lorenzo_padding_layer;
        T *pred_buffer = (T *) malloc(
                (size.block_size + params.lorenzo_padding_layer) * (size.block_size + params.lorenzo_padding_layer) *
                (size.block_size + params.lorenzo_padding_layer) * sizeof(T));

        memset(pred_buffer, 0,
               (size.block_size + params.lorenzo_padding_layer) * (size.block_size + params.lorenzo_padding_layer) *
               (size.block_size + params.lorenzo_padding_layer) * sizeof(T));

        T *random_buffer = (T *) malloc(
                (size.block_size + params.lorenzo_padding_layer) * (size.block_size + params.lorenzo_padding_layer) *
                (size.block_size + params.lorenzo_padding_layer) * sizeof(T));

        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(-precision, precision);
        for (int i = 0; i < (size.block_size + params.lorenzo_padding_layer) * (size.block_size + params.lorenzo_padding_layer) *
                            (size.block_size + params.lorenzo_padding_layer); i++) {
            random_buffer[i] = dis(gen);
        }


        {
            size_t r1 = size.d1 - params.lorenzo_padding_layer;
            size_t r2 = size.d2 - params.lorenzo_padding_layer;
            size_t r3 = size.d3 - params.lorenzo_padding_layer;
            const T *data_pos = data;// + layer * r2*r3 + layer * r3 + layer;
            for (size_t ix = 0; ix < size.sample_nx; ix++) {
                for (size_t iy = 0; iy < size.sample_ny; iy++) {
                    for (size_t iz = 0; iz < size.sample_nz; iz++) {
                        // extract data
                        const T *cur_block_data_pos =
                                data_pos + ix * (size.sample_distance + size.block_size) * r2 * r3 +
                                iy * (size.sample_distance + size.block_size) * r3 +
                                iz * (size.sample_distance + size.block_size);

                        //TODO Kai comment this block if to test sampling accuracy without boarder
                        {
                            for (int ii = 0; ii < size.block_size + params.lorenzo_padding_layer; ii++) {
                                for (int jj = 0; jj < size.block_size + params.lorenzo_padding_layer; jj++) {
                                    for (int kk = 0; kk < size.block_size + params.lorenzo_padding_layer; kk++) {
                                        pred_buffer[ii * buffer_dim0_offset + jj * buffer_dim1_offset + kk] =
                                                cur_block_data_pos[ii * size.dim0_offset + jj * size.dim1_offset + kk] +
                                                random_buffer[ii * buffer_dim0_offset + jj * buffer_dim1_offset + kk];
                                    }
                                }
                            }
                        }

                        cur_block_data_pos += params.lorenzo_padding_layer * (size.dim0_offset + size.dim1_offset + 1);

                        bool enable_regression = params.use_regression_linear && size.block_size >= 2;
//                bool enable_regression = params.use_regression_linear && min_size >= 1;

                        if (enable_regression) {
                            compute_regression_coeffcients_3d(cur_block_data_pos, size.block_size, size.block_size,
                                                              size.block_size,
                                                              size.dim0_offset, size.dim1_offset,
                                                              reg_params_pos);
                        }
                        int selection_result = meta_blockwise_selection_3d<T>(cur_block_data_pos, mean_info, size.dim0_offset,
                                                                              size.dim1_offset,
                                                                              size.block_size, precision, reg_params_pos,
                                                                              params.prediction_dim,
                                                                              params.use_lorenzo, params.use_lorenzo_2layer,
                                                                              enable_regression);
                        *indicator_pos = selection_result;
                        if (selection_result == SELECTOR_REGRESSION) {
                            // regression
                            compress_regression_coefficient_3d(RegCoeffNum3d, reg_precisions, reg_recip_precisions,
                                                               reg_params_pos,
                                                               reg_params_type_pos,
                                                               reg_unpredictable_data_pos);
                            compress_regression_3d_predict<T>(cur_block_data_pos, reg_params_pos, pred_buffer,
                                                              precision,
                                                              recip_precision, capacity, intv_radius,
                                                              size.block_size, size.block_size, size.block_size,
                                                              buffer_dim0_offset,
                                                              buffer_dim1_offset, size.dim0_offset, size.dim1_offset,
                                                              type_pos, unpred_count_buffer, unpred_data_buffer,
                                                              est_unpred_count_per_index,
                                                              params.lorenzo_padding_layer);
                            reg_count++;
                            reg_params_pos += RegCoeffNum3d;
                            reg_params_type_pos += RegCoeffNum3d;
                        } else {
                            // Lorenzo
                            lorenzo_pred_and_quant(mean_info, cur_block_data_pos, pred_buffer, precision, recip_precision,
                                                   capacity_lorenzo,
                                                   intv_radius,
                                                   size.block_size, size.block_size, size.block_size, buffer_dim0_offset,
                                                   buffer_dim1_offset,
                                                   size.dim0_offset,
                                                   size.dim1_offset, type_pos, unpred_count_buffer, unpred_data_buffer,
                                                   est_unpred_count_per_index,
                                                   params.lorenzo_padding_layer, (selection_result == SELECTOR_LORENZO_2LAYER));
                            if (selection_result == SELECTOR_LORENZO_2LAYER) {
                                lorenzo_2layer_count++;
                            } else {
                                lorenzo_count++;
                            }
                        }
                        indicator_pos++;
                    }
                }
            }
        }
        free(pred_buffer);
        free(reg_params);
        free(random_buffer);
//    printf("num_type = %ld\n", type_pos - type);
//    printf("block %ld; lorenzo %ld, lorenzo_2layer %ld, regression %ld, poly regression %ld\n", size.num_blocks,
//           lorenzo_count, lorenzo_2layer_count, reg_count, reg_poly_count);
        compress_info.lorenzo_count = lorenzo_count;
        compress_info.lorenzo2_count = lorenzo_2layer_count;
        compress_info.regression_count = reg_count;
        compress_info.block_count = size.num_blocks;

        unsigned char *compressed = NULL;
        // TODO: change to a better estimated size
        size_t est_size = size.num_elements * sizeof(T) * 4;
        compressed = (unsigned char *) malloc(est_size);
        unsigned char *compressed_pos = compressed;
        write_variable_to_dst(compressed_pos, params);
        write_variable_to_dst(compressed_pos, precision);
        write_variable_to_dst(compressed_pos, intv_radius);
        write_variable_to_dst(compressed_pos, mean_info);
        write_variable_to_dst(compressed_pos, reg_count);
        write_array_to_dst(compressed_pos, unpred_count_buffer, size.block_size * size.block_size);
        T *unpred_data_buffer_pos = unpred_data_buffer;
        for (int i = 0; i < size.block_size; i++) {
            for (int j = 0; j < size.block_size; j++) {
                write_array_to_dst(compressed_pos, unpred_data_buffer_pos, unpred_count_buffer[i * size.block_size + j]);
                unpred_data_buffer_pos += est_unpred_count_per_index;
            }
        }
        free(unpred_data_buffer);
        free(unpred_count_buffer);


        Huffman_encode_tree_and_data(SELECTOR_RADIUS, indicator, size.num_blocks, compressed_pos);

//	convertIntArray2ByteArray_fast_1b_to_result_sz(indicator, size.num_blocks, compressed_pos);

        if (reg_count) {
            encode_regression_coefficients(reg_params_type, reg_unpredictable_data, RegCoeffNum3d * reg_count,
                                           reg_unpredictable_data_pos - reg_unpredictable_data, compressed_pos);

        }



//    auto compressed_pos1 = compressed_pos;
        Huffman_encode_tree_and_data(2 * capacity, type, size.num_elements, compressed_pos);
//    printf("after huffman %ld before huffman%ld\n ", compressed_pos - compressed_pos1, size.num_elements * sizeof(T));


        compressed_size = compressed_pos - compressed;
        free(indicator);
        free(reg_params_type);
        free(reg_unpredictable_data);
        free(type);
        return compressed;
    }

    template<typename T>
    T *
    meta_decompress_3d(const unsigned char *compressed, size_t r1, size_t r2, size_t r3) {
        const unsigned char *compressed_pos = compressed;
        meta_params params;
        read_variable_from_src(compressed_pos, params);
        DSize_3d size(r1, r2, r3, params.block_size);
        double precision = 0;
        read_variable_from_src(compressed_pos, precision);
        int intv_radius = 0;
        read_variable_from_src(compressed_pos, intv_radius);
        meanInfo<T> mean_info;
        read_variable_from_src(compressed_pos, mean_info);
        size_t reg_count = 0;
        read_variable_from_src(compressed_pos, reg_count);
        // prepare unpred buffer for vectorization
        int est_unpred_count_per_index = size.num_blocks * size.block_size * 1;
        // if(!params.block_independant) est_unpred_count_per_index /= 20;
        int *unpred_count_buffer = read_array_from_src<int>(compressed_pos, size.block_size * size.block_size);
        T *unpred_data_buffer = (T *) malloc(size.block_size * size.block_size * est_unpred_count_per_index * sizeof(T));
        T *unpred_data_buffer_pos = unpred_data_buffer;
        for (int i = 0; i < size.block_size; i++) {
            for (int j = 0; j < size.block_size; j++) {
                memcpy(unpred_data_buffer_pos, compressed_pos, unpred_count_buffer[i * size.block_size + j] * sizeof(T));
                compressed_pos += unpred_count_buffer[i * size.block_size + j] * sizeof(T);
                unpred_data_buffer_pos += est_unpred_count_per_index;
            }
        }
        memset(unpred_count_buffer, 0, size.block_size * size.block_size * sizeof(int));
//	unsigned char * indicator = convertByteArray2IntArray_fast_1b_sz(size.num_blocks, compressed_pos, (size.num_blocks - 1)/8 + 1);
        int *indicator = Huffman_decode_tree_and_data(SELECTOR_RADIUS, size.num_blocks, compressed_pos);

        float *reg_params = nullptr;
        const float *reg_params_pos = nullptr;
        T precision_t = (T) precision;
        if (reg_count) {
            reg_params = decode_regression_coefficients(compressed_pos, reg_count, size.block_size, precision_t, params);
            reg_params_pos = (const float *) (reg_params + RegCoeffNum3d);
        }
        int *type = Huffman_decode_tree_and_data(4 * intv_radius, size.num_elements, compressed_pos);
        T *dec_data = (T *) malloc(size.num_elements * sizeof(T));
//    dec_data_sp_float = (float *) dec_data;

        const int *type_pos = type;
        const int *indicator_pos = indicator;
//        const float *reg_params_pos = reg_params;
        // add one more ghost layer
        size_t buffer_dim0_offset = (size.d2 + params.lorenzo_padding_layer) * (size.d3 + params.lorenzo_padding_layer);
        size_t buffer_dim1_offset = size.d3 + params.lorenzo_padding_layer;
        T *pred_buffer = (T *) malloc(
                (size.block_size + params.lorenzo_padding_layer) * (size.d2 + params.lorenzo_padding_layer) *
                (size.d3 + params.lorenzo_padding_layer) * sizeof(T));
        memset(pred_buffer, 0, (size.block_size + params.lorenzo_padding_layer) * (size.d2 + params.lorenzo_padding_layer) *
                               (size.d3 + params.lorenzo_padding_layer) * sizeof(T));
        auto *lorenzo_pred_and_decomp = decompress_lorenzo_3d_prediction<T>;
        if (params.prediction_dim == 2) lorenzo_pred_and_decomp = decompress_lorenzo_3d_as2d_prediction<T>;
        else if (params.prediction_dim == 1) lorenzo_pred_and_decomp = decompress_lorenzo_3d_as1d_prediction<T>;
        T *x_data_pos = dec_data;
        for (size_t i = 0; i < size.num_x; i++) {
            T *y_data_pos = x_data_pos;
            T *pred_buffer_pos = pred_buffer;
            for (size_t j = 0; j < size.num_y; j++) {
                T *z_data_pos = y_data_pos;
                for (size_t k = 0; k < size.num_z; k++) {
                    int size_x = ((i + 1) * size.block_size < size.d1) ? size.block_size : size.d1 - i * size.block_size;
                    int size_y = ((j + 1) * size.block_size < size.d2) ? size.block_size : size.d2 - j * size.block_size;
                    int size_z = ((k + 1) * size.block_size < size.d3) ? size.block_size : size.d3 - k * size.block_size;
                    if (*indicator_pos == SELECTOR_REGRESSION) {
                        // regression
                        decompress_regression_3d_prediction<T>(reg_params_pos, pred_buffer_pos, precision, intv_radius,
                                                               size_x, size_y, size_z, buffer_dim0_offset, buffer_dim1_offset,
                                                               size.dim0_offset, size.dim1_offset, type_pos, unpred_count_buffer,
                                                               unpred_data_buffer, est_unpred_count_per_index, z_data_pos,
                                                               params.lorenzo_padding_layer);
                        reg_params_pos += RegCoeffNum3d;
                    } else {
                        // Lorenzo
                        lorenzo_pred_and_decomp(mean_info, pred_buffer_pos, precision, intv_radius, size_x, size_y, size_z,
                                                buffer_dim0_offset, buffer_dim1_offset, size.dim0_offset, size.dim1_offset,
                                                type_pos,
                                                unpred_count_buffer, unpred_data_buffer, est_unpred_count_per_index, z_data_pos,
                                                params.lorenzo_padding_layer, *indicator_pos == SELECTOR_LORENZO_2LAYER);
                    }
                    pred_buffer_pos += size.block_size;
                    indicator_pos++;
                    z_data_pos += size_z;
                }
                y_data_pos += size.block_size * size.dim1_offset;
                pred_buffer_pos += size.block_size * buffer_dim1_offset - size.block_size * size.num_z;
            }
            memcpy(pred_buffer, pred_buffer + size.block_size * buffer_dim0_offset,
                   params.lorenzo_padding_layer * buffer_dim0_offset * sizeof(T));
            x_data_pos += size.block_size * size.dim0_offset;
        }
        free(pred_buffer);

        free(unpred_count_buffer);
        free(unpred_data_buffer);
        free(indicator);
        free(reg_params);
        free(type);
        return dec_data;
    }

    template<typename T>
    meta_compress_info meta_compress_decompress_3d(T *data, size_t num_elements, int r1, int r2, int r3, float precision,
                                                   meta_params params, bool use_decompress) {
        size_t result_size = 0;
        meta_compress_info compressInfo;

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);

        unsigned char *result = meta_compress_3d<T>(data, r1, r2, r3, precision, result_size, params, compressInfo);
        unsigned char *result_after_lossless = NULL;
//    size_t lossless_outsize = result_size;
        size_t lossless_outsize = meta_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size, &result_after_lossless);
        free(result);

        clock_gettime(CLOCK_REALTIME, &end);
        compressInfo.compress_time = (float) (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        cout << "LEVEL2 Meta compression time = " << compressInfo.compress_time << "s" << endl;

        float ratio = (num_elements * sizeof(T)) * 1.0 / lossless_outsize;
        compressInfo.ori_bytes = num_elements * sizeof(T);
        compressInfo.compress_bytes = lossless_outsize;
        compressInfo.ratio = compressInfo.ori_bytes * 1.0 / compressInfo.compress_bytes;
        cout << "Compressed size = " << lossless_outsize << endl;
        cout << "Compression ratio= " << ratio << endl;


        if (use_decompress) {
            clock_gettime(CLOCK_REALTIME, &start);
            size_t lossless_output = meta_lossless_decompress(ZSTD_COMPRESSOR, result_after_lossless, lossless_outsize, &result,
                                                              result_size);
            T *dec_data = meta_decompress_3d<T>(result, r1, r2, r3);
            clock_gettime(CLOCK_REALTIME, &end);
            compressInfo.decompress_time =
                    (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
            cout << "LEVEL2 Meta decompression time = " << compressInfo.decompress_time << "s" << endl;
            free(result_after_lossless);
            // writefile("dec_data.dat", dec_data, num_elements);

            SZ::verify<T>(data, dec_data, num_elements, compressInfo.psnr, compressInfo.nrmse);

            free(result);
            free(dec_data);
        } else {
            free(result_after_lossless);
        }

        return compressInfo;
    }


    template<typename T>
    meta_compress_info
    meta_compress_3d_sampling_wrapper(const T *data, size_t num_elements, int r1, int r2, int r3, float precision,
                                      meta_params params) {
        size_t result_size = 0;
        meta_compress_info compressInfo;

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);

        unsigned char *result = meta_compress_3d_sampling<T>(data, r1, r2, r3, precision, result_size, params, compressInfo);
        if (params.lossless) {
            unsigned char *result_after_lossless = NULL;
            size_t lossless_outsize = meta_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size, &result_after_lossless);
            compressInfo.compress_bytes = lossless_outsize;
            free(result_after_lossless);
        } else {
            compressInfo.compress_bytes = result_size;
        }
        clock_gettime(CLOCK_REALTIME, &end);
        compressInfo.compress_time = (float) (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
//    cout << "Compression time: " << std::setw(10) << std::setprecision(6) << compressInfo.compress_time << "s" << endl;

        compressInfo.ratio = compressInfo.ori_bytes * 1.0 / compressInfo.compress_bytes;
//    cout << "Compressed size = " << compressInfo.compress_bytes << endl;
//    cout << "!!!!!!!!!!!!!!!!!!!!! ratio  !!!!!!!!!!!!!= " << compressInfo.ratio << endl;

        free(result);

        return compressInfo;
    }

    template<typename T>
    unsigned char *
    meta_autotuning_3d(T *data, size_t r1, size_t r2, size_t r3, double relative_eb, size_t &compressed_size,
                       bool baseline = false, bool log = false, float sample_ratio = 0.05) {
        size_t num_elements = r1 * r2 * r3;
        float max = data[0];
        float min = data[0];
        for (int i = 1; i < num_elements; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }
        float precision = relative_eb * (max - min);

        if (baseline) {
            meta_params baseline_param(false, 6, 3, 0, true, false, true, precision);
            auto baseline_compress_info = meta_compress_decompress_3d<T>(data, num_elements, r1, r2, r3, precision,
                                                                         baseline_param, true);
            fprintf(stdout,
                    "Baseline: reb:%.1e, ratio:%.2f, compress_time:%.3f, PSNR:%.2f, NRMSE %.10e Ori_bytes %ld, Compressed_bytes %ld\n",
                    relative_eb, baseline_compress_info.ratio, baseline_compress_info.compress_time, baseline_compress_info.psnr,
                    baseline_compress_info.nrmse, baseline_compress_info.ori_bytes,
                    baseline_compress_info.compress_bytes);
        }

        char best_param_str[1000];
        char buffer[1000];

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);
        int sample_num = 0;

        int capacity = 65536 * 2;
        optimize_quant_invl_3d(data, r1, r2, r3, precision, capacity);
//    printf("tuning capacity: %d\n", capacity);

        float best_ratio = 0;
        meta_params best_params_stage1;
        {
            for (auto use_lorenzo:{true}) {
                for (auto use_lorenzo_2layer:{true, false}) {
//            for (auto use_lorenzo_2layer:{true}) {
                    if (use_lorenzo || use_lorenzo_2layer) {
                        list<int> pred_dim_set = {3, 2};
                        for (auto pred_dim: pred_dim_set) {
                            meta_params params(false, 6, pred_dim, 0, use_lorenzo,
                                               use_lorenzo_2layer, false, precision);
                            params.sample_ratio = sample_ratio * 1.2;
                            params.capacity = capacity;
                            params.lossless = false;
                            auto compress_info = meta_compress_3d_sampling_wrapper<T>(data, num_elements, r1, r2, r3, precision,
                                                                                      params);
                            sample_num++;

                            sprintf(buffer,
                                    "lorenzo:%d, lorenzo2:%d, pred_dim:%d\n",
                                    use_lorenzo, use_lorenzo_2layer, pred_dim);
                            if (log) {
                                fprintf(stdout, "stage:1, ratio:%.2f, reb:%.1e, compress_time:%.3f, %s",
                                        compress_info.ratio, relative_eb, compress_info.compress_time, buffer);
                            }
                            if (compress_info.ratio > best_ratio * 1.02) {
                                best_ratio = compress_info.ratio;
                                best_params_stage1 = params;
                                memcpy(best_param_str, buffer, 1000 * sizeof(char));
                            }
                        }
                    }
                }
            }
            if (log) {
                fprintf(stdout, "Stage:1, Best Ratio %.2f, Params %s", best_ratio, best_param_str);
            }
        }


        meta_params best_params_stage2 = best_params_stage1;
        meta_compress_info best_compress_info;
        best_ratio = 0;
        list<int> block_size_set = {5, 6, 7, 8};
        for (auto block_size:block_size_set) {
            list<double> reg_eb_base_set = {1};
            list<double> reg_eb_1_set = {block_size * 1.0};
            for (auto reg_eb_base:reg_eb_base_set) {
                for (auto reg_eb_1:reg_eb_1_set) {
                    meta_params params(false, block_size, best_params_stage1.prediction_dim, 0,
                                       best_params_stage1.use_lorenzo,
                                       best_params_stage1.use_lorenzo_2layer,
                                       true, precision, reg_eb_base,
                                       reg_eb_1);
                    params.sample_ratio = sample_ratio * 1.2;
                    params.capacity = capacity;
                    auto compress_info = meta_compress_3d_sampling_wrapper<T>(data, num_elements, r1, r2, r3,
                                                                              precision, params);
                    sample_num++;

                    sprintf(buffer,
                            "lorenzo:%d, lorenzo2:%d, regression:%d,"
                            "block_size:%d, pred_dim:%d, reg_base:%.1f, reg_1: %.1f,"
                            "lorenzo: %.0f, lorenzo2: %.0f, regression:%.0f\n",
                            best_params_stage1.use_lorenzo, best_params_stage1.use_lorenzo_2layer, true,
                            block_size, best_params_stage1.prediction_dim, reg_eb_base, reg_eb_1,
                            compress_info.lorenzo_count * 100.0 / compress_info.block_count,
                            compress_info.lorenzo2_count * 100.0 / compress_info.block_count,
                            compress_info.regression_count * 100.0 / compress_info.block_count);
                    if (log) {
                        fprintf(stdout, "stage:2, reb:%.1e, ratio:%.2f, compress_time:%.3f, %s",
                                relative_eb,
                                compress_info.ratio, compress_info.compress_time, buffer);
                    }
                    if (compress_info.ratio > best_ratio * 1.01) {
                        best_ratio = compress_info.ratio;
                        best_params_stage2 = params;
                        memcpy(best_param_str, buffer, 1000 * sizeof(char));
                        best_compress_info = compress_info;
                    }
                }
            }
        }
        if ((best_compress_info.lorenzo_count * 1.0 + best_compress_info.lorenzo2_count) / best_compress_info.block_count > 0.9) {
            meta_params params = best_params_stage2;
            params.use_lorenzo = best_params_stage1.use_lorenzo;
            params.use_lorenzo_2layer = best_params_stage1.use_lorenzo_2layer;
            params.prediction_dim = best_params_stage1.prediction_dim;
            auto compress_info = meta_compress_3d_sampling_wrapper<T>(data, num_elements, r1, r2, r3, precision, params);
            sample_num++;

            sprintf(buffer,
                    "lorenzo:%d, lorenzo2:%d, regression:%d,"
                    "block_size:%d, pred_dim:%d, reg_base:%.1f, reg_1: %.1f"
                    "lorenzo: %.0f, lorenzo2: %.0f, regression:%.0f\n",
                    params.use_lorenzo, params.use_lorenzo_2layer, false,
                    params.block_size, params.prediction_dim, 0.0, 0.0,
                    compress_info.lorenzo_count * 100.0 / compress_info.block_count,
                    compress_info.lorenzo2_count * 100.0 / compress_info.block_count,
                    compress_info.regression_count * 100.0 / compress_info.block_count);
            if (log) {
                fprintf(stdout, "stage:2, reb:%.1e, ratio:%.2f, compress_time:%.3f, %s", relative_eb,
                        compress_info.ratio, compress_info.compress_time, buffer);
            }
            if (compress_info.ratio > best_ratio * 1.01) {
                best_ratio = compress_info.ratio;
                best_params_stage2 = params;
                memcpy(best_param_str, buffer, 1000 * sizeof(char));
            }
        }
        if (log) {
            fprintf(stdout, "Stage:2, Best Ratio %.2f, Params %s", best_ratio, best_param_str);
        }

        meta_params best_params_stage3;
        if (relative_eb < 1.01e-6 && best_ratio > 5) {
            best_ratio = 0;
            list<int> capacity_set = {capacity, 16384};
            for (auto capacity1:capacity_set) {
                best_params_stage2.sample_ratio = sample_ratio;
                best_params_stage2.capacity = capacity1;
                auto compress_info = meta_compress_3d_sampling_wrapper<T>(data, num_elements, r1, r2, r3, precision,
                                                                          best_params_stage2);
                sample_num++;
                if (log) {
                    fprintf(stdout,
                            "stage:3, reb:%.1e, ratio:%.2f, compress_time:%.3f, capacity:%d, %s",
                            relative_eb, compress_info.ratio, compress_info.compress_time, capacity1, best_param_str);
                }
                if (compress_info.ratio > best_ratio * 1.01) {
                    best_ratio = compress_info.ratio;
                    best_params_stage3 = best_params_stage2;
                }
            }
        } else {
            best_params_stage3 = best_params_stage2;
            best_params_stage3.capacity = capacity;
        }


        clock_gettime(CLOCK_REALTIME, &end);
        float sample_time = (float) (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        auto compress_info = meta_compress_decompress_3d(data, num_elements, r1, r2, r3, precision, best_params_stage3,
                                                         true);
        fprintf(stdout,
                "FINAL: reb:%.1e, ratio %.2f, compress_time:%.3f, capacity:%d, PSNR:%.2f, NRMSE %.10e, sample_time:%.1f, sample_num:%d, %s\n",
                relative_eb, compress_info.ratio, compress_info.compress_time, best_params_stage3.capacity,
                compress_info.psnr,
                compress_info.nrmse, sample_time, sample_num,
                best_param_str);
        return nullptr;
    }


}


#endif