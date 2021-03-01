#ifndef _meta_compress_block_processing_hpp
#define _meta_compress_block_processing_hpp

#include "meta_def.hpp"

namespace SZMETA {


    template<typename T>
    inline void
    compute_regression_coeffcients_3d(const T *data_pos, int size_x, int size_y, int size_z, size_t dim0_offset,
                                      size_t dim1_offset, float *reg_params_pos) {
        /*Calculate regression coefficients*/
        const T *cur_data_pos = data_pos;
        float fx = 0.0;
        float fy = 0.0;
        float fz = 0.0;
        float f = 0;
        float sum_x, sum_y;
        T curData;
        for (int i = 0; i < size_x; i++) {
            sum_x = 0;
            for (int j = 0; j < size_y; j++) {
                sum_y = 0;
                for (int k = 0; k < size_z; k++) {
                    curData = *cur_data_pos;
                    sum_y += curData;
                    fz += curData * k;
                    cur_data_pos++;
                }
                fy += sum_y * j;
                sum_x += sum_y;
                cur_data_pos += (dim1_offset - size_z);
            }
            fx += sum_x * i;
            f += sum_x;
            cur_data_pos += (dim0_offset - size_y * dim1_offset);
        }
        float coeff = 1.0 / (size_x * size_y * size_z);
        reg_params_pos[0] = (2 * fx / (size_x - 1) - f) * 6 * coeff / (size_x + 1);
        reg_params_pos[1] = (2 * fy / (size_y - 1) - f) * 6 * coeff / (size_y + 1);
        reg_params_pos[2] = (2 * fz / (size_z - 1) - f) * 6 * coeff / (size_z + 1);
        reg_params_pos[3] = f * coeff - ((size_x - 1) * reg_params_pos[0] / 2 + (size_y - 1) * reg_params_pos[1] / 2 +
                                         (size_z - 1) * reg_params_pos[2] / 2);
    }


    template<typename T>
    inline T
    regression_predict_3d(const float *reg_params_pos, int x, int y, int z) {
        return reg_params_pos[0] * x + reg_params_pos[1] * y + reg_params_pos[2] * z + reg_params_pos[3];
    }


    template<typename T>
    inline T
    lorenzo_predict_1d(const T *data_pos, size_t dim0_offset) {
        return data_pos[-1];
    }

    template<typename T>
    inline T
    lorenzo_predict_1d_2layer(const T *data_pos, size_t dim0_offset) {
        return 2 * data_pos[-1] - data_pos[-2];
    }


    template<typename T>
    inline T
    lorenzo_predict_2d(const T *data_pos, size_t dim0_offset, size_t dim1_offset) {
        return data_pos[-1] + data_pos[-dim0_offset] - data_pos[-1 - dim0_offset];
    }

    template<typename T>
    inline T
    lorenzo_predict_2d_2layer(const T *data_pos, size_t dim0_offset, size_t dim1_offset) {
        return 2 * data_pos[-dim0_offset]
               - data_pos[-2 * dim0_offset]
               + 2 * data_pos[-1]
               - 4 * data_pos[-1 - dim0_offset]
               + 2 * data_pos[-1 - 2 * dim0_offset]
               - data_pos[-2]
               + 2 * data_pos[-2 - dim0_offset]
               - data_pos[-2 - 2 * dim0_offset];
    }

    template<typename T>
    inline T
    lorenzo_predict_3d(const T *data_pos, size_t dim0_offset, size_t dim1_offset) {
        return data_pos[-1] + data_pos[-dim1_offset] + data_pos[-dim0_offset]
               - data_pos[-dim1_offset - 1] - data_pos[-dim0_offset - 1]
               - data_pos[-dim0_offset - dim1_offset] + data_pos[-dim0_offset - dim1_offset - 1];
    }

    template<typename T>
    inline T
    lorenzo_predict_3d_2layer(const T *data_pos, size_t dim0_offset, size_t dim1_offset) {
        return 2 * data_pos[-1]
               - data_pos[-2]
               + 2 * data_pos[-dim1_offset]
               - 4 * data_pos[-dim1_offset - 1]
               + 2 * data_pos[-dim1_offset - 2]
               - data_pos[-2 * dim1_offset]
               + 2 * data_pos[-2 * dim1_offset - 1]
               - data_pos[-2 * dim1_offset - 2]
               + 2 * data_pos[-dim0_offset]
               - 4 * data_pos[-dim0_offset - 1]
               + 2 * data_pos[-dim0_offset - 2]
               - 4 * data_pos[-dim0_offset - dim1_offset]
               + 8 * data_pos[-dim0_offset - dim1_offset - 1]
               - 4 * data_pos[-dim0_offset - dim1_offset - 2]
               + 2 * data_pos[-dim0_offset - 2 * dim1_offset]
               - 4 * data_pos[-dim0_offset - 2 * dim1_offset - 1]
               + 2 * data_pos[-dim0_offset - 2 * dim1_offset - 2]
               - data_pos[-2 * dim0_offset]
               + 2 * data_pos[-2 * dim0_offset - 1]
               - data_pos[-2 * dim0_offset - 2]
               + 2 * data_pos[-2 * dim0_offset - dim1_offset]
               - 4 * data_pos[-2 * dim0_offset - dim1_offset - 1]
               + 2 * data_pos[-2 * dim0_offset - dim1_offset - 2]
               - data_pos[-2 * dim0_offset - 2 * dim1_offset]
               + 2 * data_pos[-2 * dim0_offset - 2 * dim1_offset - 1]
               - data_pos[-2 * dim0_offset - 2 * dim1_offset - 2];
    }

    template<typename T>
    inline void
    compress_regression_3d_predict(const T *data_pos, const float *reg_params_pos, T *buffer, T precision,
                                   T recip_precision, int capacity,
                                   int intv_radius, int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                   size_t buffer_dim1_offset,
                                   size_t dim0_offset, size_t dim1_offset, int *&type_pos,
                                   int *unpred_count_buffer, T *unpred_buffer, size_t offset, int lorenzo_layer) {
        for (int i = 0; i < size_x; i++) {
//        const T *cur_data_pos = data_pos + i * dim0_offset;
            // T * buffer_pos = buffer + (i+1)*buffer_dim0_offset + buffer_dim1_offset + 1;
            T *buffer_pos =
                    buffer + (i + lorenzo_layer) * buffer_dim0_offset + lorenzo_layer * buffer_dim1_offset + lorenzo_layer;
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    T cur_data = data_pos[i * dim0_offset + j * dim1_offset + k];
                    T pred = (T) (reg_params_pos[0] * (float) i + reg_params_pos[1] * (float) j + reg_params_pos[2] * (float) k +
                                  reg_params_pos[3]);

                    T diff = cur_data - pred;
                    int quant_index = (int) (fabs(diff) * recip_precision) + 1;
                    if (quant_index < capacity) {
                        quant_index >>= 1;
                        int half_index = quant_index;
                        quant_index <<= 1;
                        if (diff < 0) {
                            quant_index = -quant_index;
                            type_pos[j * size_z + k] = intv_radius - half_index;
                        } else type_pos[j * size_z + k] = intv_radius + half_index;
                        T decompressed_data = pred + (T) quant_index * precision;
                        if (fabs(cur_data - decompressed_data) > precision) {
                            int index = j * size_z + k;
                            type_pos[index] = 0;
                            unpred_buffer[index * offset + unpred_count_buffer[index]] = buffer_pos[j * buffer_dim1_offset +
                                                                                                    k] = cur_data;
                            unpred_count_buffer[index]++;
                        } else buffer_pos[j * buffer_dim1_offset + k] = decompressed_data;
                    } else {
                        int index = j * size_z + k;
                        type_pos[index] = 0;
                        unpred_buffer[index * offset + unpred_count_buffer[index]] = buffer_pos[j * buffer_dim1_offset +
                                                                                                k] = cur_data;
                        unpred_count_buffer[index]++;
                    }
//			 	printf("element %.3f, predicted %.3f, quant %d\n", cur_data, pred, type_pos[j*size_z+k]);
                }
            }
            type_pos += size_y * size_z;
        }
//    exit(0);

    }

    template<typename T>
    inline void
    compress_lorenzo_3d_as1d_predict(const meanInfo<T> &mean_info, const T *data_pos, T *buffer, T precision,
                                     T recip_precision, int capacity, int intv_radius,
                                     int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                     size_t buffer_dim1_offset,
                                     size_t dim0_offset, size_t dim1_offset, int *&type_pos, int *unpred_count_buffer,
                                     T *unpred_buffer, size_t offset, int layer, bool use_2layer) {
        const T *cur_data_pos = data_pos;
//	T * buffer_pos = buffer + buffer_dim0_offset + buffer_dim1_offset + 1;
        T *buffer_pos = buffer + layer * (buffer_dim0_offset + buffer_dim1_offset + 1);

        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    if (mean_info.use_mean && (fabs(cur_data_pos[k] - mean_info.mean) < precision)) {
                        type_pos[k] = 1;
                        buffer_pos[k] = mean_info.mean;
                    } else {
                        T *cur_buffer_pos = buffer_pos + k;
                        T cur_data = cur_data_pos[k];
                        T pred;
                        if (use_2layer) {
                            pred = 2 * cur_buffer_pos[-1] - cur_buffer_pos[-2];
                        } else {
                            pred = cur_buffer_pos[-1];
                        }
                        T diff = cur_data - pred;
                        int quant_index = (int) (fabs(diff) * recip_precision) + 1;
                        if (quant_index < capacity) {
                            quant_index >>= 1;
                            int half_index = quant_index;
                            quant_index <<= 1;
                            if (diff < 0) {
                                quant_index = -quant_index;
                                half_index = -half_index;
                            }
                            type_pos[k] = half_index + intv_radius;
                            T decompressed_data = pred + (T) quant_index * precision;
                            if (fabs(cur_data - decompressed_data) > precision) {
                                int index = j * size_z + k;
                                unpred_buffer[index * offset + unpred_count_buffer[index]] = *cur_buffer_pos = cur_data;
                                unpred_count_buffer[index]++;
                                type_pos[k] = 0;
                            } else *cur_buffer_pos = decompressed_data;
                        } else {
                            int index = j * size_z + k;
                            unpred_buffer[index * offset + unpred_count_buffer[index]] = *cur_buffer_pos = cur_data;
                            unpred_count_buffer[index]++;
                            type_pos[k] = 0;
                        }
                    }
                }
                type_pos += size_z;
                buffer_pos += buffer_dim1_offset;
                cur_data_pos += dim1_offset;
            }
            buffer_pos += buffer_dim0_offset - size_y * buffer_dim1_offset;
            cur_data_pos += dim0_offset - size_y * dim1_offset;
        }
    }

    template<typename T>
    inline void
    compress_lorenzo_3d_as2d_predict(const meanInfo<T> &mean_info, const T *data_pos, T *buffer, T precision,
                                     T recip_precision, int capacity, int intv_radius,
                                     int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                     size_t buffer_dim1_offset,
                                     size_t dim0_offset, size_t dim1_offset, int *&type_pos, int *unpred_count_buffer,
                                     T *unpred_buffer, size_t offset, int layer, bool use_2layer) {
        const T *cur_data_pos = data_pos;
//	T * buffer_pos = buffer + buffer_dim0_offset + buffer_dim1_offset + 1;
        T *buffer_pos = buffer + layer * (buffer_dim0_offset + buffer_dim1_offset + 1);

        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    if (mean_info.use_mean && (fabs(cur_data_pos[k] - mean_info.mean) < precision)) {
                        type_pos[k] = 1;
                        buffer_pos[k] = mean_info.mean;
                    } else {
                        T *cur_buffer_pos = buffer_pos + k;
                        T cur_data = cur_data_pos[k];
                        T pred;
                        if (use_2layer) {
                            pred = 2 * cur_buffer_pos[-buffer_dim0_offset]
                                   - cur_buffer_pos[-2 * buffer_dim0_offset]
                                   + 2 * cur_buffer_pos[-1]
                                   - 4 * cur_buffer_pos[-1 - buffer_dim0_offset]
                                   + 2 * cur_buffer_pos[-1 - 2 * buffer_dim0_offset]
                                   - cur_buffer_pos[-2]
                                   + 2 * cur_buffer_pos[-2 - buffer_dim0_offset]
                                   - cur_buffer_pos[-2 - 2 * buffer_dim0_offset];
                        } else {
                            pred = cur_buffer_pos[-1] + cur_buffer_pos[-buffer_dim0_offset] -
                                   cur_buffer_pos[-1 - buffer_dim0_offset];

                        }
                        T diff = cur_data - pred;
                        int quant_index = (int) (fabs(diff) * recip_precision) + 1;
                        if (quant_index < capacity) {
                            quant_index >>= 1;
                            int half_index = quant_index;
                            quant_index <<= 1;
                            if (diff < 0) {
                                quant_index = -quant_index;
                                half_index = -half_index;
                            }
                            type_pos[k] = half_index + intv_radius;
                            T decompressed_data = pred + (T) quant_index * precision;
                            if (fabs(cur_data - decompressed_data) > precision) {
                                int index = j * size_z + k;
                                unpred_buffer[index * offset + unpred_count_buffer[index]] = *cur_buffer_pos = cur_data;
                                unpred_count_buffer[index]++;
                                type_pos[k] = 0;
                            } else *cur_buffer_pos = decompressed_data;
//                        if(cur_data_pos + k - ori_data_sp_float == -1){
//                            printf("compress index out of bound, data=%.4f, diff=%.4f, dec_data=%.4f, type=%d, quant_radius=%d, pred=%.4f\n", cur_data, fabs(cur_data - decompressed_data), decompressed_data, type_pos[k], intv_radius, pred);
//                            exit(0);
//                        }
                        } else {
                            int index = j * size_z + k;
                            unpred_buffer[index * offset + unpred_count_buffer[index]] = *cur_buffer_pos = cur_data;
                            unpred_count_buffer[index]++;
                            type_pos[k] = 0;
                        }
                    }
                }
                type_pos += size_z;
                buffer_pos += buffer_dim1_offset;
                cur_data_pos += dim1_offset;
            }
            buffer_pos += buffer_dim0_offset - size_y * buffer_dim1_offset;
            cur_data_pos += dim0_offset - size_y * dim1_offset;
        }
    }

    template<typename T>
    inline void
    compress_lorenzo_3d_predict(const meanInfo<T> &mean_info, const T *data_pos, T *buffer, T precision,
                                T recip_precision, int capacity, int intv_radius,
                                int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                size_t buffer_dim1_offset,
                                size_t dim0_offset, size_t dim1_offset, int *&type_pos, int *unpred_count_buffer,
                                T *unpred_buffer, size_t offset, int padding_layer,
                                bool use_2layer) {
        const T *cur_data_pos = data_pos;
        T *buffer_pos = buffer + padding_layer * (buffer_dim0_offset + buffer_dim1_offset + 1);
        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    if (mean_info.use_mean && (fabs(cur_data_pos[k] - mean_info.mean) < precision)) {
                        type_pos[k] = 1;
                        buffer_pos[k] = mean_info.mean;
                    } else {
                        T *cur_buffer_pos = buffer_pos + k;
                        T cur_data = cur_data_pos[k];
                        T pred;
                        if (use_2layer) {
                            pred = 2 * cur_buffer_pos[-1]
                                   - cur_buffer_pos[-2]
                                   + 2 * cur_buffer_pos[-buffer_dim1_offset]
                                   - 4 * cur_buffer_pos[-buffer_dim1_offset - 1]
                                   + 2 * cur_buffer_pos[-buffer_dim1_offset - 2]
                                   - cur_buffer_pos[-2 * buffer_dim1_offset]
                                   + 2 * cur_buffer_pos[-2 * buffer_dim1_offset - 1]
                                   - cur_buffer_pos[-2 * buffer_dim1_offset - 2]
                                   + 2 * cur_buffer_pos[-buffer_dim0_offset]
                                   - 4 * cur_buffer_pos[-buffer_dim0_offset - 1]
                                   + 2 * cur_buffer_pos[-buffer_dim0_offset - 2]
                                   - 4 * cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset]
                                   + 8 * cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset - 1]
                                   - 4 * cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset - 2]
                                   + 2 * cur_buffer_pos[-buffer_dim0_offset - 2 * buffer_dim1_offset]
                                   - 4 * cur_buffer_pos[-buffer_dim0_offset - 2 * buffer_dim1_offset - 1]
                                   + 2 * cur_buffer_pos[-buffer_dim0_offset - 2 * buffer_dim1_offset - 2]
                                   - cur_buffer_pos[-2 * buffer_dim0_offset]
                                   + 2 * cur_buffer_pos[-2 * buffer_dim0_offset - 1]
                                   - cur_buffer_pos[-2 * buffer_dim0_offset - 2]
                                   + 2 * cur_buffer_pos[-2 * buffer_dim0_offset - buffer_dim1_offset]
                                   - 4 * cur_buffer_pos[-2 * buffer_dim0_offset - buffer_dim1_offset - 1]
                                   + 2 * cur_buffer_pos[-2 * buffer_dim0_offset - buffer_dim1_offset - 2]
                                   - cur_buffer_pos[-2 * buffer_dim0_offset - 2 * buffer_dim1_offset]
                                   + 2 * cur_buffer_pos[-2 * buffer_dim0_offset - 2 * buffer_dim1_offset - 1]
                                   - cur_buffer_pos[-2 * buffer_dim0_offset - 2 * buffer_dim1_offset - 2];
                        } else {
                            pred = cur_buffer_pos[-1] + cur_buffer_pos[-buffer_dim1_offset] + cur_buffer_pos[-buffer_dim0_offset]
                                   - cur_buffer_pos[-buffer_dim1_offset - 1] - cur_buffer_pos[-buffer_dim0_offset - 1]
                                   - cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset] +
                                   cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset - 1];
                        }
                        T diff = cur_data - pred;
                        int quant_index = (int) (fabs(diff) * recip_precision) + 1;
                        if (quant_index < capacity) {
                            quant_index >>= 1;
                            int half_index = quant_index;
                            quant_index <<= 1;
                            if (diff < 0) {
                                quant_index = -quant_index;
                                half_index = -half_index;
                            }
                            type_pos[k] = half_index + intv_radius;
                            T decompressed_data = pred + (T) quant_index * precision;
                            if (fabs(cur_data - decompressed_data) > precision) {
                                int index = j * size_z + k;
                                unpred_buffer[index * offset + unpred_count_buffer[index]] = *cur_buffer_pos = cur_data;
                                unpred_count_buffer[index]++;
                                type_pos[k] = 0;
                            } else *cur_buffer_pos = decompressed_data;
//                        if(cur_data_pos + k - ori_data_sp_float == -1){
//                            printf("compress index out of bound, data=%.4f, diff=%.4f, dec_data=%.4f, type=%d, quant_radius=%d, pred=%.4f\n", cur_data, fabs(cur_data - decompressed_data), decompressed_data, type_pos[k], intv_radius, pred);
//                            exit(0);
//                        }
                        } else {
                            int index = j * size_z + k;
                            unpred_buffer[index * offset + unpred_count_buffer[index]] = *cur_buffer_pos = cur_data;
                            unpred_count_buffer[index]++;
                            type_pos[k] = 0;
                        }
                    }
                }
                type_pos += size_z;
                buffer_pos += buffer_dim1_offset;
                cur_data_pos += dim1_offset;
            }
            buffer_pos += buffer_dim0_offset - size_y * buffer_dim1_offset;
            cur_data_pos += dim0_offset - size_y * dim1_offset;
        }
    }


    template<typename T>
    inline void
    decompress_regression_3d_prediction(const float *reg_params_pos, T *buffer, T precision, int intv_radius,
                                        int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset,
                                        size_t dim0_offset, size_t dim1_offset, const int *&type_pos, int *unpred_count_buffer,
                                        const T *unpred_data_buffer, const int offset, T *dec_data_pos, int lorenzo_layer) {
        T *cur_data_pos = dec_data_pos;
        T *buffer_pos = buffer + lorenzo_layer * (buffer_dim0_offset + buffer_dim1_offset + 1);
        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    int index = j * size_z + k;
                    int type_val = type_pos[index];
                    if (type_val == 0) {
                        cur_data_pos[j * dim1_offset + k] = buffer_pos[j * buffer_dim1_offset + k] = unpred_data_buffer[
                                index * offset + unpred_count_buffer[index]];
                        unpred_count_buffer[index]++;
                    } else {

                        T pred = (T) (reg_params_pos[0] * (float) i + reg_params_pos[1] * (float) j +
                                      reg_params_pos[2] * (float) k +
                                      reg_params_pos[3]);
                        cur_data_pos[j * dim1_offset + k] = buffer_pos[j * buffer_dim1_offset + k] =
                                pred + (T) (2 * (type_val - intv_radius)) * precision;
                    }
                }
            }
            type_pos += size_y * size_z;
            cur_data_pos += dim0_offset;
            buffer_pos += buffer_dim0_offset;
        }
    }

// block-independant use_lorenzo pred & decompress
    template<typename T>
    inline void
    decompress_lorenzo_3d_as1d_prediction(const meanInfo<T> &mean_info, T *buffer, T precision, int intv_radius,
                                          int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                          size_t buffer_dim1_offset,
                                          size_t dim0_offset, size_t dim1_offset,
                                          const int *&type_pos, int *unpred_count_buffer, const T *unpred_data_buffer,
                                          const int offset, T *dec_data_pos, const int layer,
                                          bool use_2layer) {
        T *cur_data_pos = dec_data_pos;
//	T * buffer_pos = buffer + buffer_dim0_offset + buffer_dim1_offset + 1;
        T *buffer_pos = buffer + layer * (buffer_dim0_offset + buffer_dim1_offset + 1);

        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    int index = j * size_z + k;
                    int type_val = type_pos[index];
                    if ((mean_info.use_mean) && (type_val == 1)) cur_data_pos[k] = buffer_pos[k] = mean_info.mean;
                    else {
                        T *cur_buffer_pos = buffer_pos + k;
                        if (type_val == 0) {
                            cur_data_pos[k] = *cur_buffer_pos = unpred_data_buffer[index * offset + unpred_count_buffer[index]];
                            unpred_count_buffer[index]++;
                        } else {
                            T pred;
                            if (use_2layer) {
                                pred = 2 * cur_buffer_pos[-1] - cur_buffer_pos[-2];
                            } else {
                                pred = cur_buffer_pos[-1];
                            }
                            cur_data_pos[k] = *cur_buffer_pos = pred + (T) (2 * (type_val - intv_radius)) * precision;
                        }
                    }
                }
                buffer_pos += buffer_dim1_offset;
                cur_data_pos += dim1_offset;
            }
            type_pos += size_y * size_z;
            buffer_pos += buffer_dim0_offset - size_y * buffer_dim1_offset;
            cur_data_pos += dim0_offset - size_y * dim1_offset;
        }
    }

    template<typename T>
    inline void
    decompress_lorenzo_3d_as2d_prediction(const meanInfo<T> &mean_info, T *buffer, T precision, int intv_radius,
                                          int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                          size_t buffer_dim1_offset,
                                          size_t dim0_offset, size_t dim1_offset,
                                          const int *&type_pos, int *unpred_count_buffer, const T *unpred_data_buffer,
                                          const int offset, T *dec_data_pos, const int layer,
                                          bool use_2layer) {
        T *cur_data_pos = dec_data_pos;
//	T * buffer_pos = buffer + buffer_dim0_offset + buffer_dim1_offset + 1;
        T *buffer_pos = buffer + layer * (buffer_dim0_offset + buffer_dim1_offset + 1);

        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    int index = j * size_z + k;
                    int type_val = type_pos[index];
                    if ((mean_info.use_mean) && (type_val == 1)) cur_data_pos[k] = buffer_pos[k] = mean_info.mean;
                    else {
                        T *cur_buffer_pos = buffer_pos + k;
                        if (type_val == 0) {
                            cur_data_pos[k] = *cur_buffer_pos = unpred_data_buffer[index * offset + unpred_count_buffer[index]];
                            unpred_count_buffer[index]++;
                        } else {
                            T pred;
                            if (use_2layer) {
                                pred = 2 * cur_buffer_pos[-buffer_dim0_offset]
                                       - cur_buffer_pos[-2 * buffer_dim0_offset]
                                       + 2 * cur_buffer_pos[-1]
                                       - 4 * cur_buffer_pos[-1 - buffer_dim0_offset]
                                       + 2 * cur_buffer_pos[-1 - 2 * buffer_dim0_offset]
                                       - cur_buffer_pos[-2]
                                       + 2 * cur_buffer_pos[-2 - buffer_dim0_offset]
                                       - cur_buffer_pos[-2 - 2 * buffer_dim0_offset];
                            } else {
                                pred = cur_buffer_pos[-1] + cur_buffer_pos[-buffer_dim0_offset] -
                                       cur_buffer_pos[-1 - buffer_dim0_offset];
                            }
                            cur_data_pos[k] = *cur_buffer_pos = pred + (T) (2 * (type_val - intv_radius)) * precision;
//                        if(cur_data_pos + k - dec_data_sp_float == -1){
//                            printf("index out of bound, dec_data=%.4f, type=%d, pred=%.4f\n", *cur_buffer_pos, type_val, pred);
//                            exit(0);
//                        }
                        }
                    }
                }
                buffer_pos += buffer_dim1_offset;
                cur_data_pos += dim1_offset;
            }
            type_pos += size_y * size_z;
            buffer_pos += buffer_dim0_offset - size_y * buffer_dim1_offset;
            cur_data_pos += dim0_offset - size_y * dim1_offset;
        }
    }

    template<typename T>
    inline void
    decompress_lorenzo_3d_prediction(const meanInfo<T> &mean_info, T *buffer, T precision, int intv_radius,
                                     int size_x, int size_y, int size_z, size_t buffer_dim0_offset, size_t buffer_dim1_offset,
                                     size_t dim0_offset, size_t dim1_offset,
                                     const int *&type_pos, int *unpred_count_buffer, const T *unpred_data_buffer,
                                     const int offset,
                                     T *dec_data_pos, const int layer,
                                     bool use_2layer) {
        T *cur_data_pos = dec_data_pos;
        T *buffer_pos = buffer + layer * (buffer_dim0_offset + buffer_dim1_offset + 1);
        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    int index = j * size_z + k;
                    int type_val = type_pos[index];
                    if ((mean_info.use_mean) && (type_val == 1)) cur_data_pos[k] = buffer_pos[k] = mean_info.mean;
                    else {
                        T *cur_buffer_pos = buffer_pos + k;
                        if (type_val == 0) {
                            cur_data_pos[k] = *cur_buffer_pos = unpred_data_buffer[index * offset + unpred_count_buffer[index]];
                            unpred_count_buffer[index]++;
                        } else {
                            T pred;
                            if (use_2layer) {
                                pred = 2 * cur_buffer_pos[-1]
                                       - cur_buffer_pos[-2]
                                       + 2 * cur_buffer_pos[-buffer_dim1_offset]
                                       - 4 * cur_buffer_pos[-buffer_dim1_offset - 1]
                                       + 2 * cur_buffer_pos[-buffer_dim1_offset - 2]
                                       - cur_buffer_pos[-2 * buffer_dim1_offset]
                                       + 2 * cur_buffer_pos[-2 * buffer_dim1_offset - 1]
                                       - cur_buffer_pos[-2 * buffer_dim1_offset - 2]
                                       + 2 * cur_buffer_pos[-buffer_dim0_offset]
                                       - 4 * cur_buffer_pos[-buffer_dim0_offset - 1]
                                       + 2 * cur_buffer_pos[-buffer_dim0_offset - 2]
                                       - 4 * cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset]
                                       + 8 * cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset - 1]
                                       - 4 * cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset - 2]
                                       + 2 * cur_buffer_pos[-buffer_dim0_offset - 2 * buffer_dim1_offset]
                                       - 4 * cur_buffer_pos[-buffer_dim0_offset - 2 * buffer_dim1_offset - 1]
                                       + 2 * cur_buffer_pos[-buffer_dim0_offset - 2 * buffer_dim1_offset - 2]
                                       - cur_buffer_pos[-2 * buffer_dim0_offset]
                                       + 2 * cur_buffer_pos[-2 * buffer_dim0_offset - 1]
                                       - cur_buffer_pos[-2 * buffer_dim0_offset - 2]
                                       + 2 * cur_buffer_pos[-2 * buffer_dim0_offset - buffer_dim1_offset]
                                       - 4 * cur_buffer_pos[-2 * buffer_dim0_offset - buffer_dim1_offset - 1]
                                       + 2 * cur_buffer_pos[-2 * buffer_dim0_offset - buffer_dim1_offset - 2]
                                       - cur_buffer_pos[-2 * buffer_dim0_offset - 2 * buffer_dim1_offset]
                                       + 2 * cur_buffer_pos[-2 * buffer_dim0_offset - 2 * buffer_dim1_offset - 1]
                                       - cur_buffer_pos[-2 * buffer_dim0_offset - 2 * buffer_dim1_offset - 2];
                            } else {
                                pred = cur_buffer_pos[-1] + cur_buffer_pos[-buffer_dim1_offset] +
                                       cur_buffer_pos[-buffer_dim0_offset]
                                       - cur_buffer_pos[-buffer_dim1_offset - 1] - cur_buffer_pos[-buffer_dim0_offset - 1]
                                       - cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset] +
                                       cur_buffer_pos[-buffer_dim0_offset - buffer_dim1_offset - 1];
                            }
                            cur_data_pos[k] = *cur_buffer_pos = pred + (T) (2 * (type_val - intv_radius)) * precision;
//                        if(cur_data_pos + k - dec_data_sp_float == -1){
//                            printf("index out of bound, dec_data=%.4f, type=%d, pred=%.4f\n", *cur_buffer_pos, type_val, pred);
//                            exit(0);
//                        }
                        }
                    }
                }
                buffer_pos += buffer_dim1_offset;
                cur_data_pos += dim1_offset;
            }
            type_pos += size_y * size_z;
            buffer_pos += buffer_dim0_offset - size_y * buffer_dim1_offset;
            cur_data_pos += dim0_offset - size_y * dim1_offset;
        }
    }


}
#endif