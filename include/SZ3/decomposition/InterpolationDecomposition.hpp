#ifndef _SZ_INTERPOLATION_DECOMPOSITION_HPP
#define _SZ_INTERPOLATION_DECOMPOSITION_HPP

#include <cmath>
#include <cstring>

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"

namespace SZ3 {
template <class T, uint N, class Quantizer>
class InterpolationDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    InterpolationDecomposition(const Config &conf, Quantizer quantizer) : quantizer(quantizer) {
        static_assert(std::is_base_of<concepts::QuantizerInterface<T, int>, Quantizer>::value,
                      "must implement the quantizer interface");
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        init();

        this->quant_inds = quant_inds.data();
        double eb = quantizer.get_eb();

        if (anchorStride == 0){
            *dec_data = quantizer.recover(0, this->quant_inds[quant_index++]);
        } 
        else{
            recover_anchor_grid(dec_data);                   
            interpolation_level--;           
        }

        for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
            if (alpha < 0){
                if (level >= 3) {
                    quantizer.set_eb(eb * eb_ratio);
                } else {
                    quantizer.set_eb(eb);
                }
            }
            else if (alpha >= 1){
                double cur_ratio = pow(alpha, level - 1);
                if (cur_ratio > beta){
                    cur_ratio = beta;
                }        
                quantizer.set_eb(eb / cur_ratio);
            }
            size_t stride = 1U << (level - 1);
            auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(
                dec_data, std::begin(global_dimensions), std::end(global_dimensions), stride * blocksize, 0);
            auto inter_begin = inter_block_range->begin();
            auto inter_end = inter_block_range->end();
            auto cur_interpolator_id = level >=3 ? 0 : interpolator_id;
            for (auto block = inter_begin; block != inter_end; ++block) {
                auto end_idx = block.get_global_index();
                for (int i = 0; i < N; i++) {
                    end_idx[i] += stride * blocksize;
                    if (end_idx[i] > global_dimensions[i] - 1) {
                        end_idx[i] = global_dimensions[i] - 1;
                    }
                }
                block_interpolation(dec_data, block.get_global_index(), end_idx, PB_recover,
                                    interpolators[interpolator_id], direction_sequence_id, stride);
            }
        }
        quantizer.postdecompress_data();
        return dec_data;
    }

    // compress given the error bound
    std::vector<int> compress(const Config &conf, T *data) override {
        std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
        blocksize = 32;
        interpolator_id = conf.interpAlgo;
        direction_sequence_id = conf.interpDirection;
        anchorStride = conf.interp_anchorStride;
        alpha = conf.interp_alpha;
        beta = conf.interp_beta;

        init();
        std::vector<int> quant_inds_vec(num_elements);
        quant_inds = quant_inds_vec.data();
        double eb = quantizer.get_eb();

        /*
        if(conf.tuning){
            auto range = std::make_shared<multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                                     std::end(global_dimensions), blocksize, 0);
            for (auto block = range->begin(); block != range->end(); ++block) {
                auto block_global_idx = block.get_global_index();
                auto interp_end_idx = block.get_global_index();
                uint max_interp_level = 1;
                for (auto i = 0; i < static_cast<int>(N); i++) {
                    size_t block_dim = (block_global_idx[i] + blocksize > global_dimensions[i])
                                           ? global_dimensions[i] - block_global_idx[i]
                                           : blocksize;
                    interp_end_idx[i] += block_dim - 1;
                    if (max_interp_level < ceil(log2(block_dim))) {
                        max_interp_level = static_cast<uint>(ceil(log2(block_dim)));
                    }
                }
                quant_inds[quant_index++] = quantizer.quantize_and_overwrite(*block, 0);

                for (uint level = max_interp_level; level > 0 && level <= max_interp_level; level--) {
                    uint stride_ip = 1U << (level - 1);
                    block_interpolation(data, block.get_global_index(), interp_end_idx, PB_predict_overwrite,
                                        interpolators[interpolator_id], direction_sequence_id, stride_ip);
                }
            }
        }*/
       // else{
        if (anchorStride == 0){
            quant_inds[quant_index++] = quantizer.quantize_and_overwrite(*data, 0);
        }
        else {
            build_anchor_grid(data);
            interpolation_level--;
        }
        for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
            double cur_eb = eb;
            //if (!conf.tuning){
                if (alpha < 0){
                    if (level >= 3){
                        cur_eb = eb * eb_ratio;
                    } else {
                        cur_eb = eb;
                    }
                }
                else if (alpha >= 1){              
                    double cur_ratio = pow(alpha, level - 1);
                    if (cur_ratio > beta){
                        cur_ratio = beta;
                    }            
                    cur_eb = eb / cur_ratio;
                }
            //}
            quantizer.set_eb(cur_eb);
            size_t stride = 1U << (level - 1);

            auto interp_block_size = blocksize * stride;


            auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(
                data, std::begin(global_dimensions), std::end(global_dimensions), interp_block_size, 0);

            auto inter_begin = inter_block_range->begin();
            auto inter_end = inter_block_range->end();
            auto cur_interpolator_id = level >=3 ? 0 : interpolator_id;
            for (auto block = inter_begin; block != inter_end; ++block) {
                auto end_idx = block.get_global_index();
                for (int i = 0; i < N; i++) {
                    end_idx[i] += interp_block_size;
                    if (end_idx[i] > global_dimensions[i] - 1) {
                        end_idx[i] = global_dimensions[i] - 1;
                    }
                }

                block_interpolation(data, block.get_global_index(), end_idx, PB_predict_overwrite,
                                    interpolators[cur_interpolator_id], direction_sequence_id, stride);
            }
        }
        //}
        quantizer.postcompress_data();
        return quant_inds_vec;
    }

    void save(uchar *&c) override {
        write(global_dimensions.data(), N, c);
        write(blocksize, c);
        write(interpolator_id, c);
        write(direction_sequence_id, c);
        write(anchorStride, c);
        write(alpha, c);
        write(beta, c);


        quantizer.save(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        read(global_dimensions.data(), N, c, remaining_length);
        read(blocksize, c, remaining_length);
        read(interpolator_id, c, remaining_length);
        read(direction_sequence_id, c, remaining_length);
        read(anchorStride, c, remaining_length);
        read(alpha, c, remaining_length);
        read(beta, c, remaining_length);

        quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

   private:
    enum PredictorBehavior { PB_predict_overwrite, PB_predict, PB_recover };

    void init() {
        quant_index = 0;
        assert(blocksize % 2 == 0 && "Interpolation block size should be even numbers");
        assert((anchorStride & anchorStride-1) == 0 && "Anchor stride should be 0 or 2's exponentials");
        num_elements = 1;
        interpolation_level = -1;
        for (int i = 0; i < N; i++) {
            if (interpolation_level < ceil(log2(global_dimensions[i]))) {
                interpolation_level = static_cast<uint>(ceil(log2(global_dimensions[i])));
            }
            num_elements *= global_dimensions[i];
        }

        if (anchorStride > 0){
            int max_interpolation_level=static_cast<uint>(log2(anchorStride))+1;
            if (max_interpolation_level<=interpolation_level){ 
                interpolation_level=max_interpolation_level;
            }
        }


        dimension_offsets[N - 1] = 1;
        for (int i = N - 2; i >= 0; i--) {
            dimension_offsets[i] = dimension_offsets[i + 1] * global_dimensions[i + 1];
        }

        dimension_sequences = std::vector<std::array<int, N>>();
        auto sequence = std::array<int, N>();
        for (int i = 0; i < N; i++) {
            sequence[i] = i;
        }
        do {
            dimension_sequences.push_back(sequence);
        } while (std::next_permutation(sequence.begin(), sequence.end()));
    }

    void build_anchor_grid(T *data){
        assert(anchorStride > 0);
        if (N == 1){
            for (size_t x = 0; x < global_dimensions[0]; x += anchorStride){
                quant_inds[quant_index++] = quantizer.force_save_unpred(*(data + x));
            }
        }
        else if (N == 2){
            for (size_t x = 0; x < global_dimensions[0]; x += anchorStride){
                for (size_t y = 0; y < global_dimensions[1]; y += anchorStride){
                    quant_inds[quant_index++] = quantizer.force_save_unpred(*(data + x * global_dimensions[1] + y));
                }
            }
        }
        else if (N == 3){
            for (size_t x = 0; x < global_dimensions[0]; x += anchorStride){
                for (size_t y = 0; y < global_dimensions[1]; y += anchorStride){
                    for(size_t z = 0; z < global_dimensions[2]; z += anchorStride){
                        quant_inds[quant_index++] =  quantizer.force_save_unpred(*(data + x * dimension_offsets[0] + y * dimension_offsets[1] + z));
                    }           
                }
            }
        }
        else if (N == 4){
            for (size_t x = 0; x < global_dimensions[0]; x += anchorStride){
                for (size_t y = 0; y < global_dimensions[1]; y += anchorStride){
                    for(size_t z = 0; z < global_dimensions[2]; z += anchorStride){
                        for(size_t w = 0; w < global_dimensions[3]; w += anchorStride){
                            quant_inds[quant_index++] = quantizer.force_save_unpred(*(data + x * dimension_offsets[0] + y * dimension_offsets[1] + z * dimension_offsets[2] + w));
                        }
                    }           
                }
            }
        }

    }

    void recover_anchor_grid(T *decData){
        assert(anchorStride > 0);

        if (N == 1){
            for (size_t x = 0; x < global_dimensions[0]; x += anchorStride){
                decData[x] = quantizer.recover_unpred();
                quant_index++; //not really necessary
            }
        }
        else if (N == 2){
            for (size_t x = 0; x < global_dimensions[0]; x += anchorStride){
                for (size_t y = 0; y < global_dimensions[1]; y += anchorStride){
                    decData[x * dimension_offsets[0] + y] = quantizer.recover_unpred();
                    quant_index++;
                    //decData[x * dimension_offsets[0] + y] = quantizer.recover(0, quant_inds[quant_index++]);
                }
            }
        }
        else if (N == 3){
            for (size_t x = 0; x < global_dimensions[0]; x += anchorStride){
                for (size_t y = 0; y < global_dimensions[1]; y += anchorStride){
                    for(size_t z = 0; z < global_dimensions[2]; z += anchorStride){
                        decData[x * dimension_offsets[0] + y * dimension_offsets[1] + z] = quantizer.recover_unpred();
                        quant_index++;
                    }           
                }
            }
        }
        else if (N == 4){
            for (size_t x = 0; x < global_dimensions[0]; x += anchorStride){
                for (size_t y = 0; y < global_dimensions[1]; y += anchorStride){
                    for(size_t z = 0; z < global_dimensions[2]; z += anchorStride){
                        for(size_t w = 0; w < global_dimensions[3]; w += anchorStride){
                            decData[x * dimension_offsets[0] + y * dimension_offsets[1] + z * dimension_offsets[2] + w] = quantizer.recover_unpred();
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

    inline void recover(size_t idx, T &d, T pred) { d = quantizer.recover(pred, quant_inds[quant_index++]); }

    double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride, const std::string &interp_func,
                                  const PredictorBehavior pb) {
        size_t n = (end - begin) / stride + 1;
        if (n <= 1) {
            return 0;
        }
        double predict_error = 0;

        size_t stride3x = 3 * stride;
        size_t stride5x = 5 * stride;
        if (interp_func == "linear" || n < 5) {
            if (pb == PB_predict_overwrite) {
                for (size_t i = 1; i + 1 < n; i += 2) {
                    T *d = data + begin + i * stride;
                    quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if (n < 4) {
                        quantize(d - data, *d, *(d - stride));
                    } else {
                        quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                    }
                }
            } else {
                for (size_t i = 1; i + 1 < n; i += 2) {
                    T *d = data + begin + i * stride;
                    recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if (n < 4) {
                        recover(d - data, *d, *(d - stride));
                    } else {
                        recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                    }
                }
            }
        } else if(interp_func == "cubic"){
            if (pb == PB_predict_overwrite) {
                T *d;
                size_t i;
                for (i = 3; i + 3 < n; i += 2) {
                    d = data + begin + i * stride;
                    quantize(d - data, *d,
                             interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                }
                d = data + begin + stride;
                quantize(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                d = data + begin + i * stride;
                quantize(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                if (n % 2 == 0) {
                    d = data + begin + (n - 1) * stride;
                    quantize(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                }

            } else {
                T *d;

                size_t i;
                for (i = 3; i + 3 < n; i += 2) {
                    d = data + begin + i * stride;
                    recover(d - data, *d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                }
                d = data + begin + stride;

                recover(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                d = data + begin + i * stride;
                recover(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                if (n % 2 == 0) {
                    d = data + begin + (n - 1) * stride;
                    recover(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                }
            }
        }
        else{
            if (pb == PB_predict_overwrite) {
                T *d;
                size_t i;
                for (i = 3; i + 3 < n; i += 2) {
                    d = data + begin + i * stride;
                    quantize(d - data, *d,
                             interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                }
                d = data + begin + stride;
                quantize(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                d = data + begin + i * stride;
                quantize(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                if (n % 2 == 0) {
                    d = data + begin + (n - 1) * stride;
                    quantize(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                }

            } else {
                T *d;

                size_t i;
                for (i = 3; i + 3 < n; i += 2) {
                    d = data + begin + i * stride;
                    recover(d - data, *d, interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                }
                d = data + begin + stride;

                recover(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                d = data + begin + i * stride;
                recover(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                if (n % 2 == 0) {
                    d = data + begin + (n - 1) * stride;
                    recover(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                }
            }
        }

        return predict_error;
    }
    
    double block_interpolation_1d_fastest_dim_first(
        T *data,
        const std::array<size_t, N> &begin_idx,
        const std::array<size_t, N> &end_idx,
        const size_t &direction,
        std::array<size_t, N> &steps,
        const size_t &math_stride,
        const std::string &interp_func,
        const PredictorBehavior pb) {
        for (size_t i = 0; i < N; i++) {
            if (end_idx[i] < begin_idx[i])
                return 0;
        }
        size_t math_begin_idx = begin_idx[direction], math_end_idx = end_idx[direction];
        size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;
        if (n <= 1) {
            return 0;
        }
        double predict_error = 0.0;
        size_t begin = 0;
        for (size_t i = 0; i < N; i++)
            begin += dimension_offsets[i] * begin_idx[i];
        size_t stride = math_stride * dimension_offsets[direction];
        std::array<size_t, N> begins, ends, strides;
        for (size_t i = 0; i < N; i++) {
            begins[i] = 0;
            ends[i] = end_idx[i] - begin_idx[i] + 1;
            strides[i] = dimension_offsets[i];
        }
        strides[direction] = stride;
        size_t stride2x = 2 * stride;
        if (pb == PB_predict_overwrite) {
            if (interp_func == "linear") {
                begins[direction] = 1;
                ends[direction] = n - 1;
                steps[direction] = 2;
                if constexpr (N == 1){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        T *d = data + begin + i * strides[0]; 
                        quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                        
                    }
                }
                else if constexpr (N == 2){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            T *d = data + begin + i * strides[0] + j * strides[1];
                            quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                        }
                    }
                }
                else if constexpr (N == 3){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                            }
                        }
                    }
                }
                else {
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                    T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                    quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                                }
                            }
                        }
                    }
                }
                if (n % 2 == 0) {
                    begins[direction] = n - 1;
                    ends[direction] = n;


                    if constexpr (N == 1){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            T *d = data + begin + i * strides[0]; 
                            if (n < 3)
                                quantize(d - data, *d, *(d - stride));
                            else
                                quantize(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                            
                        }
                    }
                    else if constexpr (N == 2){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                T *d = data + begin + i * strides[0] + j * strides[1];
                                if (n < 3)
                                    quantize(d - data, *d, *(d - stride));
                                else
                                    quantize(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                            }
                        }
                    }
                    else if constexpr (N == 3){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                    if (n < 3)
                                        quantize(d - data, *d, *(d - stride));
                                    else
                                        quantize(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                                }
                            }
                        }
                    }
                    else {
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                        T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                        if (n < 3)
                                            quantize(d - data, *d, *(d - stride));
                                        else
                                            quantize(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                                    }
                                }
                            }
                        }
                    }
                }
            } else if (interp_func == "cubic") {
                size_t stride3x = 3 * stride;
                T *d;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                if constexpr (N == 1){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        d = data + begin + i * strides[0]; 
                        quantize(d - data, *d,
                                interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        
                    }
                }
                else if constexpr (N == 2){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            d = data + begin + i * strides[0] + j * strides[1];
                            quantize(d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                    }
                }
                else if constexpr (N == 3){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                quantize(d - data, *d,
                                        interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                            }
                        }
                    }
                }
                else {
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                    d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                    quantize(d - data, *d,
                                            interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                }
                            }
                        }
                    }
                }
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    if constexpr (N == 1){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            d = data + begin + i * strides[0]; 
                            if (ii >= 3) {
                                if (ii + 3 < n)
                                    quantize(d - data, *d,
                                             interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                else if (ii + 1 < n)
                                    quantize(d - data, *d,
                                             interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                else
                                    quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                            } else {
                                if (ii + 3 < n)
                                    quantize(d - data, *d,
                                             interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                else if (ii + 1 < n)
                                    quantize(d - data, *d,
                                             interp_linear(*(d - stride), *(d + stride)));
                                else
                                    quantize(d - data, *d, *(d - stride));
                            }
                            
                        }
                    }
                    else if constexpr (N == 2){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                d = data + begin + i * strides[0] + j * strides[1];
                                if (ii >= 3) {
                                    if (ii + 3 < n)
                                        quantize(d - data, *d,
                                                 interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        quantize(d - data, *d,
                                                 interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                    else
                                        quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                } else {
                                    if (ii + 3 < n)
                                        quantize(d - data, *d,
                                                 interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        quantize(d - data, *d,
                                                 interp_linear(*(d - stride), *(d + stride)));
                                    else
                                        quantize(d - data, *d, *(d - stride));
                                }
                            }
                        }
                    }
                    else if constexpr (N == 3){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                    if (ii >= 3) {
                                        if (ii + 3 < n)
                                            quantize(d - data, *d,
                                                     interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                        else if (ii + 1 < n)
                                            quantize(d - data, *d,
                                                     interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                        else
                                            quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                    } else {
                                        if (ii + 3 < n)
                                            quantize(d - data, *d,
                                                     interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                        else if (ii + 1 < n)
                                            quantize(d - data, *d,
                                                     interp_linear(*(d - stride), *(d + stride)));
                                        else
                                            quantize(d - data, *d, *(d - stride));
                                    }
                                }
                            }
                        }
                    }
                    else {
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                        d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                        if (ii >= 3) {
                                            if (ii + 3 < n)
                                                quantize(d - data, *d,
                                                         interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                            else if (ii + 1 < n)
                                                quantize(d - data, *d,
                                                         interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                            else
                                                quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                        } else {
                                            if (ii + 3 < n)
                                                quantize(d - data, *d,
                                                         interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                            else if (ii + 1 < n)
                                                quantize(d - data, *d,
                                                         interp_linear(*(d - stride), *(d + stride)));
                                            else
                                                quantize(d - data, *d, *(d - stride));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                size_t stride3x = 3 * stride;
                T *d;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                if constexpr (N == 1){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        d = data + begin + i * strides[0]; 
                        quantize(d - data, *d,
                                interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        
                    }
                }
                else if constexpr (N == 2){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            d = data + begin + i * strides[0] + j * strides[1];
                            quantize(d - data, *d,
                                    interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                    }
                }
                else if constexpr (N == 3){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                quantize(d - data, *d,
                                        interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                            }
                        }
                    }
                }
                else {
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                    d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                    quantize(d - data, *d,
                                            interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                }
                            }
                        }
                    }
                }
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    if constexpr (N == 1){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            d = data + begin + i * strides[0]; 
                            if (ii >= 3) {
                                if (ii + 3 < n)
                                    quantize(d - data, *d,
                                             interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                else if (ii + 1 < n)
                                    quantize(d - data, *d,
                                             interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                else
                                    quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                            } else {
                                if (ii + 3 < n)
                                    quantize(d - data, *d,
                                             interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                else if (ii + 1 < n)
                                    quantize(d - data, *d,
                                             interp_linear(*(d - stride), *(d + stride)));
                                else
                                    quantize(d - data, *d, *(d - stride));
                            }
                            
                        }
                    }
                    else if constexpr (N == 2){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                d = data + begin + i * strides[0] + j * strides[1];
                                if (ii >= 3) {
                                    if (ii + 3 < n)
                                        quantize(d - data, *d,
                                                 interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        quantize(d - data, *d,
                                                 interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                    else
                                        quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                } else {
                                    if (ii + 3 < n)
                                        quantize(d - data, *d,
                                                 interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        quantize(d - data, *d,
                                                 interp_linear(*(d - stride), *(d + stride)));
                                    else
                                        quantize(d - data, *d, *(d - stride));
                                }
                            }
                        }
                    }
                    else if constexpr (N == 3){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                    if (ii >= 3) {
                                        if (ii + 3 < n)
                                            quantize(d - data, *d,
                                                     interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                        else if (ii + 1 < n)
                                            quantize(d - data, *d,
                                                     interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                        else
                                            quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                    } else {
                                        if (ii + 3 < n)
                                            quantize(d - data, *d,
                                                     interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                        else if (ii + 1 < n)
                                            quantize(d - data, *d,
                                                     interp_linear(*(d - stride), *(d + stride)));
                                        else
                                            quantize(d - data, *d, *(d - stride));
                                    }
                                }
                            }
                        }
                    }
                    else {
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                        d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                        if (ii >= 3) {
                                            if (ii + 3 < n)
                                                quantize(d - data, *d,
                                                         interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                            else if (ii + 1 < n)
                                                quantize(d - data, *d,
                                                         interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                            else
                                                quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                        } else {
                                            if (ii + 3 < n)
                                                quantize(d - data, *d,
                                                         interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                            else if (ii + 1 < n)
                                                quantize(d - data, *d,
                                                         interp_linear(*(d - stride), *(d + stride)));
                                            else
                                                quantize(d - data, *d, *(d - stride));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (interp_func == "linear") {
                begins[direction] = 1;
                ends[direction] = n - 1;
                steps[direction] = 2;
                if constexpr (N == 1){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        T *d = data + begin + i * strides[0]; 
                        recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                        
                    }
                }
                else if constexpr (N == 2){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            T *d = data + begin + i * strides[0] + j * strides[1];
                            recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                        }
                    }
                }
                else if constexpr (N == 3){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                            }
                        }
                    }
                }
                else {
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                    T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                    recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                                }
                            }
                        }
                    }
                }
                if (n % 2 == 0) {
                    begins[direction] = n - 1;
                    ends[direction] = n;


                    if constexpr (N == 1){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            T *d = data + begin + i * strides[0]; 
                            if (n < 3)
                                recover(d - data, *d, *(d - stride));
                            else
                                recover(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                            
                        }
                    }
                    else if constexpr (N == 2){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                T *d = data + begin + i * strides[0] + j * strides[1];
                                if (n < 3)
                                    recover(d - data, *d, *(d - stride));
                                else
                                    recover(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                            }
                        }
                    }
                    else if constexpr (N == 3){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                    if (n < 3)
                                        recover(d - data, *d, *(d - stride));
                                    else
                                        recover(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                                }
                            }
                        }
                    }
                    else {
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                        T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                        if (n < 3)
                                            recover(d - data, *d, *(d - stride));
                                        else
                                            recover(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                                    }
                                }
                            }
                        }
                    }   
                }
            } else if (interp_func == "cubic") {
                size_t stride3x = 3 * stride;
                T *d;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                if constexpr (N == 1){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        d = data + begin + i * strides[0]; 
                        recover(d - data, *d,
                                interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        
                    }
                }
                else if constexpr (N == 2){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            d = data + begin + i * strides[0] + j * strides[1];
                            recover(d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                    }
                }
                else if constexpr (N == 3){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                recover(d - data, *d,
                                        interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                            }
                        }
                    }
                }
                else {
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                    d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                    recover(d - data, *d,
                                            interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                }
                            }
                        }
                    }
                }
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    if constexpr (N == 1){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            d = data + begin + i * strides[0]; 
                            if (ii >= 3) {
                                if (ii + 3 < n)
                                    recover(d - data, *d,
                                             interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                else if (ii + 1 < n)
                                    recover(d - data, *d,
                                             interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                else
                                    recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                            } else {
                                if (ii + 3 < n)
                                    recover(d - data, *d,
                                             interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                else if (ii + 1 < n)
                                    recover(d - data, *d,
                                             interp_linear(*(d - stride), *(d + stride)));
                                else
                                    recover(d - data, *d, *(d - stride));
                            }
                            
                        }
                    }
                    else if constexpr (N == 2){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                d = data + begin + i * strides[0] + j * strides[1];
                                if (ii >= 3) {
                                    if (ii + 3 < n)
                                        recover(d - data, *d,
                                                 interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        recover(d - data, *d,
                                                 interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                    else
                                        recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                } else {
                                    if (ii + 3 < n)
                                        recover(d - data, *d,
                                                 interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        recover(d - data, *d,
                                                 interp_linear(*(d - stride), *(d + stride)));
                                    else
                                        recover(d - data, *d, *(d - stride));
                                }
                            }
                        }
                    }
                    else if constexpr (N == 3){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                    if (ii >= 3) {
                                        if (ii + 3 < n)
                                            recover(d - data, *d,
                                                     interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                        else if (ii + 1 < n)
                                            recover(d - data, *d,
                                                     interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                        else
                                            recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                    } else {
                                        if (ii + 3 < n)
                                            recover(d - data, *d,
                                                     interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                        else if (ii + 1 < n)
                                            recover(d - data, *d,
                                                     interp_linear(*(d - stride), *(d + stride)));
                                        else
                                            recover(d - data, *d, *(d - stride));
                                    }
                                }
                            }
                        }
                    }
                    else {
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                        d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                        if (ii >= 3) {
                                            if (ii + 3 < n)
                                                recover(d - data, *d,
                                                         interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                            else if (ii + 1 < n)
                                                recover(d - data, *d,
                                                         interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                            else
                                                recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                        } else {
                                            if (ii + 3 < n)
                                                recover(d - data, *d,
                                                         interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                            else if (ii + 1 < n)
                                                recover(d - data, *d,
                                                         interp_linear(*(d - stride), *(d + stride)));
                                            else
                                                recover(d - data, *d, *(d - stride));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                size_t stride3x = 3 * stride;
                T *d;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                if constexpr (N == 1){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        d = data + begin + i * strides[0]; 
                        recover(d - data, *d,
                                interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        
                    }
                }
                else if constexpr (N == 2){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            d = data + begin + i * strides[0] + j * strides[1];
                            recover(d - data, *d,
                                    interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                    }
                }
                else if constexpr (N == 3){
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                recover(d - data, *d,
                                        interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                            }
                        }
                    }
                }
                else {
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                    d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                    recover(d - data, *d,
                                            interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                }
                            }
                        }
                    }
                }
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    if constexpr (N == 1){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            d = data + begin + i * strides[0]; 
                            if (ii >= 3) {
                                if (ii + 3 < n)
                                    recover(d - data, *d,
                                             interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                else if (ii + 1 < n)
                                    recover(d - data, *d,
                                             interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                else
                                    recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                            } else {
                                if (ii + 3 < n)
                                    recover(d - data, *d,
                                             interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                else if (ii + 1 < n)
                                    recover(d - data, *d,
                                             interp_linear(*(d - stride), *(d + stride)));
                                else
                                    recover(d - data, *d, *(d - stride));
                            }
                            
                        }
                    }
                    else if constexpr (N == 2){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                d = data + begin + i * strides[0] + j * strides[1];
                                if (ii >= 3) {
                                    if (ii + 3 < n)
                                        recover(d - data, *d,
                                                 interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        recover(d - data, *d,
                                                 interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                    else
                                        recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                } else {
                                    if (ii + 3 < n)
                                        recover(d - data, *d,
                                                 interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        recover(d - data, *d,
                                                 interp_linear(*(d - stride), *(d + stride)));
                                    else
                                        recover(d - data, *d, *(d - stride));
                                }
                            }
                        }
                    }
                    else if constexpr (N == 3){
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                    if (ii >= 3) {
                                        if (ii + 3 < n)
                                            recover(d - data, *d,
                                                     interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                        else if (ii + 1 < n)
                                            recover(d - data, *d,
                                                     interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                        else
                                            recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                    } else {
                                        if (ii + 3 < n)
                                            recover(d - data, *d,
                                                     interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                        else if (ii + 1 < n)
                                            recover(d - data, *d,
                                                     interp_linear(*(d - stride), *(d + stride)));
                                        else
                                            recover(d - data, *d, *(d - stride));
                                    }
                                }
                            }
                        }
                    }
                    else {
                        for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                            for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                                for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                    for (size_t l = begins[3]; l < ends[3]; l += steps[3]) {
                                        d = data + begin + i * strides[0] + j * strides[1] + k * strides[2] + l * strides[3];
                                        if (ii >= 3) {
                                            if (ii + 3 < n)
                                                recover(d - data, *d,
                                                         interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                            else if (ii + 1 < n)
                                                recover(d - data, *d,
                                                         interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                            else
                                                recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                        } else {
                                            if (ii + 3 < n)
                                                recover(d - data, *d,
                                                         interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                            else if (ii + 1 < n)
                                                recover(d - data, *d,
                                                         interp_linear(*(d - stride), *(d + stride)));
                                            else
                                                recover(d - data, *d, *(d - stride));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return predict_error;
    }
    template <uint NN = N>
    typename std::enable_if<NN == 1, double>::type block_interpolation(T *data, std::array<size_t, N> begin,
                                                                       std::array<size_t, N> end,
                                                                       const PredictorBehavior pb,
                                                                       const std::string &interp_func,
                                                                       const int direction, size_t stride = 1) {
        return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb);
    }

    template <uint NN = N>
    typename std::enable_if<NN == 2, double>::type block_interpolation(T *data, std::array<size_t, N> begin,
                                                                       std::array<size_t, N> end,
                                                                       const PredictorBehavior pb,
                                                                       const std::string &interp_func,
                                                                       const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        const std::array<int, N> dims = dimension_sequences[direction];
        std::array<size_t, N> steps;
        std::array<size_t, N> begin_idx = begin, end_idx = end;
        steps[dims[0]] =1 ;
        begin_idx[dims[1]] = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
        steps[dims[1]] = stride2x;

        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[0], steps,
                                                            stride, interp_func, pb);

        begin_idx[dims[1]] = begin[dims[1]];
        begin_idx[dims[0]] = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
        steps[dims[0]] = stride;
       
        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[1], steps,
                                                            stride, interp_func, pb);
        return predict_error;
    }

    template <uint NN = N>
    typename std::enable_if<NN == 3, double>::type block_interpolation(T *data, std::array<size_t, N> begin,
                                                                       std::array<size_t, N> end,
                                                                       const PredictorBehavior pb,
                                                                       const std::string &interp_func,
                                                                       const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        const std::array<int, N> dims = dimension_sequences[direction];
        std::array<size_t, N> steps;
        std::array<size_t, N> begin_idx = begin, end_idx = end;
        steps[dims[0]] = 1;
        begin_idx[dims[1]] = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
        begin_idx[dims[2]] = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
        steps[dims[1]] = stride2x;
        steps[dims[2]] = stride2x;

        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[0], steps,
                                                            stride, interp_func, pb);
        
        begin_idx[dims[1]] = begin[dims[1]];
        begin_idx[dims[0]] = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
        steps[dims[0]] = stride;

        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[1], steps,
                                                            stride, interp_func, pb);

        
        begin_idx[dims[2]] = begin[dims[2]];
        begin_idx[dims[1]] = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
        steps[dims[1]] = stride;

        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[2], steps,
                                                            stride, interp_func, pb);
        return predict_error;
    }

    template <uint NN = N>
    typename std::enable_if<NN == 4, double>::type block_interpolation(T *data, std::array<size_t, N> begin,
                                                                       std::array<size_t, N> end,
                                                                       const PredictorBehavior pb,
                                                                       const std::string &interp_func,
                                                                       const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        max_error = 0;
        const std::array<int, N> dims = dimension_sequences[direction];
        std::array<size_t, N> steps;
        std::array<size_t, N> begin_idx = begin, end_idx = end;
        steps[dims[0]] = 1;
        begin_idx[dims[1]] = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
        begin_idx[dims[2]] = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
        begin_idx[dims[3]] = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
        steps[dims[1]] = stride2x;
        steps[dims[2]] = stride2x;
        steps[dims[3]] = stride2x;

        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[0], steps,
                                                            stride, interp_func, pb);
        
        begin_idx[dims[1]] = begin[dims[1]];
        begin_idx[dims[0]] = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
        steps[dims[0]] = stride;

        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[1], steps,
                                                            stride, interp_func, pb);

        
        begin_idx[dims[2]] = begin[dims[2]];
        begin_idx[dims[1]] = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
        steps[dims[1]] = stride;

        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[2], steps,
                                                            stride, interp_func, pb);

        begin_idx[dims[3]] = begin[dims[3]];
        begin_idx[dims[2]] = (begin[dims[2]] ? begin[dims[2]] + stride : 0);
        steps[dims[2]] = stride;

        predict_error += block_interpolation_1d_fastest_dim_first(data, begin_idx,
                                                            end_idx, dims[3], steps,
                                                            stride, interp_func, pb);
        return predict_error;
    }

    int interpolation_level = -1;
    uint blocksize;
    int interpolator_id;
    double eb_ratio = 0.5;
    std::vector<std::string> interpolators = {"linear", "cubic", "cubic_natural"};
    int *quant_inds;
    size_t quant_index = 0;
    double max_error;
    Quantizer quantizer;
    size_t num_elements;
    std::array<size_t, N> global_dimensions;
    std::array<size_t, N> dimension_offsets;
    std::vector<std::array<int, N>> dimension_sequences;
    int direction_sequence_id;

    int anchorStride = 0;
    double alpha = -1;
    double beta = -1;
};

template <class T, uint N, class Quantizer>
InterpolationDecomposition<T, N, Quantizer> make_decomposition_interpolation(const Config &conf, Quantizer quantizer) {
    return InterpolationDecomposition<T, N, Quantizer>(conf, quantizer);
}

}  // namespace SZ3

#endif
