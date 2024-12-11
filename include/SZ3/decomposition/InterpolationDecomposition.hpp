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
        //            lossless.postdecompress_data(buffer);
        double eb = quantizer.get_eb();

        *dec_data = quantizer.recover(0, this->quant_inds[quant_index++]);

        for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
            if (level >= 3) {
                quantizer.set_eb(eb * eb_ratio);
            } else {
                quantizer.set_eb(eb);
            }
            size_t stride = 1U << (level - 1);
            auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(
                dec_data, std::begin(global_dimensions), std::end(global_dimensions), stride * blocksize, 0);
            auto inter_begin = inter_block_range->begin();
            auto inter_end = inter_block_range->end();
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
        //            timer.stop("Interpolation Decompress");

        return dec_data;
    }

    // compress given the error bound
    std::vector<int> compress(const Config &conf, T *data) override {
        std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
        blocksize = 32;
        interpolator_id = conf.interpAlgo;
        direction_sequence_id = conf.interpDirection;

        init();

        std::vector<int> quant_inds_vec(num_elements);
        quant_inds = quant_inds_vec.data();

        double eb = quantizer.get_eb();

        //            quant_inds.push_back(quantizer.quantize_and_overwrite(*data, 0));
        quant_inds[quant_index++] = quantizer.quantize_and_overwrite(*data, 0);

        //            Timer timer;
        //            timer.start();

        for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
            if (level >= 3) {
                quantizer.set_eb(eb * eb_ratio);
            } else {
                quantizer.set_eb(eb);
            }
            size_t stride = 1U << (level - 1);

            auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(
                data, std::begin(global_dimensions), std::end(global_dimensions), blocksize * stride, 0);

            auto inter_begin = inter_block_range->begin();
            auto inter_end = inter_block_range->end();

            for (auto block = inter_begin; block != inter_end; ++block) {
                auto end_idx = block.get_global_index();
                for (int i = 0; i < N; i++) {
                    end_idx[i] += blocksize * stride;
                    if (end_idx[i] > global_dimensions[i] - 1) {
                        end_idx[i] = global_dimensions[i] - 1;
                    }
                }

                block_interpolation(data, block.get_global_index(), end_idx, PB_predict_overwrite,
                                    interpolators[interpolator_id], direction_sequence_id, stride);
            }
        }

        quantizer.postcompress_data();
        return quant_inds_vec;
    }

    void save(uchar *&c) override {
        write(global_dimensions.data(), N, c);
        write(blocksize, c);
        write(interpolator_id, c);
        write(direction_sequence_id, c);

        quantizer.save(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        read(global_dimensions.data(), N, c, remaining_length);
        read(blocksize, c, remaining_length);
        read(interpolator_id, c, remaining_length);
        read(direction_sequence_id, c, remaining_length);

        quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

   private:
    enum PredictorBehavior { PB_predict_overwrite, PB_predict, PB_recover };

    void init() {
        quant_index = 0;
        assert(blocksize % 2 == 0 && "Interpolation block size should be even numbers");
        num_elements = 1;
        interpolation_level = -1;
        for (int i = 0; i < N; i++) {
            if (interpolation_level < ceil(log2(global_dimensions[i]))) {
                interpolation_level = static_cast<uint>(ceil(log2(global_dimensions[i])));
            }
            num_elements *= global_dimensions[i];
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
        } else {
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
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
            size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
            predict_error += block_interpolation_1d(
                data, begin_offset, begin_offset + (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                stride * dimension_offsets[dims[0]], interp_func, pb);
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
            predict_error += block_interpolation_1d(
                data, begin_offset, begin_offset + (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                stride * dimension_offsets[dims[1]], interp_func, pb);
        }
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
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                      k * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(
                    data, begin_offset, begin_offset + (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                    stride * dimension_offsets[dims[0]], interp_func, pb);
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                      k * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(
                    data, begin_offset, begin_offset + (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                    stride * dimension_offsets[dims[1]], interp_func, pb);
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                      begin[dims[2]] * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(
                    data, begin_offset, begin_offset + (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
                    stride * dimension_offsets[dims[2]], interp_func, pb);
            }
        }
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
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0); t <= end[dims[3]]; t += stride2x) {
                    size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]] + t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(
                        data, begin_offset, begin_offset + (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                        stride * dimension_offsets[dims[0]], interp_func, pb);
                }
            }
        }
        max_error = 0;
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0); t <= end[dims[3]]; t += stride2x) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]] + t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(
                        data, begin_offset, begin_offset + (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                        stride * dimension_offsets[dims[1]], interp_func, pb);
                }
            }
        }
        max_error = 0;
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0); t <= end[dims[3]]; t += stride2x) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          begin[dims[2]] * dimension_offsets[dims[2]] + t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(
                        data, begin_offset, begin_offset + (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
                        stride * dimension_offsets[dims[2]], interp_func, pb);
                }
            }
        }

        max_error = 0;
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0); k <= end[dims[2]]; k += stride) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]] + begin[dims[3]] * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(
                        data, begin_offset, begin_offset + (end[dims[3]] - begin[dims[3]]) * dimension_offsets[dims[3]],
                        stride * dimension_offsets[dims[3]], interp_func, pb);
                }
            }
        }
        return predict_error;
    }

    int interpolation_level = -1;
    uint blocksize;
    int interpolator_id;
    double eb_ratio = 0.5;
    std::vector<std::string> interpolators = {"linear", "cubic"};
    int *quant_inds;
    size_t quant_index = 0;
    double max_error;
    Quantizer quantizer;
    size_t num_elements;
    std::array<size_t, N> global_dimensions;
    std::array<size_t, N> dimension_offsets;
    std::vector<std::array<int, N>> dimension_sequences;
    int direction_sequence_id;
};

template <class T, uint N, class Quantizer>
InterpolationDecomposition<T, N, Quantizer> make_decomposition_interpolation(const Config &conf, Quantizer quantizer) {
    return InterpolationDecomposition<T, N, Quantizer>(conf, quantizer);
}

}  // namespace SZ3

#endif
