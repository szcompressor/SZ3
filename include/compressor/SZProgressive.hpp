#ifndef _SZ_SZ_PROG_INTERPOLATION_V2_HPP
#define _SZ_SZ_PROG_INTERPOLATION_V2_HPP

#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/MemoryUtil.hpp"
#include "utils/Config.hpp"
#include "utils/FileUtil.hpp"
#include "utils/Interpolators.hpp"
#include "utils/Timer.hpp"
#include "def.hpp"
#include <cstring>
#include <cmath>

namespace SZ {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZProgressive {
    public:


        SZProgressive(Quantizer quantizer, Encoder encoder, Lossless lossless,
                      const std::array<size_t, N> dims,
                      int interpolator,
                      int direction,
                      size_t interp_dim_limit,
                      int level_independent_,
                      size_t block_size_,
                      int level_fill_) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                global_dimensions(dims),
                interpolators({"linear", "cubic"}),
                interpolator_id(interpolator), direction_sequence_id(direction),
                interp_dim_limit(interp_dim_limit),
                level_independent(level_independent_),
                block_size(block_size_),
                level_fill(level_fill_) {
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");

            assert(interp_dim_limit % 2 == 0 &&
                   "Interpolation dimension should be even numbers to avoid extrapolation");
            num_elements = 1;
            levels = -1;
            for (int i = 0; i < N; i++) {
                if (levels < ceil(log2(dims[i]))) {
                    levels = (uint) ceil(log2(dims[i]));
                }
                num_elements *= dims[i];
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


        T *decompress(uchar const *lossless_data, const std::vector<size_t> &lossless_size) {

            int lossless_id = 0;
            size_t remaining_length = lossless_size[lossless_id];
            uchar const *data_header = lossless_data;

            read(global_dimensions.data(), N, data_header, remaining_length);
            num_elements = std::accumulate(global_dimensions.begin(), global_dimensions.end(), (size_t) 1, std::multiplies<>());

            read(interp_dim_limit, data_header, remaining_length);
            lossless_data += lossless_size[lossless_id];
            lossless_id++;

            T *dec_data = new T[num_elements];
            size_t quant_inds_count = 0;

            for (uint level = levels; level > level_fill; level--) {
                Timer timer(true);
                size_t stride = 1U << (level - 1);

                if (level > level_independent) {
                    lossless_decode(lossless_data, lossless_size, lossless_id++);

                    if (level == levels) {
                        *dec_data = quantizer.recover(0, quant_inds[quant_index++]);
                    }

                    auto block_range = std::make_shared<
                            SZ::multi_dimensional_range<T, N>>(dec_data,
                                                               std::begin(global_dimensions), std::end(global_dimensions),
                                                               stride * interp_dim_limit, 0);
                    for (auto block = block_range->begin(); block != block_range->end(); ++block) {
                        auto end_idx = block.get_global_index();
                        for (int i = 0; i < N; i++) {
                            end_idx[i] += stride * interp_dim_limit;
                            if (end_idx[i] > global_dimensions[i] - 1) {
                                end_idx[i] = global_dimensions[i] - 1;
                            }
                        }
                        block_interpolation(dec_data, block.get_global_index(), end_idx, PB_recover,
                                            interpolators[interpolator_id], direction_sequence_id, stride, true);
                    }
                    quantizer.postdecompress_data();
                    std::cout << "Level = " << level << " , quant size = " << quant_inds.size() << ", Time = " << timer.stop() << std::endl;
                    quant_inds_count += quant_index;
                } else {
                    auto block_range = std::make_shared<
                            SZ::multi_dimensional_range<T, N>>(dec_data, std::begin(global_dimensions),
                                                               std::end(global_dimensions),
                                                               block_size, 0);
                    size_t nBlock = 0;
                    for (auto block = block_range->begin(); block != block_range->end(); ++block, nBlock++) {
                        lossless_decode(lossless_data, lossless_size, lossless_id++);

                        auto end_idx = block.get_global_index();
                        for (int i = 0; i < N; i++) {
                            end_idx[i] += block_size;
                            if (end_idx[i] > global_dimensions[i] - 1) {
                                end_idx[i] = global_dimensions[i] - 1;
                            }
                        }
                        auto first = block_range->begin();
                        *first = quantizer.recover(0, quant_inds[quant_index++]);

                        block_interpolation(dec_data, block.get_global_index(), end_idx, PB_recover,
                                            interpolators[interpolator_id], direction_sequence_id, stride, false);

                        quantizer.postdecompress_data();
                        quant_inds_count += quant_index;
                    }
                    std::cout << "Level = " << level << " , #blocks = " << nBlock <<", Time = "<< timer.stop() << std::endl;
                }

            }

            for (uint level = level_fill; level > 0; level--) {
                Timer timer(true);
                size_t stride = 1U << (level - 1);

                auto block_range = std::make_shared<
                        SZ::multi_dimensional_range<T, N>>(dec_data,
                                                           std::begin(global_dimensions), std::end(global_dimensions),
                                                           stride * interp_dim_limit, 0);
                for (auto block = block_range->begin(); block != block_range->end(); ++block) {
                    auto end_idx = block.get_global_index();
                    for (int i = 0; i < N; i++) {
                        end_idx[i] += stride * interp_dim_limit;
                        if (end_idx[i] > global_dimensions[i] - 1) {
                            end_idx[i] = global_dimensions[i] - 1;
                        }
                    }
                    block_interpolation(dec_data, block.get_global_index(), end_idx, PB_fill,
                                        interpolators[interpolator_id], direction_sequence_id, stride, true);
                }
                std::cout << "Fill Level = " << level << " " << std::endl;
            }
//            assert(num_elements == quant_inds_count);
            return dec_data;
        }

        uchar *compress(T *data, std::vector<size_t> &compressed_size) {
            return compress(data, compressed_size, false);
        }

        // compress given the error bound
        uchar *compress(T *data, std::vector<size_t> &lossless_size, bool deleteData) {
            quant_inds.reserve(num_elements);
//            quant_inds.resize(num_elements);
            size_t interp_compressed_size = 0;
//            debug.resize(num_elements, 0);
//            preds.resize(num_elements, 0);
            size_t quant_inds_total = 0;

            T eb = quantizer.get_eb();
            std::cout << "Absolute error bound = " << eb << std::endl;
//            quantizer.set_eb(eb * eb_ratio);


            Timer timer(true);
            uchar *lossless_data = (num_elements < 1000000) ?
                                   new uchar[4 * num_elements * sizeof(T)] :
                                   new uchar[size_t(1.2 * num_elements) * sizeof(T)];
            uchar *lossless_data_pos = lossless_data;

            write(global_dimensions.data(), N, lossless_data_pos);
            write(interp_dim_limit, lossless_data_pos);
            lossless_size.push_back(lossless_data_pos - lossless_data);

            for (uint level = levels; level > 0; level--) {
                quantizer.set_eb((level >= 3) ? eb * eb_ratio : eb);

                uint stride = 1U << (level - 1);
//                std::cout << "Level = " << level << ", stride = " << stride << std::endl;

                if (level > level_independent) {
                    timer.start();
                    if (level == levels) {
                        quant_inds.push_back(quantizer.quantize_and_overwrite(*data, 0));
                    }
                    auto block_range = std::make_shared<
                            SZ::multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                               std::end(global_dimensions),
                                                               interp_dim_limit * stride, 0);

                    for (auto block = block_range->begin(); block != block_range->end(); ++block) {
                        auto end_idx = block.get_global_index();
                        for (int i = 0; i < N; i++) {
                            end_idx[i] += interp_dim_limit * stride;
                            if (end_idx[i] > global_dimensions[i] - 1) {
                                end_idx[i] = global_dimensions[i] - 1;
                            }
                        }

                        block_interpolation(data, block.get_global_index(), end_idx, PB_predict_overwrite,
                                            interpolators[interpolator_id], direction_sequence_id, stride, true);

                    }
                    printf("level = %d , quant size = %lu , time=%.3f\n", level, quant_inds.size(), timer.stop());

                    quant_inds_total += quant_inds.size();
                    encode_lossless(lossless_data_pos, lossless_size);

                } else {
                    timer.start();
                    lossless_time = 0;
                    encode_time = 0;
                    auto block_range = std::make_shared<
                            SZ::multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                               std::end(global_dimensions),
                                                               block_size, 0);
                    size_t nBlock = 0;
                    for (auto block = block_range->begin(); block != block_range->end(); ++block, nBlock++) {
                        auto end_idx = block.get_global_index();
                        for (int i = 0; i < N; i++) {
                            end_idx[i] += block_size;
                            if (end_idx[i] > global_dimensions[i] - 1) {
                                end_idx[i] = global_dimensions[i] - 1;
                            }
                        }
                        auto first = block_range->begin();
                        quant_inds.push_back(quantizer.quantize_and_overwrite(*first, 0));
                        block_interpolation(data, block.get_global_index(), end_idx, PB_predict_overwrite,
                                            interpolators[interpolator_id], direction_sequence_id, stride, false);
                        quant_inds_total += quant_inds.size();
                        encode_lossless(lossless_data_pos, lossless_size);
                    }
                    printf("level = %d , #block = %lu, time=%.3f, encode_time=%.3f lossless_time=%.3f\n", level, nBlock, timer.stop(), encode_time,
                           lossless_time);

                }
            }
            if (deleteData) {
                delete[]data;
            }
            std::cout << "total element = " << num_elements << ", quantization element = " << quant_inds_total << std::endl;
            assert(quant_inds_total >= num_elements);

//            compressed_size += interp_compressed_size;
            return lossless_data;
        }

    private:

        void lossless_decode(uchar const *&lossless_data_pos, const std::vector<size_t> &lossless_size, int lossless_id) {

            size_t remaining_length = lossless_size[lossless_id];

            uchar *compressed_data = lossless.decompress(lossless_data_pos, remaining_length);
            uchar const *compressed_data_pos = compressed_data;

            quantizer.load(compressed_data_pos, remaining_length);
            double eb = quantizer.get_eb();

            size_t quant_size;
            read(quant_size, compressed_data_pos, remaining_length);
            //                printf("%lu\n", quant_size);
            if (quant_size < 128) {
                quant_inds.resize(quant_size);
                read(quant_inds.data(), quant_size, compressed_data_pos, remaining_length);
            } else {
                encoder.load(compressed_data_pos, remaining_length);
                quant_inds = encoder.decode(compressed_data_pos, quant_size);
                encoder.postprocess_decode();
            }
            quant_index = 0;

            lossless.postdecompress_data(compressed_data);
            lossless_data_pos += lossless_size[lossless_id];
        }

        void encode_lossless(uchar *&lossless_data_pos, std::vector<size_t> &lossless_size) {
            uchar *compressed_data = (quant_inds.size() < 1000000) ?
                                     new uchar[10 * quant_inds.size() * sizeof(T)] :
                                     new uchar[size_t(1.2 * quant_inds.size()) * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;

            quantizer.save(compressed_data_pos);
            quantizer.postcompress_data();
            quantizer.clear();
            write((size_t) quant_inds.size(), compressed_data_pos);
            if (quant_inds.size() < 128) {
                write(quant_inds.data(), quant_inds.size(), compressed_data_pos);
            } else {
                encoder.preprocess_encode(quant_inds, 0);
                encoder.save(compressed_data_pos);
                encoder.encode(quant_inds, compressed_data_pos);
                encoder.postprocess_encode();
            }

            size_t size = 0;
            uchar *lossless_data_cur_level = lossless.compress(compressed_data,
                                                               compressed_data_pos - compressed_data,
                                                               size);
            lossless.postcompress_data(compressed_data);
            memcpy(lossless_data_pos, lossless_data_cur_level, size);
            lossless_data_pos += size;
            lossless_size.push_back(size);
            delete[]lossless_data_cur_level;

            quant_inds.clear();

        }

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover, PB_fill
        };

        inline void quantize1(size_t idx, T &d, T pred) {
            if (idx >= 2 * 449 * 449 * 235 && idx < 3 * 449 * 449 * 235 &&
                fabs(d - pred) > max_error) {
                max_error = fabs(d - pred);
            }
            auto quant = quantizer.quantize_and_overwrite(d, pred);
//            quant_inds.push_back(quant);
            quant_inds[idx] = quant;
        }

        inline void quantize(size_t idx, T &d, T pred) {
//            preds[idx] = pred;
//            quant_inds[idx] = quantizer.quantize_and_overwrite(d, pred);
            quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
        }

        inline void recover(size_t idx, T &d, T pred) {
            d = quantizer.recover(pred, quant_inds[quant_index++]);
//            if (fabs(ori_data[idx] - d) > 1e-3) {
//                printf("%lu %.5f %.5f %d %lu\n", idx, ori_data[idx], d, quant_inds[quant_index-1],quant_index-1);
//            }
        };

        inline void fill(size_t idx, T &d, T pred) {
            d = pred;
        };

        double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride,
                                      const std::string &interp_func,
                                      const PredictorBehavior pb) {
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }
            double predict_error = 0;


            size_t stride3x = 3 * stride;
            size_t stride5x = 5 * stride;
//            printf("stride %lu %lu %lu\n", stride, stride3x, stride5x);
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
                } else if (pb == PB_recover) {
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
                } else {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        fill(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
//                        fill(d - data, *d, *(d - stride));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            fill(d - data, *d, *(d - stride));
                        } else {
                            fill(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
//                            fill(d - data, *d, *(d - stride));
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

                } else if (pb == PB_recover) {
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
                } else {
                    T *d;

                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        fill(d - data, *d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    d = data + begin + stride;

                    fill(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    d = data + begin + i * stride;
                    fill(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        fill(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    }
                }
            }

            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 1, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride, bool overlap) {
            return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb);
        }

        template<uint NN = N>
        typename std::enable_if<NN == 2, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride, bool overlap) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = ((overlap && begin[dims[1]]) ? begin[dims[1]] + stride2x : begin[dims[1]]); j <= end[dims[1]]; j += stride2x) {
                size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset + (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                        stride * dimension_offsets[dims[0]], interp_func, pb);
            }
            for (size_t i = ((overlap && begin[dims[0]]) ? begin[dims[0]] + stride : begin[dims[0]]); i <= end[dims[0]]; i += stride) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset + (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                        stride * dimension_offsets[dims[1]], interp_func, pb);
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 3, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride, bool overlap) {
            double predict_error = 0;
            size_t stride2x = stride * 2;

            std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = ((overlap && begin[dims[1]]) ? begin[dims[1]] + stride2x : begin[dims[1]]); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = ((overlap && begin[dims[2]]) ? begin[dims[2]] + stride2x : begin[dims[2]]); k <= end[dims[2]]; k += stride2x) {
                    size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]]
                                          + k * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset + (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[0]], interp_func, pb);
                }
            }


            for (size_t i = ((overlap && begin[dims[0]]) ? begin[dims[0]] + stride : begin[dims[0]]); i <= end[dims[0]]; i += stride) {
                for (size_t k = ((overlap && begin[dims[2]]) ? begin[dims[2]] + stride2x : begin[dims[2]]); k <= end[dims[2]]; k += stride2x) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset + (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb);
                }
            }

            for (size_t i = ((overlap && begin[dims[0]]) ? begin[dims[0]] + stride : begin[dims[0]]); i <= end[dims[0]]; i += stride) {
                for (size_t j = ((overlap && begin[dims[1]]) ? begin[dims[1]] + stride : begin[dims[1]]); j <= end[dims[1]]; j += stride) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          begin[dims[2]] * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset + (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
                                                            stride * dimension_offsets[dims[2]], interp_func, pb);
                }
            }


            return predict_error;
        }


        template<uint NN = N>
        typename std::enable_if<NN == 4, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride, bool overlap) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            max_error = 0;
            std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                stride * dimension_offsets[dims[0]], interp_func, pb);
                    }
                }
            }
//            printf("%.8f ", max_error);
            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[1]], interp_func, pb);
                    }
                }
            }
//            printf("%.8f ", max_error);
            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              begin[dims[2]] * dimension_offsets[dims[2]] +
                                              t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[2]], interp_func, pb);
                    }
                }
            }

//            printf("%.8f ", max_error);
            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0); k <= end[dims[2]]; k += stride) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                begin[dims[3]] * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[3]] - begin[dims[3]]) *
                                                                dimension_offsets[dims[3]],
                                                                stride * dimension_offsets[dims[3]], interp_func, pb);
                    }
                }
            }
//            printf("%.8f \n", max_error);
            return predict_error;
        }

        std::vector<T> ori_data;
        double encode_time = 0, lossless_time = 0;
        int levels = -1, level_independent = -1, level_fill = 0;
        size_t interp_dim_limit, block_size;
        int interpolator_id;
        double eb_ratio = 0.5;
        std::vector<std::string> interpolators;
        std::vector<int> quant_inds;
        size_t quant_index = 0; // for decompress
//        std::vector<int> debug;
//        std::vector<T> preds;
        double max_error;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
        std::array<size_t, N> dimension_offsets;
        std::vector<std::array<int, N>> dimension_sequences;
        int direction_sequence_id;
    };


};


#endif

