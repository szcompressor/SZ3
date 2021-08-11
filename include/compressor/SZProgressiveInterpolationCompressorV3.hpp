#ifndef _SZ_SZ_PROG_INTERPOLATION_V3_HPP
#define _SZ_SZ_PROG_INTERPOLATION_V3_HPP

#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/MemoryUtil.hpp"
#include "utils/Config.hpp"
#include "utils/FileUtil.h"
#include "utils/Interpolators.hpp"
#include "utils/Timer.hpp"
#include "def.hpp"
#include <cstring>
#include <cmath>

namespace SZ {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZProgressiveInterpolationCompressorV3 {
    public:


        SZProgressiveInterpolationCompressorV3(Quantizer quantizer, Encoder encoder, Lossless lossless,
                                               const std::array<size_t, N> dims,
                                               int interpolator,
                                               int direction,
                                               size_t interp_dim_limit,
                                               size_t interp_block_size,
                                               int level_fill_) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                global_dimensions(dims),
                interpolators({"linear", "cubic"}),
                interpolator_id(interpolator), direction_sequence_id(direction),
                interp_block_limit(interp_dim_limit),
                level_fill(level_fill_) {
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");

            assert(interp_dim_limit % 2 == 0 &&
                   "Interpolation dimension should be even numbers to avoid extrapolation");
            level_independent = 2;
            core_stride = round(pow(2, level_independent));
            core_blocksize = interp_block_size / core_stride;
            block_size = core_blocksize * core_stride;
            printf("core_stride = %lu\ncore_blocksize = %lu\nblock_size = %lu\n", core_stride, core_blocksize, block_size);

            printf("core_dimension = ");
            levels = -1;
            for (int i = 0; i < N; i++) {
                if (levels < ceil(log2(dims[i]))) {
                    levels = (uint) ceil(log2(dims[i]));
                }
                core_dimensions[i] = ceil(1.0 * dims[i] / core_stride);
                printf("%lu ", core_dimensions[i]);
            }
            printf("\nlevels = %d\n", levels);


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
            read(core_dimensions.data(), N, data_header, remaining_length);
            read(interp_block_limit, data_header, remaining_length);
            size_t num_elements = std::accumulate(global_dimensions.begin(), global_dimensions.end(), 1, std::multiplies<T>());
            size_t core_num_elements = std::accumulate(core_dimensions.begin(), core_dimensions.end(), 1, std::multiplies<>());
            size_t block_num_elements = round(pow(block_size + 1, N));

            T *dec_data = new T[num_elements];
            std::vector<T> core_data(core_num_elements);

            read(core_data[0], data_header, remaining_length);
            lossless_data += lossless_size[lossless_id++];

            size_t quant_inds_count = 1;
            Timer timer;
            {
                timer.start();
                dimension_offsets[N - 1] = 1;
                for (int i = N - 2; i >= 0; i--) {
                    dimension_offsets[i] = dimension_offsets[i + 1] * core_dimensions[i + 1];
                }
                for (uint level = levels; level > level_independent; level--) {
                    uint stride = 1U << ((level - level_independent) - 1);

                    if (level > level_fill) {
                        lossless_decode(lossless_data, lossless_size, lossless_id++);
                    }

                    auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                            core_data.data(), std::begin(core_dimensions), std::end(core_dimensions), stride * interp_block_limit, 0);
                    for (auto block = block_range->begin(); block != block_range->end(); ++block) {
                        auto end_idx = block.get_global_index();
                        for (int i = 0; i < N; i++) {
                            end_idx[i] += stride * interp_block_limit;
                            if (end_idx[i] > core_dimensions[i] - 1) {
                                end_idx[i] = core_dimensions[i] - 1;
                            }
                        }
                        block_interpolation(core_data.data(), block.get_global_index(), end_idx,
                                            (level > level_fill ? PB_recover : PB_fill),
                                            interpolators[interpolator_id], direction_sequence_id, stride, true);
                    }
                    quantizer.postdecompress_data();
                    std::cout << "Level = " << level << " , quant size = " << quant_inds.size()
                              << " , time = " << timer.stop() << std::endl;
                    quant_inds_count += quant_index;
                }
            }

            {
                std::vector<T> block_data(block_num_elements);
                auto global_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        dec_data, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
                auto global_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        dec_data, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

                auto core_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        core_data.data(), std::begin(core_dimensions), std::end(core_dimensions), core_blocksize, 0);
                auto core_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        core_data.data(), std::begin(core_dimensions), std::end(core_dimensions), 1, 0);

                std::vector<size_t> quant_size(level_independent + 1, 0);
                double comptime = 0, copytime1 = 0, copytime2 = 0;

                size_t nBlock = 0;
                std::array<size_t, N> block_dims, core_block_dims, block_start_idx, block_end_idx;
                block_start_idx.fill(0);
                auto block = global_range_inter->begin();
                auto core_block = core_range_inter->begin();
                for (; block != global_range_inter->end(); ++block, nBlock++, ++core_block) {
                    for (int i = 0; i < N; i++) {
                        block_dims[i] = std::min(global_dimensions[i] - block.get_local_index(i) * block_size, block_size + 1);
                        core_block_dims[i] = std::min(core_dimensions[i] - core_block.get_local_index(i) * core_blocksize, core_blocksize + 1);
                        block_end_idx[i] = block_dims[i] - 1;
                    }
                    dimension_offsets[N - 1] = 1;
                    for (int i = N - 2; i >= 0; i--) {
                        dimension_offsets[i] = dimension_offsets[i + 1] * block_dims[i + 1];
                    }

                    {
                        timer.start();
                        auto block_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                                block_data.data(), std::begin(block_dims), std::end(block_dims), core_stride, 0);
                        core_range->set_dimensions(core_block_dims.begin(), core_block_dims.end());
                        core_range->set_offsets(core_block.get_offset());
                        core_range->set_starting_position(core_block.get_local_index());

                        auto block_iter = block_range_inter->begin();
                        auto core_iter = core_range->begin();
                        for (; core_iter != core_range->end(); ++core_iter, ++block_iter) {
                            *block_iter = *core_iter;
                        }
                        copytime1 += timer.stop();
                    }

                    timer.start();
                    for (uint level = level_independent; level > level_fill; level--) {
                        size_t stride = 1U << (level - 1);
                        lossless_decode(lossless_data, lossless_size, lossless_id++);

                        block_interpolation(block_data.data(), block_start_idx, block_end_idx, PB_recover,
                                            interpolators[interpolator_id], direction_sequence_id, stride, false);

                        quantizer.postdecompress_data();
                        quant_inds_count += quant_index;
                        quant_size[level] += quant_index;
                    }
                    for (uint level = level_fill; level > 0; level--) {
                        lossless_data += lossless_size[lossless_id++];
                        size_t stride = 1U << (level - 1);
                        block_interpolation(block_data.data(), block_start_idx, block_end_idx, PB_fill,
                                            interpolators[interpolator_id], direction_sequence_id, stride, false);
                    }
                    comptime += timer.stop();

                    {
                        timer.start();
                        auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                                block_data.data(), std::begin(block_dims), std::end(block_dims), 1, 0);
                        global_range->set_dimensions(block_dims.begin(), block_dims.end());
                        global_range->set_offsets(block.get_offset());
                        global_range->set_starting_position(block.get_local_index());
                        auto block_iter = block_range->begin();
                        auto global_iter = global_range->begin();
                        for (; global_iter != global_range->end(); ++global_iter, ++block_iter) {
                            *global_iter = *block_iter;
                        }
                        copytime2 += timer.stop();
                    }
                }
                for (uint level = level_independent; level > 0; level--) {
                    printf("level = %d , quant size = %lu\n", level, quant_size[level]);
                }
                printf("level %d to 1, comptime = %.3f, copytime = %.3f %.3f\n", level_independent, comptime, copytime1, copytime2);
            }

            return dec_data;
        }

        uchar *compress(T *data, std::vector<size_t> &compressed_size) {
            return compress(data, compressed_size, false);
        }

        // compress given the error bound
        uchar *compress(T *data, std::vector<size_t> &lossless_size, bool deleteData) {
            size_t num_elements = std::accumulate(global_dimensions.begin(), global_dimensions.end(), 1, std::multiplies<>());
            size_t core_num_elements = std::accumulate(core_dimensions.begin(), core_dimensions.end(), 1, std::multiplies<>());
            size_t block_num_elements = round(pow(block_size + 1, N));
            size_t quant_inds_total = 1;
            T eb = quantizer.get_eb();
            std::cout << "Absolute error bound = " << eb << std::endl;

            std::vector<T> core_data(core_num_elements);
            uchar *lossless_data = (num_elements < 1000000) ?
                                   new uchar[4 * num_elements * sizeof(T)] :
                                   new uchar[size_t(1.2 * num_elements) * sizeof(T)];
            uchar *lossless_data_pos = lossless_data;
            write(global_dimensions.data(), N, lossless_data_pos);
            write(core_dimensions.data(), N, lossless_data_pos);
            write(interp_block_limit, lossless_data_pos);
            write(data[0], lossless_data_pos);
            lossless_size.push_back(lossless_data_pos - lossless_data);
            Timer timer;

            {

                auto core_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        core_data.data(), std::begin(core_dimensions), std::end(core_dimensions), 1, 0);
                auto global_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        data, std::begin(global_dimensions), std::end(global_dimensions), core_stride, 0);

                auto core_iter = core_range->begin();
                auto global_iter = global_range_inter->begin();
                for (; core_iter != core_range->end(); ++core_iter, ++global_iter) {
                    *core_iter = *global_iter;
                }
//                timer.stop("prepare core");
            }
            {
                quant_inds.reserve(core_num_elements);
                //quantizer.set_eb(eb * eb_ratio);

                dimension_offsets[N - 1] = 1;
                for (int i = N - 2; i >= 0; i--) {
                    dimension_offsets[i] = dimension_offsets[i + 1] * core_dimensions[i + 1];
                }

                for (uint level = levels; level > level_independent; level--) {
                    timer.start();
                    quantizer.set_eb((level >= 3) ? eb * eb_ratio : eb);
                    uint stride = 1U << ((level - level_independent) - 1);
                    //                std::cout << "Level = " << level << ", stride = " << stride << std::endl;

                    //                if (level == levels) {
                    //                    quant_inds.push_back(quantizer.quantize_and_overwrite(*core_data, 0));
                    //                }
                    auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                            core_data.data(), std::begin(core_dimensions), std::end(core_dimensions), interp_block_limit * stride, 0);

                    for (auto block = block_range->begin(); block != block_range->end(); ++block) {
                        auto end_idx = block.get_global_index();
                        for (int i = 0; i < N; i++) {
                            end_idx[i] += interp_block_limit * stride;
                            if (end_idx[i] > core_dimensions[i] - 1) {
                                end_idx[i] = core_dimensions[i] - 1;
                            }
                        }

                        block_interpolation(core_data.data(), block.get_global_index(), end_idx, PB_predict_overwrite,
                                            interpolators[interpolator_id], direction_sequence_id, stride, true);

                    }
                    size_t quantsize = quant_inds.size();
                    quant_inds_total += quantsize;
                    encode_lossless(lossless_data_pos, lossless_size);
                    printf("level = %d , quant size = %lu , time=%.3f\n", level, quantsize, timer.stop());
                }
            }
//            std::cout << "quantization element = " << quant_inds_total << std::endl;

            {

                quant_inds.clear();
                quant_inds.reserve(block_num_elements);

                double comptime = 0, copytime1 = 0, copytime2 = 0;

                std::vector<T> block_data(block_num_elements);
                auto global_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        data, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
                auto global_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        data, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

                auto core_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        core_data.data(), std::begin(core_dimensions), std::end(core_dimensions), core_blocksize, 0);
                auto core_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        core_data.data(), std::begin(core_dimensions), std::end(core_dimensions), 1, 0);

                std::vector<size_t> quant_size(level_independent + 1, 0);

                size_t nBlock = 0;
                std::array<size_t, N> block_dims, core_block_dims, block_start_idx, block_end_idx;
                block_start_idx.fill(0);
                auto block = global_range_inter->begin();
                auto core_block = core_range_inter->begin();
                for (; block != global_range_inter->end(); ++block, nBlock++, ++core_block) {
//                    block.print();
                    for (int i = 0; i < N; i++) {
                        block_dims[i] = std::min(global_dimensions[i] - block.get_local_index(i) * block_size, block_size + 1);
                        core_block_dims[i] = std::min(core_dimensions[i] - core_block.get_local_index(i) * core_blocksize, core_blocksize + 1);
                        block_end_idx[i] = block_dims[i] - 1;
                    }
                    dimension_offsets[N - 1] = 1;
                    for (int i = N - 2; i >= 0; i--) {
                        dimension_offsets[i] = dimension_offsets[i + 1] * block_dims[i + 1];
                    }


                    {
                        timer.start();
                        auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                                block_data.data(), std::begin(block_dims), std::end(block_dims), 1, 0);
                        global_range->set_dimensions(block_dims.begin(), block_dims.end());
                        global_range->set_offsets(block.get_offset());
                        global_range->set_starting_position(block.get_local_index());
                        auto block_iter = block_range->begin();
                        auto global_iter = global_range->begin();
                        for (; global_iter != global_range->end(); ++global_iter, ++block_iter) {
                            *block_iter = *global_iter;
                        }
                        copytime1 += timer.stop();
                    }
                    {
                        timer.start();
                        auto block_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                                block_data.data(), std::begin(block_dims), std::end(block_dims), core_stride, 0);
                        core_range->set_dimensions(core_block_dims.begin(), core_block_dims.end());
                        core_range->set_offsets(core_block.get_offset());
                        core_range->set_starting_position(core_block.get_local_index());

                        auto block_iter = block_range_inter->begin();
                        auto core_iter = core_range->begin();
                        for (; core_iter != core_range->end(); ++core_iter, ++block_iter) {
                            *block_iter = *core_iter;
                        }
                        copytime2 += timer.stop();
                    }

                    timer.start();
                    for (uint level = level_independent; level > 0; level--) {
                        quantizer.set_eb((level >= 3) ? eb * eb_ratio : eb);
                        uint stride = 1U << (level - 1);
                        block_interpolation(block_data.data(), block_start_idx, block_end_idx, PB_predict_overwrite,
                                            interpolators[interpolator_id], direction_sequence_id, stride, false);
                        quant_inds_total += quant_inds.size();
                        quant_size[level] += quant_inds.size();

                        encode_lossless(lossless_data_pos, lossless_size);

//                        printf("level = %d , #block = %lu, quant_bins=%lu, time=%.3f, encode_time=%.3f lossless_time=%.3f\n",
//                               level, nBlock, tmp, timer.stop(), encode_time, lossless_time);
                    }
                    comptime += timer.stop();

                }
                for (uint level = level_independent; level > 0; level--) {
                    printf("level = %d , quant size = %lu\n", level, quant_size[level]);
                }
                printf("level %d to 1, comptime = %.3f, copytime = %.3f %.3f\n", level_independent, comptime, copytime1, copytime2);
            }


            if (deleteData) {
                delete[]data;
            }
            std::cout << "total element = " << num_elements << ", quantization element = " << quant_inds_total << std::endl;
            assert(quant_inds_total >= num_elements);
//            timer.stop("Predition & Quantization");


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
                int min;
                read(min, compressed_data_pos, remaining_length);
                encoder.load(compressed_data_pos, remaining_length);
                quant_inds = encoder.decode(compressed_data_pos, quant_size);
                for (auto &q:quant_inds) {
                    q += min;
                }
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
            write((size_t) quant_inds.size(), compressed_data_pos);
            if (quant_inds.size() < 128) {
                write(quant_inds.data(), quant_inds.size(), compressed_data_pos);
            } else {
                auto mm = std::minmax_element(quant_inds.begin(), quant_inds.end());
                int max = *mm.second, min = *mm.first;
//                printf("quant_max=%d, quant_min=%d\n", max, min);
                for (auto &q:quant_inds) {
                    q -= min;
                }
//                encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
                encoder.preprocess_encode(quant_inds, max - min + 1);

                write(min, compressed_data_pos);
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

        inline void recover(T &d, T pred) {
            d = quantizer.recover(pred, quant_inds[quant_index++]);
        };

        inline void fill(T &d, T pred) {
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
                        recover(*d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            recover(*d, *(d - stride));
                        } else {
                            recover(*d, interp_linear1(*(d - stride3x), *(d - stride)));
                        }
                    }
                } else {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        fill(*d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            fill(*d, *(d - stride));
                        } else {
                            fill(*d, interp_linear1(*(d - stride3x), *(d - stride)));
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
                        recover(*d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    d = data + begin + stride;

                    recover(*d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    d = data + begin + i * stride;
                    recover(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        recover(*d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    }
                } else {
                    T *d;

                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        fill(*d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    d = data + begin + stride;

                    fill(*d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    d = data + begin + i * stride;
                    fill(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        fill(*d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
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
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
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

        int levels = -1, level_independent = -1, level_fill = 0;
        size_t interp_block_limit, block_size, core_blocksize, core_stride;
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
        std::array<size_t, N> global_dimensions, core_dimensions;
        std::array<size_t, N> dimension_offsets;
        std::vector<std::array<int, N>> dimension_sequences;
        int direction_sequence_id;
    };


};


#endif

