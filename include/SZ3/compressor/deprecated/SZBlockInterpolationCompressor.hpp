#ifndef _SZ_BLOCK_INTERPOLATION_COMPRESSOR_HPP
#define _SZ_BLOCK_INTERPOLATION_COMPRESSOR_HPP

#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/def.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include <cstring>
#include <cmath>

namespace SZ {
    template<class T, uint N, class Predictor, class Quantizer, class Encoder, class Lossless>
    class SZBlockInterpolationCompressor {
    public:


        SZBlockInterpolationCompressor(const Config &conf,
                                       Predictor predictor, Quantizer quantizer, Encoder encoder, Lossless lossless,
                                       int interpolator, int direction, int interpo_level) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(conf.absErrorBound)),
                predictor(predictor), quantizer(quantizer), encoder(encoder), lossless(lossless),
                block_size(conf.block_size), stride(conf.stride),
                num_elements(conf.num),
//                interpolators({"linear", "cubic", "cubic2", "akima", "pchip"}),
                interpolators({"linear", "cubic"}),
                interpolator_op(interpolator), direction_op(direction), interp_level(interpo_level) {

            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());

            static_assert(std::is_base_of<concepts::PredictorInterface<T, N>, Predictor>::value,
                          "must implement the predictor interface");
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");

        }


        T *decompress(uchar *compressed_data, const size_t length, bool pre_de_lossless = false) {
            size_t remaining_length = length;
            uchar *lossless_decompressed;
            uchar const *compressed_data_pos;
            if (pre_de_lossless) {
                compressed_data_pos = compressed_data;
            } else {
                lossless_decompressed = lossless.decompress(compressed_data, remaining_length);
                compressed_data_pos = lossless_decompressed;
                double eb;
                read(eb, compressed_data_pos, remaining_length);
            }

            read(global_dimensions.data(), N, compressed_data_pos, remaining_length);
            num_elements = 1;
            for (const auto &d: global_dimensions) {
                num_elements *= d;
                std::cout << d << " ";
            }
            std::cout << std::endl;
            read(block_size, compressed_data_pos, remaining_length);
            stride = block_size;
            predictor.load(compressed_data_pos, remaining_length);
            quantizer.load(compressed_data_pos, remaining_length);
            encoder.load(compressed_data_pos, remaining_length);
            quant_inds = encoder.decode(compressed_data_pos, num_elements);

            encoder.postprocess_decode();

            std::vector<int> block_selection;
            size_t selection_size = *reinterpret_cast<const size_t *>(compressed_data_pos);
            compressed_data_pos += sizeof(size_t);
            remaining_length -= sizeof(size_t);
            HuffmanEncoder<int> selection_encoder;
            selection_encoder.load(compressed_data_pos, remaining_length);
            block_selection = selection_encoder.decode(compressed_data_pos, selection_size);
            selection_encoder.postprocess_decode();

            if (!pre_de_lossless) {
                lossless.postdecompress_data(lossless_decompressed);
            }

            auto dec_data = std::make_unique<T[]>(num_elements);
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.get(),
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions),
                                                                                         block_size,
                                                                                         0);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.get(),
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions),
                                                                                         1, 0);

            predictor.predecompress_data(inter_block_range->begin());
            fallback_predictor.predecompress_data(inter_block_range->begin());
            quantizer.predecompress_data();

//            debug.resize(num_elements, 0);

            auto inter_begin = inter_block_range->begin();
            auto inter_end = inter_block_range->end();
            std::array<size_t, N> intra_block_dims;
            size_t block_idx = 0;

            for (auto block = inter_begin; block != inter_end; block++) {
//                auto block_global_idx = block.get_global_index();
                auto interp_end_idx = block.get_global_index();
                uint max_interp_level = 1;
                for (int i = 0; i < N; i++) {
                    interp_end_idx[i] += intra_block_dims[i] - 1;
                    if (max_interp_level < ceil(log2(intra_block_dims[i]))) {
                        max_interp_level = (uint) ceil(log2(intra_block_dims[i]));
                    }
                }
                intra_block_range->update_block_range(block, block_size);

                if (block_selection[block_idx++] == 0) {
                    concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                    if (!predictor.predecompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        *element = quantizer.recover(predictor_withfallback->predict(element), quant_inds[quant_index++]);
                    }
                } else {
                    if (!interp_level) {
                        auto first = intra_block_range->begin();
                        *first = quantizer.recover(fallback_predictor.predict(first), quant_inds[quant_index++]);
                    } else {
                        max_interp_level = interp_level;
                        size_t stride_sz = 1U << max_interp_level;
                        auto interp_stationary_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                                dec_data.get(), std::begin(global_dimensions), std::end(global_dimensions), stride_sz, 0);
                        std::array<size_t, N> interp_stationary_dims;
                        for (int i = 0; i < N; i++) {
                            interp_stationary_dims[i] = ceil(1.0 * intra_block_dims[i] / stride_sz);
                        }
                        interp_stationary_range->update_block_range(block, interp_stationary_dims);
//                        interp_stationary_range->set_dimensions(interp_stationary_dims.begin(), interp_stationary_dims.end());
//                        interp_stationary_range->set_offsets(block.get_offset());
//                        interp_stationary_range->set_starting_position(block.get_local_index());
                        concepts::PredictorInterface<T, N> *interp_stationary_predictor = &predictor;
                        if (!predictor.predecompress_block(intra_block_range)) {
                            interp_stationary_predictor = &fallback_predictor;
                        }
                        for (auto element = interp_stationary_range->begin();
                             element != interp_stationary_range->end(); ++element) {
                            *element = quantizer.recover(interp_stationary_predictor->predict(element),
                                                         quant_inds[quant_index++]);
//                            debug[element.get_offset()]++;
                        }
//                            std::cout << "quan: " << quant_inds.size() << std::endl;
                    }
                    for (uint level = max_interp_level; level > 0 && level <= max_interp_level; level--) {
                        size_t stride_ip = 1U << (level - 1);
                        block_interpolation(dec_data.get(), block.get_global_index(), interp_end_idx, PB_recover,
                                            interpolators[interpolator_op], direction_op, stride_ip);
                    }
                }
            }

            assert(quant_index == num_elements);
            predictor.postdecompress_data(inter_block_range->begin());

//            fallback_predictor.postdecompress_data(inter_block_range->begin());
            quantizer.postdecompress_data();

            return dec_data.release();
        }


        // compress given the error bound
        uchar *compress(T *data, size_t &compressed_size) {
            quant_inds.clear();
            uchar *compressed_data = new uchar[2 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            std::vector<int> block_selection;
            data2 = std::vector<T>(data, data + num_elements);

//            debug.resize(num_elements, 0);

            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions),
                                                                                         block_size, 0);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1, 0);
            std::array<size_t, N> intra_block_dims;
            predictor.precompress_data(inter_block_range->begin());
            quantizer.precompress_data();
            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    auto block_global_idx = block.get_global_index();
                    auto interp_end_idx = block.get_global_index();
                    uint max_interp_level = 1;
                    for (int i = 0; i < N; i++) {
                        intra_block_dims[i] = (block_global_idx[i] + block_size > global_dimensions[i]) ?
                                              global_dimensions[i] - block_global_idx[i] : block_size;
                        interp_end_idx[i] += intra_block_dims[i] - 1;
                        if (max_interp_level < ceil(log2(intra_block_dims[i]))) {
                            max_interp_level = (uint) ceil(log2(intra_block_dims[i]));
                        }
                    }

                    intra_block_range->update_block_range(block, block_size);
                    concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                    if (!predictor.precompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    double sz_predict_error = 1;
//                    for (auto element = intra_block_range->begin(); element != intra_block_range->end(); ++element) {
//                        sz_predict_error += predictor_withfallback->estimate_error(element);
//                    }

                    double interp_predict_error = 0;

//                    for (uint level = max_interp_level; level > 0 && level <= max_interp_level; level--) {
//                        uint stride_ip = 1U << (level - 1);
//                        interp_predict_error += block_interpolation(data, block.get_global_index(), interp_end_idx, PB_predict,
//                                                                    interpolators[interpolator_op], direction_op, stride_ip);
//                    }

                    if (sz_predict_error < interp_predict_error) {
                        predictor_withfallback->precompress_block_commit();
                        for (auto element = intra_block_range->begin(); element != intra_block_range->end(); ++element) {
                            quant_inds.push_back(quantizer.quantize_and_overwrite(
                                    *element, predictor_withfallback->predict(element)));
                        }
                        block_selection.push_back(0);
                    } else {
                        if (!interp_level) {
                            auto first = intra_block_range->begin();
                            quant_inds.push_back(quantizer.quantize_and_overwrite(*first, fallback_predictor.predict(first)));
//                            debug[first.get_offset()]++;
                        } else {
                            max_interp_level = std::min(interp_level, max_interp_level);
                            size_t stride_sz = 1U << max_interp_level;
                            auto interp_stationary_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                                    data, std::begin(global_dimensions), std::end(global_dimensions), stride_sz, 0);
                            std::array<size_t, N> interp_stationary_dims;
                            for (int i = 0; i < N; i++) {
                                interp_stationary_dims[i] = ceil(1.0 * intra_block_dims[i] / stride_sz);
//                                if (block_selection.size() == 11) {
//                                   std::cout << "Dim " << interp_stationary_dims[i] << std::endl;
//                                }
                            }
                            interp_stationary_range->update_block_range(block, interp_stationary_dims);
//                            interp_stationary_range->set_dimensions(interp_stationary_dims.begin(), interp_stationary_dims.end());
//                            interp_stationary_range->set_offsets(block.get_offset());
//                            interp_stationary_range->set_starting_position(block.get_local_index());
                            concepts::PredictorInterface<T, N> *interp_stationary_predictor = &predictor;
                            if (!predictor.precompress_block(intra_block_range)) {
                                interp_stationary_predictor = &fallback_predictor;
                            }
                            interp_stationary_predictor->precompress_block_commit();
                            for (auto element = interp_stationary_range->begin();
                                 element != interp_stationary_range->end(); ++element) {
                                auto quan = quantizer.quantize_and_overwrite(
                                        *element, interp_stationary_predictor->predict(element));
                                quant_inds.push_back(quan);
//                                debug[element.get_offset()]++;
//                                if (block_selection.size() == 11) {
//                                    std::cout << "Ele " << element.get_offset() << std::endl;
//                                }
                            }
                        }
                        for (uint level = max_interp_level; level > 0 && level <= max_interp_level; level--) {
                            uint stride_ip = 1U << (level - 1);
                            block_interpolation(data, block.get_global_index(), interp_end_idx, PB_predict_overwrite,
                                                interpolators[interpolator_op], direction_op, stride_ip);
                        }
                        block_selection.push_back(1);
                    }
                }
            }
//            for (size_t i = 0; i < num_elements; i++) {
//                if (debug[i] != 1) {
//                    std::cout << i << "," << debug[i] << std::endl;
//                }
//            }

            assert(quant_inds.size() == num_elements);

            {
                std::vector<size_t> cnt(2, 0);
                size_t cnt_total = 0;
                for (auto &sel: block_selection) {
                    cnt[sel]++;
                    cnt_total++;
                }
//                printf("Interp Percentage = %.3f\nSZ Percentage = %.3f\n", 1.0 * cnt[1] / cnt_total, 1.0 * cnt[0] / cnt_total);
            }
            quantizer.postcompress_data();
//            predictor.print();

            write(quantizer.get_eb(), compressed_data_pos);
            write(global_dimensions.data(), N, compressed_data_pos);
            write(block_size, compressed_data_pos);
            predictor.save(compressed_data_pos);
            quantizer.save(compressed_data_pos);

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

            *reinterpret_cast<size_t *>(compressed_data_pos) = (size_t) block_selection.size();
            compressed_data_pos += sizeof(size_t);
            HuffmanEncoder<int> selection_encoder;
            selection_encoder.preprocess_encode(block_selection, 4);
            selection_encoder.save(compressed_data_pos);
            selection_encoder.encode(block_selection, compressed_data_pos);
            selection_encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);

            return lossless_data;
        }

    private:

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover
        };

        size_t offset2(std::array<size_t, N> idx) {
            size_t offset = idx[0];
            for (int i = 1; i < N; i++) {
                offset = offset * global_dimensions[i] + idx[i];
            }
            return offset;
        }

        template<class... Args>
        size_t offset(Args &&... args) {
            std::array<size_t, N> idx{static_cast<size_t>(std::forward<Args>(args))...};
            return offset2(idx);
        }

        inline void quantize(T &d, T pred) {
            quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
        }

        inline void recover(T &d, T pred) {
            d = quantizer.recover(pred, quant_inds[quant_index++]);
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
                        quantize(*d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            quantize(*d, *(d - stride));
                        } else {
                            quantize(*d, interp_linear1(*(d - stride3x), *(d - stride)));
                        }
                    }
                } else {
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
                }
            } else {
                if (pb == PB_predict_overwrite) {

                    T *d = data + begin + stride;
                    quantize(*d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    for (size_t i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        quantize(*d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    if (n % 2 == 0) {
                        d = data + begin + (n - 3) * stride;
                        quantize(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                        d += 2 * stride;
                        quantize(*d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    } else {
                        d = data + begin + (n - 2) * stride;
                        quantize(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                    }
                } else {
                    T *d = data + begin + stride;
                    recover(*d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    for (size_t i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        recover(*d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    if (n % 2 == 0) {
                        d = data + begin + (n - 3) * stride;
                        recover(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                        d += 2 * stride;
                        recover(*d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    } else {
                        d = data + begin + (n - 2) * stride;
                        recover(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                    }

                }
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 1, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride_ip = 1) {
            return block_interpolation_1d(data, offset2(begin), offset2(end), 1, interp_func, pb);
        }

        template<uint NN = N>
        typename std::enable_if<NN == 2, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride_ip = 1) {
            double predict_error = 0;
            size_t stride_ip2 = stride_ip * 2;
            if (direction == 0) {
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip2) {
                    predict_error += block_interpolation_1d(data, offset(begin[0], j), offset(end[0], j),
                                                            stride_ip * global_dimensions[1],
                                                            interp_func, pb);
                }
                for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                    predict_error += block_interpolation_1d(data, offset(i, begin[1]), offset(i, end[1]), stride_ip, interp_func,
                                                            pb);
                }
            } else {
                for (size_t i = begin[0]; i <= end[0]; i += stride_ip2) {
                    predict_error += block_interpolation_1d(data, offset(i, begin[1]), offset(i, end[1]), stride_ip, interp_func,
                                                            pb);
                }
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip) {
                    predict_error += block_interpolation_1d(data, offset(begin[0], j), offset(end[0], j),
                                                            stride_ip * global_dimensions[1],
                                                            interp_func, pb);
                }
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 4, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride_ip = 1) {
            double predict_error = 0;
            size_t stride_ip2 = stride_ip * 2;
            for (size_t j = begin[1]; j <= end[1]; j += stride_ip2) {
                for (size_t k = begin[2]; k <= end[2]; k += stride_ip2) {
                    for (size_t t = (begin[3] ? begin[3] + stride_ip2 : 0); t <= end[3]; t += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(begin[0], j, k, t), offset(end[0], j, k, t),
                                                                stride_ip * global_dimensions[1] * global_dimensions[2] *
                                                                global_dimensions[3],
                                                                interp_func, pb);
                    }
                }
            }
            for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                for (size_t k = begin[2]; k <= end[2]; k += stride_ip2) {
                    for (size_t t = (begin[3] ? begin[3] + stride_ip2 : 0); t <= end[3]; t += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(i, begin[1], k, t), offset(i, end[1], k, t),
                                                                stride_ip * global_dimensions[2] * global_dimensions[3],
                                                                interp_func, pb);
                    }
                }
            }
            for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip) {
                    for (size_t t = (begin[3] ? begin[3] + stride_ip2 : 0); t <= end[3]; t += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(i, j, begin[2], t), offset(i, j, end[2], t),
                                                                stride_ip * global_dimensions[3],
                                                                interp_func, pb);
                    }
                }
            }
            for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip) {
                    for (size_t k = begin[2]; k <= end[2]; k += stride_ip) {
                        predict_error += block_interpolation_1d(data, offset(i, j, k, begin[3]), offset(i, j, k, end[3]),
                                                                stride_ip, interp_func, pb);
                    }
                }
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 3, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride_ip = 1) {
            double predict_error = 0;
            size_t stride_ip2 = stride_ip * 2;

            if (direction == 0 || direction == 1) {
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip2) {
                    for (size_t k = begin[2]; k <= end[2]; k += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(begin[0], j, k), offset(end[0], j, k),
                                                                stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                interp_func,
                                                                pb);
                    }
                }
                if (direction == 0) {
                    for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                        for (size_t k = begin[2]; k <= end[2]; k += stride_ip2) {
                            predict_error += block_interpolation_1d(data, offset(i, begin[1], k), offset(i, end[1], k),
                                                                    stride_ip * global_dimensions[2],
                                                                    interp_func, pb);
                        }
                    }
                    for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                        for (size_t j = begin[1]; j <= end[1]; j += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(i, j, begin[2]), offset(i, j, end[2]), stride_ip,
                                                                    interp_func,
                                                                    pb);
                        }
                    }
                } else {
                    for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                        for (size_t j = begin[1]; j <= end[1]; j += stride_ip2) {
                            predict_error += block_interpolation_1d(data, offset(i, j, begin[2]), offset(i, j, end[2]), stride_ip,
                                                                    interp_func,
                                                                    pb);
                        }
                    }
                    for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                        for (size_t k = begin[2]; k <= end[2]; k += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(i, begin[1], k), offset(i, end[1], k),
                                                                    stride_ip * global_dimensions[2],
                                                                    interp_func, pb);
                        }
                    }
                }

            } else if (direction == 2 || direction == 3) {
                for (size_t k = begin[0]; k <= end[0]; k += stride_ip2) {
                    for (size_t j = begin[2]; j <= end[2]; j += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(k, begin[1], j), offset(k, end[1], j),
                                                                stride_ip * global_dimensions[2],
                                                                interp_func,
                                                                pb);
                    }
                }
                if (direction == 2) {
                    for (size_t i = begin[1]; i <= end[1]; i += stride_ip) {
                        for (size_t j = begin[2]; j <= end[2]; j += stride_ip2) {
                            predict_error += block_interpolation_1d(data, offset(begin[0], i, j), offset(end[0], i, j),
                                                                    stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                    interp_func,
                                                                    pb);
                        }
                    }
                    for (size_t k = begin[0]; k <= end[0]; k += stride_ip) {
                        for (size_t i = begin[1]; i <= end[1]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(k, i, begin[2]), offset(k, i, end[2]),
                                                                    stride_ip,
                                                                    interp_func, pb);
                        }
                    }
                } else {
                    for (size_t k = begin[0]; k <= end[0]; k += stride_ip2) {
                        for (size_t i = begin[1]; i <= end[1]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(k, i, begin[2]), offset(k, i, end[2]),
                                                                    stride_ip,
                                                                    interp_func, pb);
                        }
                    }
                    for (size_t i = begin[1]; i <= end[1]; i += stride_ip) {
                        for (size_t j = begin[2]; j <= end[2]; j += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(begin[0], i, j), offset(end[0], i, j),
                                                                    stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                    interp_func,
                                                                    pb);
                        }
                    }
                }

            } else if (direction == 4 || direction == 5) {

                for (size_t j = begin[0]; j <= end[0]; j += stride_ip2) {
                    for (size_t k = begin[1]; k <= end[1]; k += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(j, k, begin[2]), offset(j, k, end[2]),
                                                                stride_ip, interp_func, pb);
                    }
                }
                if (direction == 4) {
                    for (size_t k = begin[1]; k <= end[1]; k += stride_ip2) {
                        for (size_t i = begin[2]; i <= end[2]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(begin[0], k, i), offset(end[0], k, i),
                                                                    stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                    interp_func, pb);
                        }
                    }
                    for (size_t j = begin[0]; j <= end[0]; j += stride_ip) {
                        for (size_t i = begin[2]; i <= end[2]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(j, begin[1], i), offset(j, end[1], i),
                                                                    stride_ip * global_dimensions[2], interp_func,
                                                                    pb);
                        }
                    }
                } else {
                    for (size_t j = begin[0]; j <= end[0]; j += stride_ip2) {
                        for (size_t i = begin[2]; i <= end[2]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(j, begin[1], i), offset(j, end[1], i),
                                                                    stride_ip * global_dimensions[2], interp_func,
                                                                    pb);
                        }
                    }
                    for (size_t k = begin[1]; k <= end[1]; k += stride_ip) {
                        for (size_t i = begin[2]; i <= end[2]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(begin[0], k, i), offset(end[0], k, i),
                                                                    stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                    interp_func, pb);
                        }
                    }
                }
            }
            return predict_error;
        }

        int interpolator_op;
        int direction_op;
        double sz_eb_ratio = 1;
        uint interp_level;
        std::vector<std::string> interpolators;
        std::vector<int> quant_inds;
        std::vector<int> block_select_quant;
        size_t quant_index = 0; // for decompress
//        std::vector<int> debug;
        std::vector<T> data2;
        Predictor predictor;
        LorenzoPredictor<T, N, 1> fallback_predictor;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
    };

    template<class T, uint N, class Predictor, class Quantizer, class Encoder, class Lossless>
    SZBlockInterpolationCompressor<T, N, Predictor, Quantizer, Encoder, Lossless>
    make_sz_block_interpolation_compressor(const Config &conf, Predictor predictor, Quantizer quantizer,
                                           Encoder encoder, Lossless lossless,
                                           int interp_op,
                                           int direction_op,
                                           uint interp_level) {
        return SZBlockInterpolationCompressor<T, N, Predictor, Quantizer, Encoder, Lossless>(conf, predictor, quantizer,
                                                                                             encoder, lossless,
                                                                                             interp_op, direction_op,
                                                                                             interp_level
        );
    }

};


#endif

