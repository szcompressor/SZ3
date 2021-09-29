#ifndef _SZ_SZ_PROG_INTERPOLATION_MULTILEVEL_QUANTIZATION_HPP
#define _SZ_SZ_PROG_INTERPOLATION_MULTILEVEL_QUANTIZATION_HPP

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
#include <utils/ByteUtil.hpp>
#include <utils/ska_hash/unordered_map.hpp>
#include <utils/Verification.hpp>

namespace SZ {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZProgressiveMQuant {
    public:


        SZProgressiveMQuant(Quantizer quantizer, Encoder encoder, Lossless lossless,
                            const std::array<size_t, N> dims,
                            int interpolator,
                            int direction_id_,
                            size_t interp_dim_limit,
                            int level_independent_,
                            size_t block_size_,
                            int level_fill_) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                global_dimensions(dims),
                interpolators({"linear", "cubic"}),
                interpolator_id(interpolator),
                interp_dim_limit(interp_dim_limit),
                level_independent(level_independent_),
                block_size(block_size_),
                level_fill(level_fill_) {
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
//            static_assert(std::is_base_of<concepts::EncoderInterface<>, Encoder>::value,
//                          "must implement the encoder interface");
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
                global_begin[i] = 0;
                global_end[i] = global_dimensions[i] - 1;
            }

            dim_offsets[N - 1] = 1;
            for (int i = N - 2; i >= 0; i--) {
                dim_offsets[i] = dim_offsets[i + 1] * global_dimensions[i + 1];
            }

            std::array<int, N> base_direction;
            std::array<int, N - 1> stride_multiplication;
            for (int i = 0; i < N - 1; i++) {
                base_direction[i] = i;
                stride_multiplication[i] = 2;
            }
            base_direction[N - 1] = N - 1;

            int direction_id = 0;
            do {
                if (direction_id_ == direction_id) {
                    for (int i = 0; i < N - 1; i++) {
                        auto direction = base_direction;
                        std::rotate(direction.begin() + i, direction.begin() + i + 1, direction.end());
                        directions[i] = std::pair{direction, stride_multiplication};
                        stride_multiplication[i] = 1;
                    }
                    directions[N - 1] = std::pair{base_direction, stride_multiplication};
                    break;
                }
                direction_id++;
            } while (std::next_permutation(base_direction.begin(), base_direction.end()));

            for (int i = 0; i < N; i++) {
                printf("direction %d is ", i);
                for (int j = 0; j < N; j++) {
                    printf("%d ", directions[i].first[j]);
                }
                printf("\n");
            }
        }


        T *decompress(uchar const *lossless_data, const std::vector<size_t> &lossless_size, T *data) {
            Timer timer(true);

            int lossless_id = 0;
            size_t remaining_length = lossless_size[lossless_id];
            retrieved_size += remaining_length;
            uchar const *data_header = lossless_data;

            read(global_dimensions.data(), N, data_header, remaining_length);
            num_elements = std::accumulate(global_dimensions.begin(), global_dimensions.end(), (size_t) 1, std::multiplies<>());

            read(interp_dim_limit, data_header, remaining_length);
            lossless_data += lossless_size[lossless_id];
            lossless_id++;

            T *dec_data = new T[num_elements];
            size_t quant_inds_count = 0;

//            check_dec = true;
            for (uint level = levels; level > level_fill && level > 1; level--) {
                timer.start();
                size_t stride = 1U << (level - 1);

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
                    for (const auto &direction: directions) {
                        block_interpolation(dec_data, block.get_global_index(), end_idx, PB_recover,
                                            interpolators[interpolator_id], direction, stride, true);
                    }
                }
                quantizer.postdecompress_data();
                std::cout << "Level = " << level << " , quant size = " << quant_inds.size() << ", Time = " << timer.stop() << std::endl;
                quant_inds_count += quant_index;
            }


            for (uint level = level_fill; level > 1; level--) {
                timer.start();
                size_t stride = 1U << (level - 1);
                for (const auto &direction: directions) {
                    block_interpolation(dec_data, global_begin, global_end, PB_fill,
                                        interpolators[interpolator_id], direction, stride, true);
                }
                std::cout << "Fill Level = " << level << " " << std::endl;
            }


            if (level_independent == 0) {
                for (const auto &direction: directions) {
                    block_interpolation(dec_data, global_begin, global_end, PB_fill,
                                        interpolators[interpolator_id], direction, 1, true);
                }
                std::cout << "Fill Level = " << 1 << " " << std::endl;
            }


            std::vector<int> quant_sign;
            size_t quant_size;
            bool bitplanetogether = false;
            for (int b = 0, bitstart = 0; b < std::min<int>(bitplane.size(), level_independent); b++) {
                timer.start();

                size_t length = lossless_size[lossless_id];
                retrieved_size += length;

                uchar *compressed_data = lossless.decompress(lossless_data, length);
                uchar const *compressed_data_pos = compressed_data;

                if (b == 0) {
                    quantizer.load(compressed_data_pos, length);
                    read(quant_size, compressed_data_pos, length);
                    quant_inds.clear();
                    quant_inds.resize(quant_size, 0);
                    dec_delta.resize(num_elements, 0);

                    quant_sign = decode_int_2bits(compressed_data_pos, length);
                }

                std::vector<int> quant_ind_truncated;
                if (bitplane[b] == 2) {
                    quant_ind_truncated = decode_int_2bits(compressed_data_pos, length);
                } else {
                    encoder.load(compressed_data_pos, length);
                    quant_ind_truncated = encoder.decode(compressed_data_pos, quant_size);
                    encoder.postprocess_decode();
                }

                lossless.postdecompress_data(compressed_data);
                lossless_data += lossless_size[lossless_id];
                lossless_id++;

                quantizer.reset_unpred_index();
                quant_index = 0;

                printf("\n************Bitplane = %d *****************\n", b);
                int bitshift = 32 - bitstart - bitplane[b];
                if (bitplanetogether) {
                    for (size_t i = 0; i < quant_size; i++) {
                        if (quant_sign[i] == 0) { //unpredictable data
                            quant_inds[i] = -quantizer.get_radius();
                        } else if (quant_sign[i] == 1) { // pred >= 0
                            quant_inds[i] += (quant_ind_truncated[i] << bitshift);
                        } else { // pred < 0
                            quant_inds[i] += -(quant_ind_truncated[i] << bitshift);
                        }
                    }
                } else {

                    for (size_t i = 0; i < quant_size; i++) {
                        if (quant_sign[i] == 0) { //unpredictable data
                            quant_inds[i] = -quantizer.get_radius();
                        } else if (quant_sign[i] == 1) { // pred >= 0
                            quant_inds[i] = (quant_ind_truncated[i] << bitshift);
                        } else { // pred < 0
                            quant_inds[i] = -(quant_ind_truncated[i] << bitshift);
                        }
                    }

                    for (const auto &direction: directions) {
                        block_interpolation(dec_data, global_begin, global_end, (b == 0 ? PB_recover : PB_recover_delta),
                                            interpolators[interpolator_id], direction, 1, true);
                    }
                    SZ::verify<float>(data, dec_data, num_elements);

                }
                bitstart += bitplane[b];
                printf("retrieved = %.3f%% %lu\n", retrieved_size * 100.0 / (num_elements * sizeof(float)), retrieved_size);
                printf("Decompression time = %.3f\n", timer.stop());

            }
            if (bitplanetogether) {
                for (const auto &direction: directions) {
                    block_interpolation(dec_data, global_begin, global_end, PB_recover,
                                        interpolators[interpolator_id], direction, 1, true);
                }
            }
            printf("\nretrieved = %.3f%% %lu\n", retrieved_size * 100.0 / (num_elements * sizeof(float)), retrieved_size);
            return dec_data;
        }

        // compress given the error bound
        uchar *compress(T *data, std::vector<size_t> &lossless_size) {
            Timer timer(true);

            quant_inds.reserve(num_elements);
//            quant_inds.resize(num_elements);
            size_t interp_compressed_size = 0;
//            debug.resize(num_elements, 0);
//            preds.resize(num_elements, 0);
            size_t quant_inds_total = 0;

            T eb = quantizer.get_eb();
            std::cout << "Absolute error bound = " << eb << std::endl;
//            quantizer.set_eb(eb * eb_ratio);

            uchar *lossless_data = new uchar[size_t((num_elements < 1000000 ? 4 : 1.2) * num_elements) * sizeof(T)];
            uchar *lossless_data_pos = lossless_data;

            write(global_dimensions.data(), N, lossless_data_pos);
            write(interp_dim_limit, lossless_data_pos);
            lossless_size.push_back(lossless_data_pos - lossless_data);

            for (uint level = levels; level > 1; level--) {
                timer.start();

                quantizer.set_eb((level >= 3) ? eb * eb_ratio : eb);
                uint stride = 1U << (level - 1);
//                std::cout << "Level = " << level << ", stride = " << stride << std::endl;

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
                    for (const auto &direction: directions) {
                        block_interpolation(data, block.get_global_index(), end_idx, PB_predict_overwrite,
                                            interpolators[interpolator_id], direction, stride, true);
                    }
                }

                auto quant_size = quant_inds.size();
                quant_inds_total += quant_size;
                encode_lossless(lossless_data_pos, lossless_size);

                printf("level = %d , quant size = %lu , time=%.3f\n", level, quant_size, timer.stop());

            }

            timer.start();
            quantizer.set_eb(eb);
            for (const auto &direction: directions) {
                block_interpolation(data, global_begin, global_end, PB_predict_overwrite,
                                    interpolators[interpolator_id], direction, 1, true);
            }
            printf("level = %d , quant size = %lu , prediction time=%.3f\n", 1, quant_inds.size(), timer.stop());
            quant_inds_total += quant_inds.size();

//            for (int i = 0; i < 10; i++) {
//                printf("%d ", quant_inds[i]);
//            }
//            printf("\n");

            Timer timer2(true);
            timer.start();
            int radius = quantizer.get_radius();
            size_t bsize = bitplane.size();
            size_t qsize = quant_inds.size();
            std::vector<int> quants(qsize);
            uchar *compressed_data = new uchar[size_t((quant_inds.size() < 1000000 ? 10 : 1.2) * quant_inds.size()) * sizeof(T)];

            for (int b = 0, bitstart = 0; b < bsize; b++) {
                timer.start();
                uint32_t bitshifts = 32 - bitstart - bitplane[b];
                uint32_t bitmasks = (1 << bitplane[b]) - 1;
                bitstart += bitplane[b];
                uchar *compressed_data_pos = compressed_data;

                if (b == 0) {
                    quantizer.save(compressed_data_pos);
                    quantizer.clear();

                    write((size_t) qsize, compressed_data_pos);

                    std::vector<int> quant_sign(qsize);
                    for (size_t i = 0; i < qsize; i++) {
                        if (quant_inds[i] == -radius) {
                            quant_sign[i] = 0;
                            quant_inds[i] = 0;
                        } else if (quant_inds[i] < 0) {
                            quant_inds[i] = -quant_inds[i];
                            quant_sign[i] = 2;
                        } else {
                            quant_sign[i] = 1;
                        }
                    }
                    encode_int_2bits(quant_sign, compressed_data_pos);

                }
                for (size_t i = 0; i < qsize; i++) {
                    quants[i] = ((uint32_t) quant_inds[i]) >> bitshifts & bitmasks;
                }

                if (bitplane[b] == 2) {
                    encode_int_2bits(quants, compressed_data_pos);
                } else {
                    encoder.preprocess_encode(quants, 0);
                    encoder.save(compressed_data_pos);
                    encoder.encode(quants, compressed_data_pos);
                    encoder.postprocess_encode();
                }

                size_t size = 0;
                uchar *lossless_result = lossless.compress(
                        compressed_data, compressed_data_pos - compressed_data, size);
                memcpy(lossless_data_pos, lossless_result, size);
                lossless_data_pos += size;
                lossless_size.push_back(size);
                delete[]lossless_result;
                printf("Bitplane = %d, time = %.3f\n", b, timer.stop());
            }
            delete[]compressed_data;

//            quant_inds.clear();
            std::cout << "total element = " << num_elements << ", quantization element = " << quant_inds_total << std::endl;
            assert(quant_inds_total >= num_elements);

            timer2.stop("multilevel_quantization");
            return lossless_data;
        }

    private:

        inline void encode_int_2bits(const std::vector<int> &data, uchar *&c) {

            size_t intLen = data.size();
            size_t byteLen = intLen * 2 / 8 + (intLen % 4 == 0 ? 0 : 1);

            write(intLen, c);
            write(byteLen, c);

            size_t b, i = 0;
            int mod4 = intLen % 4;
            for (b = 0; b < (mod4 == 0 ? byteLen : byteLen - 1); b++, i += 4) {
                c[b] = (data[i] << 6) | (data[i + 1] << 4) | (data[i + 2] << 2) | data[i + 3];
            }
            if (mod4 > 0) {
                if (mod4 == 1) {
                    c[b] = (data[i] << 6);
                } else if (mod4 == 2) {
                    c[b] = (data[i] << 6) | (data[i + 1] << 4);
                } else if (mod4 == 3) {
                    c[b] = (data[i] << 6) | (data[i + 1] << 4) | (data[i + 2] << 2);
                }
            }
            c += byteLen;
        }


        std::vector<int> decode_int_2bits(const uchar *&c, size_t &remaining_length) {
            size_t byteLen, intLen;
            read(intLen, c, remaining_length);
            read(byteLen, c, remaining_length);
            std::vector<int> ints(intLen);
            size_t i = 0, b = 0;

            int mod4 = intLen % 4;
            for (; b < (mod4 == 0 ? byteLen : byteLen - 1); b++, i += 4) {
                ints[i] = (c[b] & 0xC0) >> 6;
                ints[i + 1] = (c[b] & 0x30) >> 4;
                ints[i + 2] = (c[b] & 0x0C) >> 2;
                ints[i + 3] = c[b] & 0x03;
            }
            if (mod4 > 0) {
                if (mod4 >= 1) {
                    ints[i] = (c[b] & 0xC0) >> 6;
                }
                if (mod4 >= 2) {
                    ints[i + 1] = (c[b] & 0x30) >> 4;
                }
                if (mod4 >= 3) {
                    ints[i + 2] = (c[b] & 0x0C) >> 2;
                }
            }
            c += byteLen;
            remaining_length -= byteLen;
            return ints;
        }

        void global_index(size_t offset) {
            std::array<size_t, N> global_idx{0};
            for (int i = N - 1; i >= 0; i--) {
                global_idx[i] = offset % global_dimensions[i];
                offset /= global_dimensions[i];
            }
            for (const auto &id: global_idx) {
                printf("%lu ", id);
            }
        }

        void lossless_decode(uchar const *&lossless_data_pos, const std::vector<size_t> &lossless_size, int lossless_id) {

            size_t remaining_length = lossless_size[lossless_id];
            retrieved_size += remaining_length;

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
            uchar *compressed_data = new uchar[size_t((quant_inds.size() < 1000000 ? 10 : 1.2) * quant_inds.size()) * sizeof(T)];
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
            PB_predict_overwrite, PB_predict, PB_recover, PB_fill, PB_recover_delta
        };

        inline void quantize1(size_t idx, T &d, T pred) {
            auto quant = quantizer.quantize_and_overwrite(d, pred);
//            quant_inds.push_back(quant);
            quant_inds[idx] = quant;
        }

        inline void quantize(size_t idx, T &d, T pred) {
            quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
        }

        inline void recover(size_t idx, T &d, T pred) {
            d = quantizer.recover(pred, quant_inds[quant_index++]);
        };


        inline void recover_delta(size_t idx, T &d, T pred) {
            quantizer.recover_delta(d, dec_delta[idx], pred, quant_inds[quant_index++]);
//            if (check_dec) {
//                if (fabs(d - ori_data[idx]) > 1.1e-3) {
//                    printf("ERROR %lu %lu %.5f %.5f ", quant_index-1, idx, ori_data[idx], d);
//                    print_global_index(idx);
//                    printf("\n");
//                }
//            }
        };

        inline void fill(size_t idx, T &d, T pred) {
            d = pred;
        }

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
                } else if (pb == PB_recover_delta) {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        recover_delta(d - data, *d, interp_linear(dec_delta[d - stride - data], dec_delta[d + stride - data]));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            recover_delta(d - data, *d, dec_delta[d - stride - data]);
                        } else {
                            recover_delta(d - data, *d, interp_linear1(dec_delta[d - stride3x - data], dec_delta[d - stride - data]));
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

                    d = data + begin + stride;
                    quantize(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        quantize(d - data, *d,
                                 interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }

                    d = data + begin + i * stride;
                    quantize(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        quantize(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    }

                } else if (pb == PB_recover) {
                    T *d;
                    size_t i;

                    d = data + begin + stride;
                    recover(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        recover(d - data, *d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }

                    d = data + begin + i * stride;
                    recover(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        recover(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    }
                } else if (pb == PB_recover_delta) {
                    T *d;
                    size_t i;

                    d = data + begin + stride;
                    recover_delta(d - data, *d,
                                  interp_quad_1(dec_delta[d - stride - data], dec_delta[d + stride - data], dec_delta[d + stride3x - data]));

                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        recover_delta(d - data, *d,
                                      interp_cubic(dec_delta[d - stride3x - data], dec_delta[d - stride - data], dec_delta[d + stride - data],
                                                   dec_delta[d + stride3x - data]));
                    }

                    d = data + begin + i * stride;
                    recover_delta(d - data, *d,
                                  interp_quad_2(dec_delta[d - stride3x - data], dec_delta[d - stride - data], dec_delta[d + stride - data]));

                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        recover_delta(d - data, *d,
                                      interp_quad_3(dec_delta[d - stride5x - data], dec_delta[d - stride3x - data], dec_delta[d - stride - data]));
                    }
                } else {
                    T *d;
                    size_t i;

                    d = data + begin + stride;
                    fill(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        fill(d - data, *d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }

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


        void block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, PredictorBehavior pb,
                                 const std::string &interp_func, const std::pair<std::array<int, N>, std::array<int, N - 1>> direction, uint stride,
                                 bool overlap) {

            auto dims = direction.first;
            auto s = direction.second;

            if (N == 1) {
                block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb);
            } else if (N == 2) {
                for (size_t i = begin[dims[0]] + ((overlap && begin[dims[0]]) ? stride * s[0] : 0); i <= end[dims[0]]; i += stride * s[0]) {
                    size_t begin_offset = i * dim_offsets[dims[0]] + begin[dims[1]] * dim_offsets[dims[1]];
                    block_interpolation_1d(data, begin_offset, begin_offset + (end[dims[1]] - begin[dims[1]]) * dim_offsets[dims[1]],
                                           stride * dim_offsets[dims[1]], interp_func, pb);
                }
            } else if (N == 3) {
                for (size_t i = begin[dims[0]] + ((overlap && begin[dims[0]]) ? stride * s[0] : 0); i <= end[dims[0]]; i += stride * s[0]) {
                    for (size_t j = begin[dims[1]] + ((overlap && begin[dims[1]]) ? stride * s[1] : 0); j <= end[dims[1]]; j += stride * s[1]) {
                        size_t begin_offset = i * dim_offsets[dims[0]] + j * dim_offsets[dims[1]] + begin[dims[2]] * dim_offsets[dims[2]];
                        block_interpolation_1d(data, begin_offset, begin_offset + (end[dims[2]] - begin[dims[2]]) * dim_offsets[dims[2]],
                                               stride * dim_offsets[dims[2]], interp_func, pb);
                    }
                }
            } else {
                for (size_t i = begin[dims[0]] + ((overlap && begin[dims[0]]) ? stride * s[0] : 0); i <= end[dims[0]]; i += stride * s[0]) {
                    for (size_t j = begin[dims[1]] + ((overlap && begin[dims[1]]) ? stride * s[1] : 0); j <= end[dims[1]]; j += stride * s[1]) {
                        for (size_t k = begin[dims[2]] + ((overlap && begin[dims[2]]) ? stride * s[2] : 0); k <= end[dims[2]]; k += stride * s[2]) {
                            size_t begin_offset = i * dim_offsets[dims[0]] + j * dim_offsets[dims[1]] + k * dim_offsets[dims[2]] +
                                                  begin[dims[3]] * dim_offsets[dims[3]];
                            block_interpolation_1d(data, begin_offset, begin_offset + (end[dims[3]] - begin[dims[3]]) * dim_offsets[dims[3]],
                                                   stride * dim_offsets[dims[3]], interp_func, pb);
                        }
                    }
                }
            }
        }

        int levels = -1;
        int level_independent = -1;
        int level_fill = 0;
        int interpolator_id;
        size_t interp_dim_limit, block_size;
        double eb_ratio = 0.5;
        std::vector<std::string> interpolators;
        std::vector<int> quant_inds;
        size_t quant_index = 0; // for decompress
        size_t num_elements;
        std::array<size_t, N> global_dimensions, global_begin, global_end;
        std::array<size_t, N> dim_offsets;
        std::array<std::pair<std::array<int, N>, std::array<int, N - 1>>, N> directions;
//        int direction_sequence_id;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;

//        std::vector<int> bitplane = {8, 8, 8, 2, 2, 2, 1, 1};
        std::vector<int> bitplane = {24, 4, 2, 2};
        std::vector<T> dec_delta;
        size_t retrieved_size = 0;

        //debug only
//        bool check_dec = false;
        double max_error;
//        float eb;
//        std::vector<T> ori_data;
//        std::vector<int> debug;
//        std::vector<size_t> debug_idx_list;
//        std::vector<int> debug;
//        std::vector<T> preds;

    };


};


#endif

