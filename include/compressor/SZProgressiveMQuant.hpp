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
#include "utils/ByteUtil.hpp"
#include "utils/ska_hash/unordered_map.hpp"
#include "utils/Verification.hpp"
#include "def.hpp"
#include <cstring>
#include <cmath>

namespace SZ {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZProgressiveMQuant {
    public:


        SZProgressiveMQuant(Quantizer quantizer, Encoder encoder, Lossless lossless,
                            const std::array<size_t, N> dims,
                            int interpolator,
                            int direction_id_,
                            size_t interp_dim_limit,
                            int level_progressive_,
                            size_t block_size_,
                            int level_fill_) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                global_dimensions(dims),
                interpolators({"linear", "cubic"}),
                interpolator_id(interpolator),
                interp_dim_limit(interp_dim_limit),
                level_progressive(level_progressive_),
                block_size(block_size_)
//                level_fill(level_fill_)
        {
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
            set_directions_and_stride(direction_id_);

        }


        T *decompress(uchar const *lossless_data, const std::vector<size_t> &lossless_size, T *data) {
            Timer timer(true);

            size_t buffer_len = lossless_size[0];
            retrieved_size += buffer_len;
            uchar const *buffer = lossless_data;
            read(global_dimensions.data(), N, buffer, buffer_len);
            num_elements = std::accumulate(global_dimensions.begin(), global_dimensions.end(), (size_t) 1, std::multiplies<>());
            read(interp_dim_limit, buffer, buffer_len);
            l2_diff.resize(level_progressive * N * bitgroup.size() , 0);
            read(l2_diff.data(), l2_diff.size(), buffer, buffer_len);

            //load unpredictable data
            buffer = lossless_data;
            for (int i = 0; i < lossless_size.size() - 1; i++) {
                buffer += lossless_size[i];
            }
            buffer_len = lossless_size[lossless_size.size() - 1];
            retrieved_size += buffer_len;
            buffer = lossless.decompress(buffer, buffer_len);
//            uchar const *buffer_pos = buffer;
            quantizer.load(buffer, buffer_len);

            printf("retrieved = %.3f%% %lu\n", retrieved_size * 100.0 / (num_elements * sizeof(float)), retrieved_size);

            int lossless_id = 1;
            lossless_data += lossless_size[0];

            T *dec_data = new T[num_elements];

            for (uint level = levels; level > level_progressive; level--) {
                timer.start();
                size_t stride = 1U << (level - 1);

                lossless_decode(lossless_data, lossless_size, lossless_id++);

                if (level == levels) {
                    *dec_data = quantizer.recover(0, 0, quant_inds[quant_cnt++]);
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
                        block_interpolation(dec_data, dec_data, block.get_global_index(), end_idx, &SZProgressiveMQuant::recover,
                                            interpolators[interpolator_id], direction, stride, true);
                    }
                }
                quantizer.postdecompress_data();
                std::cout << "Level = " << level << " , quant size = " << quant_inds.size() << ", Time = " << timer.stop() << std::endl;
            }

            for (uint level = level_progressive; level > 0; level--) {
                for (const auto &direction: directions) {
                    block_interpolation(dec_data, dec_data, global_begin, global_end, &SZProgressiveMQuant::recover_no_quant,
                                        interpolators[interpolator_id], direction, 1U << (level - 1), true);
                }
            }

            printf("retrieved = %.3f%% %lu\n", retrieved_size * 100.0 / (num_elements * sizeof(float)), retrieved_size);

            if (level_progressive > 0) {
                {
                    double psnr, nrmse, max_err, range;
                    verify(data, dec_data, num_elements, psnr, nrmse, max_err, range);
                }

                std::vector<ska::unordered_map<std::string, double>> result_stat;

                int lsize = N * level_progressive, bsize = bitgroup.size();
                std::vector<int> bsum(lsize), bdelta(lsize);
                bdelta[0] = 1;
                std::vector<uchar const *> data_lb(lsize * bsize);
                std::vector<size_t> size_lb(lsize * bsize);
                for (int l = 0; l < lsize; l++) {
                    for (int b = bsize - 1; b >= 0; b--) {
                        data_lb[l * bsize + b] = lossless_data;
                        size_lb[l * bsize + b] = lossless_size[lossless_id];
                        lossless_data += lossless_size[lossless_id];
                        lossless_id++;
                    }
                }
                dec_delta.clear();
                dec_delta.resize(num_elements, 0);

//                double targetl2 = 60;
                while (true) {
//                    for (uint level = level_progressive; level > 0; level--) {
//                        for (int direct = 0; direct < N; direct++) {
//                            int lid = (level_progressive - level) * N + direct;
//                            l2_diff[lid * bsize + b]);
//                        }
//                    }


                    bool changed = false;
                    ska::unordered_map<std::string, double> result;
                    for (uint level = level_progressive; level > 0; level--) {
                        for (int direct = 0; direct < N; direct++) {
                            int lid = (level_progressive - level) * N + direct;
                            if (bdelta[lid] > 0 && bsum[lid] < bsize) {
                                printf("\n-----------------------\n");
                                printf("Level = %d , direction = %d , lid = %d, bg = %d\n", level, direct, lid, bsum[lid]);
                                changed = true;
                                result["level"] = level;
                                result["direct"] = direct;
                                quant_inds.clear();
                                quant_cnt = 0;
                                int bg_end = std::min(bsize, bsum[lid] + bdelta[lid]);
                                for (int b = bsum[lid]; b < bg_end; b++) {
                                    printf("reduce l2 = %.10G\n", l2_diff[lid * bsize + b]);
                                    uchar const *bg_data = data_lb[lid * bsize + b];
                                    size_t bg_len = size_lb[lid * bsize + b];
                                    lossless_decode_bitgroup(b, bg_data, bg_len);
                                }
                                if (bsum[lid] == 0) {
                                    block_interpolation(dec_data, dec_data, global_begin, global_end, &SZProgressiveMQuant::recover,
                                                        interpolators[interpolator_id], directions[direct], 1U << (level - 1), true);
                                } else {
                                    block_interpolation(dec_data, dec_delta.data(), global_begin, global_end,
                                                        &SZProgressiveMQuant::recover_set_delta,
                                                        interpolators[interpolator_id], directions[direct], 1U << (level - 1), true);
                                }
                                bsum[lid] = bg_end;
                                double psnr, nrmse, max_err, range;
                                verify(data, dec_data, num_elements, psnr, nrmse, max_err, range);

                            } else {
                                if (changed) {
                                    block_interpolation(dec_data, dec_delta.data(), global_begin, global_end,
                                                        &SZProgressiveMQuant::recover_set_delta_no_quant,
                                                        interpolators[interpolator_id], directions[direct], 1U << (level - 1), true);
                                }
                            }
                        }
                    }

//                    bool last = true;
//                    for (int b = 0; b < lsize - 1; b++) {
//                        if (bsum[b] < bsize) {
//                            last = false;
//                            break;
//                        }
//                    }
//                    if (last) {
//                        bdelta[lsize - 1] = 1;
//                    } else {
                    for (int b = 0; b < lsize; b++) {
                        if (bdelta[b]) {
                            bdelta[b] = 0;
                            bdelta[(b + 1) % (lsize)] = (b == lsize - 2 ? 1 : 1);
                            break;
                        }
                    }
//                    }

                    if (!changed) {
                        continue;
                    }

                    dec_delta.clear();
                    dec_delta.resize(num_elements, 0);

                    float retrieved_rate = retrieved_size * 100.0 / (num_elements * sizeof(float));
                    printf("\nretrieved = %.3f%% %lu\n", retrieved_rate, retrieved_size);
                    double psnr, nrmse, max_err, range;
                    verify(data, dec_data, num_elements, psnr, nrmse, max_err, range);
                    result["retrieved_rate"] = retrieved_rate;
                    result["psnr"] = psnr;
                    result["max_err"] = max_err;
                    result["max_rel_err"] = max_err / range;
                    result_stat.push_back(result);

                    bool done = true;
                    for (const auto &b: bsum) {
                        printf("%d ", b);
                        if (b < bsize) {
                            done = false;
                        }
                    }
                    printf("\n");
                    if (done) {
                        break;
                    }
                }


                for (int i = 0; i < result_stat.size(); i++) {
                    printf("Level = %d, direction = %d, Rate= %.3f%% , PSNR= %.2f , ABS = %.3G , REL = %.3G\n",
                           (int) result_stat[i]["level"], (int) result_stat[i]["direct"],
                           result_stat[i]["retrieved_rate"], result_stat[i]["psnr"],
                           result_stat[i]["max_err"], result_stat[i]["max_rel_err"]);
                }
            }

            return dec_data;
        }

        // compress given the error bound
        uchar *compress(T *data, std::vector<size_t> &lossless_size) {
            Timer timer(true);

            quant_inds.reserve(num_elements);
//            quant_inds.resize(num_elements);
            size_t interp_compressed_size = 0;
            size_t quant_inds_total = 0;

            T eb = quantizer.get_eb();
            std::cout << "Absolute error bound = " << eb << std::endl;
//            quantizer.set_eb(eb * eb_ratio);
            l2_diff.resize(level_progressive * N * bitgroup.size(), 0);

            uchar *lossless_data = new uchar[size_t((num_elements < 1000000 ? 4 : 1.2) * num_elements) * sizeof(T)];
            uchar *lossless_data_pos = lossless_data;

            write(global_dimensions.data(), N, lossless_data_pos);
            write(interp_dim_limit, lossless_data_pos);
            uchar *error_mse_pos = lossless_data_pos;
            lossless_data_pos += l2_diff.size() * sizeof(T);
            lossless_size.push_back(lossless_data_pos - lossless_data);

            for (uint level = levels; level > level_progressive; level--) {
                timer.start();

//                quantizer.set_eb((level >= 3) ? eb * eb_ratio : eb);
                uint stride = 1U << (level - 1);

                if (level == levels) {
                    quant_inds.push_back(quantizer.quantize_and_overwrite(0, *data, 0));
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
                        block_interpolation(data, data, block.get_global_index(), end_idx, &SZProgressiveMQuant::quantize,
                                            interpolators[interpolator_id], direction, stride, true);
                    }
                }

                auto quant_size = quant_inds.size();
                quant_inds_total += quant_size;
                auto size = encode_lossless(lossless_data_pos, lossless_size);
                printf("level = %d , quant size = %lu , lossless size = %lu, time=%.3f\n", level, quant_size, size, timer.stop());

            }

            for (uint level = level_progressive; level > 0; level--) {
                timer.start();

//                quantizer.set_eb((level >= 3) ? eb * eb_ratio : eb);
                uint stride = 1U << (level - 1);
                if (level == levels) {
                    quant_inds.push_back(quantizer.quantize_and_overwrite(0, *data, 0));
                }
                for (int d = 0; d < N; d++) {
                    block_interpolation(data, data, global_begin, global_end, &SZProgressiveMQuant::quantize,
                                        interpolators[interpolator_id], directions[d], stride, true);

                    auto quant_size = quant_inds.size();
                    quant_inds_total += quant_size;
                    auto size = encode_lossless_bitplane((level_progressive - level) * N + d, lossless_data_pos, lossless_size, eb);
                    printf("level = %d , direction = %d, quant size = %lu, lossless size = %lu, time=%.3f\n\n",
                           level, d, quant_size, size, timer.stop());

                }
            }

//            quant_inds.clear();
            std::cout << "total element = " << num_elements << ", quantization element = " << quant_inds_total << std::endl;
            assert(quant_inds_total >= num_elements);

            write(l2_diff.data(), l2_diff.size(), error_mse_pos);

            uchar *buffer = new uchar[quantizer.get_unpred_size() * (sizeof(T) + sizeof(size_t)) + 40];
            uchar *buffer_pos = buffer;
            quantizer.save(buffer_pos);
            size_t size = lossless.compress(buffer, buffer_pos - buffer, lossless_data_pos);
            delete[] buffer;
            lossless_data_pos += size;
            lossless_size.push_back(size);

            return lossless_data;
        }

    private:
        typedef void (SZProgressiveMQuant::*PredictionFunc)(size_t, T &, T);

        int levels = -1;
        int level_progressive = -1;
//        int level_fill = 0;
        int interpolator_id;
        size_t interp_dim_limit, block_size;
        double eb_ratio = 0.5;
        std::vector<std::string> interpolators;
        std::vector<int> quant_inds;
        std::vector<T> error;
        std::vector<T> l2_diff;
        size_t quant_cnt = 0; // for decompress
        size_t num_elements;
        std::array<size_t, N> global_dimensions, global_begin, global_end;
        std::array<size_t, N> dim_offsets;
        std::array<std::pair<std::array<int, N>, std::array<int, N - 1>>, N> directions;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;

//        std::vector<int> bitgroup = {8, 8, 8, 2, 2, 2, 1, 1};
//TODO quantization bins in different levels have different distribution.
// a dynamic bitgroup should be used for each level
        std::vector<int> bitgroup = {16, 8, 4, 2, 1, 1};
//        std::vector<int> bitgroup = {30, 2};
//        std::vector<int> bitgroup = {16, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1};
        std::vector<T> dec_delta;
        size_t retrieved_size = 0;

        //debug only
        double max_error;
//        float eb;

        void
        lossless_decode_bitgroup(int bg, uchar const *data_pos, const size_t data_length) {
            Timer timer(true);

            size_t length = data_length;
            retrieved_size += length;

            uchar *compressed_data = lossless.decompress(data_pos, length);
            uchar const *compressed_data_pos = compressed_data;

            size_t quant_size;
            read(quant_size, compressed_data_pos, length);

            std::vector<int> quant_ind_truncated;
            if (bitgroup[bg] == 1) {
                quant_ind_truncated = decode_int_1bit(compressed_data_pos, length);
            } else if (bitgroup[bg] == 2) {
                quant_ind_truncated = decode_int_2bits(compressed_data_pos, length);
            } else {
                encoder.load(compressed_data_pos, length);
                quant_ind_truncated = encoder.decode(compressed_data_pos, quant_size);
                encoder.postprocess_decode();
            }

            lossless.postdecompress_data(compressed_data);


//                printf("\n************Bitplane = %d *****************\n", bg);
            int bitshift = 32;
            for (int bb = 0; bb <= bg; bb++) {
                bitshift -= bitgroup[bb];
            }
            quant_inds.resize(quant_size, 0);
            for (size_t i = 0; i < quant_size; i++) {
                quant_inds[i] += (((uint32_t) quant_ind_truncated[i] << bitshift) ^ 0xaaaaaaaau) - 0xaaaaaaaau;
            }
        }

        size_t encode_lossless_bitplane(int lid, uchar *&lossless_data_pos, std::vector<size_t> &lossless_size, T eb) {
            Timer timer;
            int bsize = bitgroup.size();
            size_t qsize = quant_inds.size();
            std::vector<int> quants(qsize);
            uchar *buffer = new uchar[size_t((quant_inds.size() < 1000000 ? 10 : 1.2)
                                             * quant_inds.size()) * sizeof(T)];
            for (size_t i = 0; i < qsize; i++) {
                quant_inds[i] = ((int32_t) quant_inds[i] + (uint32_t) 0xaaaaaaaau) ^ (uint32_t) 0xaaaaaaaau;
            }

            double l2_error_base = 0;
            for (size_t i = 0; i < qsize; i++) {
                l2_error_base += error[i] * error[i];
            }
            printf("l2 = %.10G \n", l2_error_base);
            size_t total_size = 0;
            int shift = 0;
            for (int b = bsize - 1; b >= 0; b--) {
                timer.start();
                uchar *buffer_pos = buffer;
                write((size_t) qsize, buffer_pos);

                double l2_error = 0;
                for (size_t i = 0; i < qsize; i++) {
                    quants[i] = quant_inds[i] & (((uint64_t) 1 << bitgroup[b]) - 1);
                    quant_inds[i] >>= bitgroup[b];
                    int qu = (((uint32_t) quants[i] << shift) ^ 0xaaaaaaaau) - 0xaaaaaaaau;
                    error[i] += qu * 2.0 * eb;
                    l2_error += error[i] * error[i];
                }
                l2_diff[lid * bsize + b] = l2_error - ((b == bsize - 1) ? l2_error_base : l2_diff[lid * bsize + b + 1]);
                printf("l2 = %.10G , diff = %.10G\n", l2_error, l2_diff[lid * bsize + b]);
                shift += bitgroup[b];

                if (bitgroup[b] == 1) {
                    encode_int_1bit(quants, buffer_pos);
                } else if (bitgroup[b] == 2) {
                    encode_int_2bits(quants, buffer_pos);
                } else {
                    //TODO huffman tree is huge if using large radius on early levels
                    // set different radius for each level
                    encoder.preprocess_encode(quants, 0);
                    encoder.save(buffer_pos);
                    encoder.encode(quants, buffer_pos);
                    encoder.postprocess_encode();
                }

                size_t size = lossless.compress(
                        buffer, buffer_pos - buffer, lossless_data_pos);
//                printf("%d %lu, ", bitgroup[b], size);
                total_size += size;
                lossless_data_pos += size;
                lossless_size.push_back(size);
            }
//            printf("\n");
            delete[]buffer;
            quant_inds.clear();
            error.clear();

            return total_size;
        }

        void lossless_decode(uchar const *&lossless_data_pos, const std::vector<size_t> &lossless_size, int lossless_id) {

            size_t remaining_length = lossless_size[lossless_id];
            retrieved_size += remaining_length;

            uchar *compressed_data = lossless.decompress(lossless_data_pos, remaining_length);
            uchar const *compressed_data_pos = compressed_data;

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
            quant_cnt = 0;

            lossless.postdecompress_data(compressed_data);
            lossless_data_pos += lossless_size[lossless_id];
        }

        size_t encode_lossless(uchar *&lossless_data_pos, std::vector<size_t> &lossless_size) {
            uchar *compressed_data = new uchar[size_t((quant_inds.size() < 1000000 ? 10 : 1.2) * quant_inds.size()) * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;

            write((size_t) quant_inds.size(), compressed_data_pos);
            if (quant_inds.size() < 128) {
                write(quant_inds.data(), quant_inds.size(), compressed_data_pos);
            } else {
                encoder.preprocess_encode(quant_inds, 0);
                encoder.save(compressed_data_pos);
                encoder.encode(quant_inds, compressed_data_pos);
                encoder.postprocess_encode();
            }

            size_t size = lossless.compress(compressed_data, compressed_data_pos - compressed_data,
                                            lossless_data_pos);
            lossless.postcompress_data(compressed_data);

            lossless_data_pos += size;
            lossless_size.push_back(size);

            quant_inds.clear();
            error.clear();

            return size;
        }

        inline void quantize(size_t idx, T &data, T pred) {
            T data0 = data;
            quant_inds.push_back(quantizer.quantize_and_overwrite(idx, data, pred));
            error.push_back(data0 - data);
        }

        inline void recover(size_t idx, T &d, T pred) {
            d = quantizer.recover(idx, pred, quant_inds[quant_cnt++]);
        };

        inline void recover_only_quant(size_t idx, T &d, T pred) {
            d = quantizer.recover(idx, 0, quant_inds[quant_cnt++]);
        };

        inline void recover_no_quant(size_t idx, T &d, T pred) {
            d = quantizer.recover(idx, pred, 0);
        };

        inline void recover_set_delta_no_quant(size_t idx, T &d, T pred) {
            quantizer.recover_and_residual(idx, d, dec_delta[idx], pred);
        };

        inline void recover_set_delta(size_t idx, T &d, T pred) {
            quantizer.recover_and_residual(idx, d, dec_delta[idx], pred, quant_inds[quant_cnt++]);
        };


        inline void fill(size_t idx, T &d, T pred) {
            d = pred;
        }


        double block_interpolation_1d(T *d, T *pd, size_t begin, size_t end, size_t stride,
                                      const std::string &interp_func, PredictionFunc func) {
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }

            size_t c;
            size_t stride3x = stride * 3, stride5x = stride * 5;
            if (interp_func == "linear" || n < 5) {
                size_t i = 1;
                for (i = 1; i + 1 < n; i += 2) {
                    c = begin + i * stride;
                    (this->*func)(c, d[c], interp_linear(pd[c - stride], pd[c + stride]));
                }
                if (n % 2 == 0) {
                    c = begin + (n - 1) * stride;
                    if (n < 4) {
                        (this->*func)(c, d[c], pd[c - stride]);
                    } else {
                        (this->*func)(c, d[c], interp_linear1(pd[c - stride3x], pd[c - stride]));
                    }
                }
            } else {
                size_t i = 1;
                c = begin + i * stride;
                (this->*func)(c, d[c], interp_quad_1(pd[c - stride], pd[c + stride], pd[c + stride3x]));
                for (i = 3; i + 3 < n; i += 2) {
                    c = begin + i * stride;
                    (this->*func)(c, d[c], interp_cubic(pd[c - stride3x], pd[c - stride], pd[c + stride], pd[c + stride3x]));
                }
                c = begin + i * stride;
                (this->*func)(c, d[c], interp_quad_2(pd[c - stride3x], pd[c - stride], pd[c + stride]));
                if (n % 2 == 0) {
                    c = begin + (n - 1) * stride;
                    (this->*func)(c, d[c], interp_quad_3(pd[c - stride5x], pd[c - stride3x], pd[c - stride]));
                }
            }
            return 0;
        }


        void block_interpolation(T *data, T *pred_data, std::array<size_t, N> begin, std::array<size_t, N> end, PredictionFunc func,
                                 const std::string &interp_func, const std::pair<std::array<int, N>, std::array<int, N - 1>> direction,
                                 uint stride, bool overlap) {

            auto dims = direction.first;
            auto s = direction.second;

            if (N == 1) {
                block_interpolation_1d(data, pred_data, begin[0], end[0], stride, interp_func, func);
            } else if (N == 2) {
                for (size_t i = begin[dims[0]] + ((overlap && begin[dims[0]]) ? stride * s[0] : 0); i <= end[dims[0]]; i += stride * s[0]) {
                    size_t begin_offset = i * dim_offsets[dims[0]] + begin[dims[1]] * dim_offsets[dims[1]];
                    block_interpolation_1d(data, pred_data, begin_offset, begin_offset + (end[dims[1]] - begin[dims[1]]) * dim_offsets[dims[1]],
                                           stride * dim_offsets[dims[1]], interp_func, func);
                }
            } else if (N == 3) {
                for (size_t i = begin[dims[0]] + ((overlap && begin[dims[0]]) ? stride * s[0] : 0); i <= end[dims[0]]; i += stride * s[0]) {
                    for (size_t j = begin[dims[1]] + ((overlap && begin[dims[1]]) ? stride * s[1] : 0); j <= end[dims[1]]; j += stride * s[1]) {
                        size_t begin_offset = i * dim_offsets[dims[0]] + j * dim_offsets[dims[1]] + begin[dims[2]] * dim_offsets[dims[2]];
                        block_interpolation_1d(data, pred_data, begin_offset, begin_offset + (end[dims[2]] - begin[dims[2]]) * dim_offsets[dims[2]],
                                               stride * dim_offsets[dims[2]], interp_func, func);
                    }
                }
            } else {
                for (size_t i = begin[dims[0]] + ((overlap && begin[dims[0]]) ? stride * s[0] : 0); i <= end[dims[0]]; i += stride * s[0]) {
                    for (size_t j = begin[dims[1]] + ((overlap && begin[dims[1]]) ? stride * s[1] : 0); j <= end[dims[1]]; j += stride * s[1]) {
                        for (size_t k = begin[dims[2]] + ((overlap && begin[dims[2]]) ? stride * s[2] : 0);
                             k <= end[dims[2]]; k += stride * s[2]) {
                            size_t begin_offset = i * dim_offsets[dims[0]] + j * dim_offsets[dims[1]] + k * dim_offsets[dims[2]] +
                                                  begin[dims[3]] * dim_offsets[dims[3]];
                            block_interpolation_1d(data, pred_data, begin_offset,
                                                   begin_offset + (end[dims[3]] - begin[dims[3]]) * dim_offsets[dims[3]],
                                                   stride * dim_offsets[dims[3]], interp_func, func);
                        }
                    }
                }
            }
        }

        void set_directions_and_stride(int direction_id_) {
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

    };


};


#endif

