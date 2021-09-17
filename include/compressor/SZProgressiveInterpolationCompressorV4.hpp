#ifndef _SZ_SZ_PROG_INTERPOLATION_V4_HPP
#define _SZ_SZ_PROG_INTERPOLATION_V4_HPP

#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/FileIterator.hpp"
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
    class SZProgressiveInterpolationCompressorV4 {
    public:


        SZProgressiveInterpolationCompressorV4(Quantizer quantizer, Encoder encoder, Lossless lossless,
                                               const std::array<size_t, N> dims,
                                               int interpolator,
                                               int direction,
                                               size_t interp_dim_limit,
                                               int level_independent_,
                                               size_t interp_block_size,
                                               int level_fill_) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                global_dimensions(dims),
                interpolators({"linear", "cubic"}),
                interpolator_id(interpolator), direction_sequence_id(direction),
                interp_block_limit(interp_dim_limit),
                level_fill(level_fill_), level_independent(level_independent_),
                core_data(0) {
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");

            assert(interp_dim_limit % 2 == 0 &&
                   "Interpolation dimension should be even numbers to avoid extrapolation");
//            level_independent = 1;
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


        void decompress(const char *cmpr_datafile, const std::vector<size_t> &lossless_size, const char *dec_datafile) {
            cmpr_ifs = std::ifstream(cmpr_datafile, std::ios::binary);
            int lossless_id = 0;
            const uchar *lossless_data = nullptr;

            size_t remaining_length = lossless_size[lossless_id];
            std::vector<uchar> header(remaining_length);
            cmpr_ifs.read(reinterpret_cast<char *>(header.data()), remaining_length);

            const uchar *header_pos = header.data();
            read(global_dimensions.data(), N, header_pos, remaining_length);
            read(core_dimensions.data(), N, header_pos, remaining_length);
            read(interp_block_limit, header_pos, remaining_length);
            size_t num_elements = std::accumulate(global_dimensions.begin(), global_dimensions.end(), (size_t) 1, std::multiplies<>());
            size_t core_num_elements = std::accumulate(core_dimensions.begin(), core_dimensions.end(), (size_t) 1, std::multiplies<>());
            size_t block_num_elements = round(pow(block_size + 1, N));
            printf("total num = %lu\ncore num = %lu\nblock num = %lu \n", num_elements, core_num_elements, block_num_elements);
            core_data.resize(core_num_elements);

            read(core_data[0], header_pos, remaining_length);
//            lossless_data += lossless_size[lossless_id++];
            lossless_id++;
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
//                    quantizer.postdecompress_data();
                    quantizer.clear();
                    std::cout << "Level = " << level << " , quant size = " << quant_inds.size()
                              << " , time = " << timer.stop() << std::endl;
                    quant_inds_count += quant_index;
                }
            }

            if (level_independent == 0) {
                writefile(dec_datafile, core_data.data(), num_elements);
            } else {
                std::ofstream fout(dec_datafile, std::ios::binary | std::ios::out);
//                T *dec_data = new T[num_elements];

                std::vector<T> block_data(block_num_elements);
                auto global_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        nullptr, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
                auto global_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        nullptr, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

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

//                        quantizer.postdecompress_data();
                        quantizer.clear();
                        quant_inds_count += quant_index;
                        quant_size[level] += quant_index;
                    }
                    for (uint level = level_fill; level > 0; level--) {
                        size_t size = lossless_size[lossless_id++];
                        cmpr_ifs.seekg(size, std::ios::cur);
//                        lossless_data += lossless_size[lossless_id++];
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
                        std::array<int, N> move_pos{0};
                        move_pos[N - 2] = 1;
                        do {
                            fout.seekp(global_iter.get_offset() * sizeof(T), std::ios::beg);
                            fout.write(reinterpret_cast<const char *>(&block_data[block_iter.get_offset()]), block_dims[N - 1] * sizeof(T));
//                            memcpy(&dec_data[global_iter.get_offset()], &block_data[block_iter.get_offset()], block_dims[N - 1] * sizeof(T));
                            block_iter.move2(move_pos);
                        } while (global_iter.move2(move_pos));
                        copytime2 += timer.stop();
                    }
                }
                for (uint level = level_independent; level > 0; level--) {
                    printf("level = %d , quant size = %lu\n", level, quant_size[level]);
                }
                printf("level %d to 1, comptime = %.3f, copytime = %.3f %.3f\n", level_independent, comptime, copytime1, copytime2);

//                delete[]core_data;
//                return dec_data;
            }
            printf("retrieved = %.3f%% %lu\n", retrieved_size * 100.0 / (num_elements * sizeof(float)), retrieved_size);

        }

        // compress given the error bound
        uchar *compress(const char *ori_datafile, const char *cmpr_datafile, std::vector<size_t> &lossless_size) {
            size_t num_elements = std::accumulate(global_dimensions.begin(), global_dimensions.end(), (size_t) 1, std::multiplies<>());
            size_t core_num_elements = std::accumulate(core_dimensions.begin(), core_dimensions.end(), (size_t) 1, std::multiplies<>());
            size_t block_num_elements = round(pow(block_size + 1, N));
            size_t quant_inds_total = 1;
            T eb = quantizer.get_eb();
            std::cout << "Absolute error bound = " << eb << std::endl;

            Timer timer;
            size_t buffer_size = std::max(core_num_elements, block_num_elements);
            lossless_buffer.reserve(1.2 * sizeof(T) * buffer_size);
            quant_inds.reserve(buffer_size);
            core_data.resize(core_num_elements);
            printf("[Memory] core_data = %.2f MB\n", sizeof(T) * core_data.capacity() / 1000.0 / 1000.0);
            printf("[Memory] quant buffer = %.2f MB\n", sizeof(int) * quant_inds.capacity() / 1000.0 / 1000.0);
            printf("[Memory] lossless buffer = %.2f MB\n", sizeof(uchar) * lossless_buffer.capacity() / 1000.0 / 1000.0);

            if (level_independent == 0) {
                assert(core_num_elements == num_elements);
                SZ::readfile<T>(ori_datafile, num_elements, core_data.data());
            } else {
                auto core_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        core_data.data(), std::begin(core_dimensions), std::end(core_dimensions), 1, 0);
                auto global_range_inter = std::make_shared<SZ::file_multi_dimensional_range<T, N>>(
                        ori_datafile, std::begin(global_dimensions), std::end(global_dimensions), core_stride, 0);

                auto core_iter = core_range->begin();
                auto global_iter = global_range_inter->begin();
                for (; core_iter != core_range->end(); ++core_iter, ++global_iter) {
                    *core_iter = *global_iter;
                }
//                timer.stop("prepare core");
            }

            cmpr_ofs = std::ofstream(cmpr_datafile, std::ios::binary);
            {
                std::vector<uchar> header(N * 2 * sizeof(size_t) + 100);
                uchar *header_pos = header.data();
                write(global_dimensions.data(), N, header_pos);
                write(core_dimensions.data(), N, header_pos);
                write(interp_block_limit, header_pos);
                write(core_data[0], header_pos);
                cmpr_ofs.write(reinterpret_cast<const char *>(header.data()), header_pos - header.data());
                lossless_size.push_back(header_pos - header.data());
            }
            uchar *lossless_pos = nullptr;

            {
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
                    encode_lossless(lossless_pos, lossless_size);
                    printf("level = %d , quant size = %lu , time=%.3f\n", level, quantsize, timer.stop());
                }
            }
//            std::cout << "quantization element = " << quant_inds_total << std::endl;

            if (level_independent > 0) {

                quant_inds.clear();

                double comptime = 0, copytime1 = 0, copytime2 = 0;

                std::vector<T> block_data(block_num_elements);
                auto global_range_inter = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                        nullptr, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
                auto global_range = std::make_shared<SZ::file_multi_dimensional_range<T, N>>(
                        ori_datafile, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

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
//                        std::array<int, N> move_pos{0};
//                        move_pos[N - 2] = 1;
//                        do {
//                            fin.seekg(global_iter.get_offset() * sizeof(T), std::ios::beg);
//                            fin.read(reinterpret_cast<char *>(&block_data[block_iter.get_offset()]), block_dims[N - 1] * sizeof(T));
//                            block_iter.move2(move_pos);
//                        } while (global_iter.move2(move_pos));
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

                        encode_lossless(lossless_pos, lossless_size);

//                        printf("level = %d , #block = %lu, quant_bins=%lu, time=%.3f, encode_time=%.3f lossless_time=%.3f\n",
//                               level, nBlock, tmp, timer.stop(), encode_time, lossless_time);
                    }
                    comptime += timer.stop();

                }
                for (uint level = level_independent; level > 0; level--) {
                    printf("level = %d , quant size = %lu\n", level, quant_size[level]);
                }
                printf("level %d to 1, comptime = %.3f, copytime = %.3f %.3f\n", level_independent, comptime, copytime1, copytime2);
//                delete[] core_data;
            }


            std::cout << "total element = " << num_elements << ", quantization element = " << quant_inds_total << std::endl;
            assert(quant_inds_total >= num_elements);
//            timer.stop("Predition & Quantization");
            return nullptr;

        }

    private:
        std::vector<uchar> lossless_buffer;
        std::vector<int> quant_inds;
        std::vector<T> core_data;

        std::ofstream cmpr_ofs;
        std::ifstream cmpr_ifs;
        int levels = -1, level_independent = -1, level_fill = 0;
        size_t interp_block_limit, block_size, core_blocksize, core_stride;
        int interpolator_id;
        double eb_ratio = 0.5;
        std::vector<std::string> interpolators;
        size_t quant_index = 0; // for decompress
        double max_error;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        std::array<size_t, N> global_dimensions, core_dimensions;
        std::array<size_t, N> dimension_offsets;
        std::vector<std::array<int, N>> dimension_sequences;
        int direction_sequence_id;
        size_t retrieved_size = 0;

        void lossless_decode(uchar const *&lossless_data_pos, const std::vector<size_t> &lossless_size, int lossless_id) {

            size_t remaining_length = lossless_size[lossless_id];
            if (lossless_buffer.size() < remaining_length) {
                lossless_buffer.resize(remaining_length);
            }
            retrieved_size += remaining_length;
            uchar *lossless_buffer_pos = lossless_buffer.data();
            cmpr_ifs.read(reinterpret_cast<char *>(lossless_buffer_pos), remaining_length);

            uchar *compressed_data = lossless.decompress(lossless_buffer_pos, remaining_length);
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
//            lossless_data_pos += lossless_size[lossless_id];
        }

        void encode_lossless(uchar *&lossless_data_pos, std::vector<size_t> &lossless_size) {
            size_t buffer_size = (quant_inds.size() < 10000 ? 4 : 1.2) * quant_inds.size() * sizeof(T);
            if (lossless_buffer.size() < buffer_size) {
                lossless_buffer.resize(buffer_size);
//                printf("resize\n");
            }
            uchar *lossless_buffer_pos = lossless_buffer.data();

            quantizer.save(lossless_buffer_pos);
//            quantizer.postcompress_data();
            quantizer.clear();
            write((size_t) quant_inds.size(), lossless_buffer_pos);
            if (quant_inds.size() < 128) {
                write(quant_inds.data(), quant_inds.size(), lossless_buffer_pos);
            } else {
                encoder.preprocess_encode(quant_inds, 0);

//                write(min, lossless_buffer_pos);
                encoder.save(lossless_buffer_pos);

                encoder.encode(quant_inds, lossless_buffer_pos);

                encoder.postprocess_encode();

            }

            size_t size = 0;
            uchar *compressed = lossless.compress(lossless_buffer.data(), lossless_buffer_pos - lossless_buffer.data(), size);
//            lossless.postcompress_data(lossless_buffer.data());

            cmpr_ofs.write(reinterpret_cast<const char *>(compressed), size);
//            memcpy(lossless_data_pos, compressed, size);
//            lossless_data_pos += size;
            lossless_size.push_back(size);
            delete[]compressed;

            quant_inds.clear();

        }

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover, PB_fill
        };

        inline void quantize1(size_t idx, T &d, T pred) {
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


    };


};


#endif

