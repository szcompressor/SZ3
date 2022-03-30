#ifndef _SZ_QOI_INTERPOLATION_COMPRESSSOR_HPP
#define _SZ_QOI_INTERPOLATION_COMPRESSSOR_HPP

#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include <cstring>
#include <cmath>

namespace SZ {
    template<class T, uint N, class Quantizer, class Quantizer_EB, class Encoder, class Lossless>
    class SZQoIInterpolationCompressor {
    public:


        SZQoIInterpolationCompressor(Quantizer quantizer, Quantizer_EB quantizer_eb, std::shared_ptr<concepts::QoIInterface<T, N>> qoi, Encoder encoder, Lossless lossless) :
                quantizer(quantizer), quantizer_eb(quantizer_eb), qoi(qoi), encoder(encoder), lossless(lossless) {

            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) {
            Timer timer(true);
            timer.start();
            size_t remaining_length = cmpSize;
            uchar *buffer = lossless.decompress(cmpData, remaining_length);
            uchar const *buffer_pos = buffer;

            // timer.stop("Lossless");
            // timer.start();

            read(global_dimensions.data(), N, buffer_pos, remaining_length);
            read(blocksize, buffer_pos, remaining_length);
            read(interpolator_id, buffer_pos, remaining_length);
            read(direction_sequence_id, buffer_pos, remaining_length);

            init();

            // std::cout << "load quantizers, offset = " << buffer_pos - buffer << std::endl;
            quantizer_eb.load(buffer_pos, remaining_length);
            quantizer.load(buffer_pos, remaining_length);
            // std::cout << "load encoder, offset = " << buffer_pos - buffer << std::endl;
            encoder.load(buffer_pos, remaining_length);
            quant_inds = encoder.decode(buffer_pos, num_elements);

            encoder.postprocess_decode();
            // std::cout << "after encoder, offset = " << buffer_pos - buffer << std::endl;

            lossless.postdecompress_data(buffer);
            // timer.stop("Decoding");
            // double eb = qoi->get_global_eb();

            // std::cout << "start with first data\n";
            // *decData = quantizer.recover(0, quant_inds[quant_index++]);
            // timer.start();
            quant_index = 0;
            recover_data(0, *decData, 0);

            // std::cout << "start decompressing data\n";
            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
                // if (level >= 3) {
                //     quantizer.set_eb(eb * eb_ratio);
                // } else {
                //     quantizer.set_eb(eb);
                // }
                size_t stride = 1U << (level - 1);
                auto inter_block_range = std::make_shared<
                        SZ::multi_dimensional_range<T, N>>(decData,
                                                           std::begin(global_dimensions), std::end(global_dimensions),
                                                           stride * blocksize, 0);
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
                    block_interpolation(decData, block.get_global_index(), end_idx, PB_recover,
                                        interpolators[interpolator_id], direction_sequence_id, stride);

                }
            }
            quantizer.postdecompress_data();
            // timer.stop("Interpolation Decompress");

            return decData;
        }

        // compress given the error bound
        uchar *compress(const Config &conf, T *data, size_t &compressed_size) {
            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
            blocksize = conf.interpBlockSize;
            interpolator_id = conf.interpAlgo;
            direction_sequence_id = conf.interpDirection;

            // Timer timer(true);
            // timer.start();
            init();

            // quant_inds.reserve(num_elements);
            // 0 ~ num_elements - 1: quant_inds_eb
            // num_elements ~ 2*num_elements - 1: quant_inds_data
            quant_inds = std::vector<int>(num_elements * 2);
            quant_index = 0;
            size_t interp_compressed_size = 0;

            double eb = qoi->get_global_eb();
//            printf("Absolute error bound = %.5f\n", eb);

            // quant_inds.push_back(quantizer.quantize_and_overwrite(*data, 0));
            quantize_data(0, data, 0);

            // Timer timer;
            // timer.start();

            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
                if (level >= 3) {
                    // quantizer.set_eb(eb * eb_ratio);
                    qoi->set_global_eb(eb * eb_ratio);
                } else {
                    // quantizer.set_eb(eb);
                    qoi->set_global_eb(eb);
                }
                uint stride = 1U << (level - 1);
//                std::cout << "Level = " << level << ", stride = " << stride << std::endl;

                auto inter_block_range = std::make_shared<
                        SZ::multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                           std::end(global_dimensions),
                                                           blocksize * stride, 0);

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
//            std::cout << "Number of data point = " << num_elements << std::endl;
//            std::cout << "quantization element = " << quant_inds.size() << std::endl;
            // assert(quant_inds.size() == num_elements);
            // timer.stop("Prediction & Quantization");
            // timer.start();
//            writefile("pred.dat", preds.data(), num_elements);
//            writefile("quant.dat", quant_inds.data(), num_elements);
            size_t bufferSize = 1.5 * (quant_inds.size() * sizeof(T) + quantizer.size_est());
            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;

            write(global_dimensions.data(), N, buffer_pos);
            write(blocksize, buffer_pos);
            write(interpolator_id, buffer_pos);
            write(direction_sequence_id, buffer_pos);

            // std::cout << "save quantizers, offset = " << buffer_pos - buffer << std::endl;
            quantizer_eb.save(buffer_pos);
            quantizer_eb.postcompress_data();
            quantizer.save(buffer_pos);
            quantizer.postcompress_data();


            // timer.start();
            // std::cout << "save encoder, offset = " << buffer_pos - buffer << std::endl;
            encoder.preprocess_encode(quant_inds, 0);
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();
            // std::cout << "after encoder, offset = " << buffer_pos - buffer << std::endl;
            // timer.stop("Encoding");
            // timer.start();
            assert(buffer_pos - buffer < bufferSize);

            // timer.start();
            uchar *lossless_data = lossless.compress(buffer,
                                                     buffer_pos - buffer,
                                                     compressed_size);
            lossless.postcompress_data(buffer);
            // timer.stop("Lossless");
            // std::cout << "quant_index = " << quant_index << ", num_elements = " << num_elements << std::endl;
            compressed_size += interp_compressed_size;
            return lossless_data;
        }

        void clear() {
            quantizer.clear();
            quantizer_eb.clear();
        }

    private:

        inline void quantize_data(size_t offset, T * data, T pred){
            auto ori_data = *data;
            auto eb = qoi->interpret_eb(data, offset);
            quant_inds[quant_index] = quantizer_eb.quantize_and_overwrite(eb);
            quant_inds[num_elements + quant_index] = quantizer.quantize_and_overwrite(
                    *data, pred, eb);
            if(!qoi->check_compliance(ori_data, *data)){
                // std::cout << "not compliant" << std::endl;
                // save as unpredictable
                eb = 0;
                *data = ori_data;
                quant_inds[quant_index] = quantizer_eb.quantize_and_overwrite(eb);
                if(quant_inds[num_elements + quant_index] != 0){
                    // avoiding push multiple data
                    quant_inds[num_elements + quant_index] = quantizer.quantize_and_overwrite(
                            *data, 0, 0);                    
                }
            }
            // update cumulative tolerance if needed 
            qoi->update_tolerance(ori_data, *data);
            quant_index ++;
        }

        inline void recover_data(size_t idx, T& data, T pred) {
            auto eb = quantizer_eb.recover(quant_inds[quant_index]);
            data = quantizer.recover(pred, quant_inds[num_elements + quant_index], eb);
            quant_index ++;
        };

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover
        };

        void init() {
            assert(blocksize % 2 == 0 && "Interpolation block size should be even numbers");
            num_elements = 1;
            interpolation_level = -1;
            for (int i = 0; i < N; i++) {
                if (interpolation_level < ceil(log2(global_dimensions[i]))) {
                    interpolation_level = (uint) ceil(log2(global_dimensions[i]));
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
                        quantize_data(d - data, d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            quantize_data(d - data, d, *(d - stride));
                        } else {
                            quantize_data(d - data, d, interp_linear1(*(d - stride3x), *(d - stride)));
                        }
                    }
                } else {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        recover_data(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            recover_data(d - data, *d, *(d - stride));
                        } else {
                            recover_data(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                        }
                    }
                }
            } else {
                if (pb == PB_predict_overwrite) {

                    T *d;
                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        quantize_data(d - data, d,
                                 interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    d = data + begin + stride;
                    quantize_data(d - data, d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    d = data + begin + i * stride;
                    quantize_data(d - data, d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        quantize_data(d - data, d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    }

                } else {
                    T *d;

                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        recover_data(d - data, *d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    d = data + begin + stride;

                    recover_data(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    d = data + begin + i * stride;
                    recover_data(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        recover_data(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    }
                }
            }

            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 1, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride = 1) {
            return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb);
        }

        template<uint NN = N>
        typename std::enable_if<NN == 2, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            const std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                        stride * dimension_offsets[dims[0]], interp_func, pb);
            }
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                        stride * dimension_offsets[dims[1]], interp_func, pb);
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 3, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            const std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[0]] - begin[dims[0]]) *
                                                            dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[0]], interp_func, pb);
                }
            }
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[1]] - begin[dims[1]]) *
                                                            dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb);
                }
            }
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          begin[dims[2]] * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[2]] - begin[dims[2]]) *
                                                            dimension_offsets[dims[2]],
                                                            stride * dimension_offsets[dims[2]], interp_func, pb);
                }
            }
            return predict_error;
        }


        template<uint NN = N>
        typename std::enable_if<NN == 4, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            max_error = 0;
            const std::array<int, N> dims = dimension_sequences[direction];
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

        int interpolation_level = -1;
        uint blocksize;
        int interpolator_id;
        double eb_ratio = 0.5;
        std::vector<std::string> interpolators = {"linear", "cubic"};
        std::vector<int> quant_inds;
        size_t quant_index = 0; // for decompress
        double max_error;
        Quantizer_EB quantizer_eb;
        Quantizer quantizer;
        std::shared_ptr<concepts::QoIInterface<T, N>> qoi;
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

