#ifndef _SZ_EXAALT_HPP
#define _SZ_EXAALT_HPP

#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/MemoryUtil.hpp"
#include "utils/Config.hpp"
#include "utils/FileUtil.h"
#include "def.hpp"
#include <cstring>

namespace SZ {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZ_Exaalt_Compressor {
    public:


        SZ_Exaalt_Compressor(const Config<T, N> &conf, Quantizer quantizer, Encoder encoder, Lossless lossless) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                global_dimensions(conf.dims), num_elements(conf.num) {
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        void set_level(float level_start_, float level_offset_, int level_num_) {
            this->level_start = level_start_;
            this->level_offset = level_offset_;
            this->level_num = level_num_;
        }

        inline int quantize_to_level(T data) {
            return round((data - level_start) / level_offset);
        }

        inline T level(int l) {
            return level_start + l * level_offset;
        }

        // compress given the error bound
        uchar *compress2(T *data, size_t &compressed_size) {

            std::vector<int> quant_inds(num_elements);
            std::vector<int> pred_inds(num_elements);
            quantizer.precompress_data();
            struct timespec start, end;

            clock_gettime(CLOCK_REALTIME, &start);


            size_t n1 = 7852;
            size_t n2 = 1037;
            for (size_t i = 0; i < n1; i++) {
                for (size_t j = 0; j < n2; j++) {
                    size_t t = i * n2 + j;
                    T p1 = (t >= 2) ? fabs(data[t - 2] - data[t]) : std::numeric_limits<T>::max();
                    T p2 = (i >= 1) ? fabs(data[(i - 1) * n2 + j] - data[t]) : std::numeric_limits<T>::max();
//                    size_t t2 = i + j * n1;
//                    if (i == 0) {
//                        if (j == 0) {
//                            quant_inds[t2] = quantizer.quantize_and_overwrite(data[t], data[0]);
//                        } else {
//                            quant_inds[t2] = quantizer.quantize_and_overwrite(data[t], data[(n1 - 1) * n2 + j - 1]);
//                        }
//                    } else {
//                        quant_inds[t2] = quantizer.quantize_and_overwrite(data[t], data[(i - 1) * n2 + j]);
//                    }if (j % 128 <= 1) {    if (i >= 1) {        quant_inds[t] = quantizer.quantize_and_overwrite(data[t], data[(i - 1) * n2 + j]);    } else {        if (t >= 2) {            quant_inds[t] = quantizer.quantize_and_overwrite(data[t], data[t - 2]);        } else {            quant_inds[t] = quantizer.quantize_and_overwrite(data[t], 0);        }    }} else {    quant_inds[t] = quantizer.quantize_and_overwrite(data[t], data[t - 2]);}
//                    if (p1 < p2 || (p1 == p2 && i >= 1)) {
//                        quant_inds[t] = quantizer.quantize_and_overwrite(data[t], data[t - 2]);
//                        pred_inds[t] = 0;
//                    } else if (p1 > p2) {
//                        quant_inds[t] = quantizer.quantize_and_overwrite(data[t], data[(i - 1) * n2 + j]);
//                        pred_inds[t] = 1;
//                    } else {
//                        quant_inds[t] = quantizer.quantize_and_overwrite(data[t], 0);
//                        pred_inds[t] = 0;
//                    }
//                    if (t > 100000 && t < 100400) {if (i >= 100 && i < 110) {    std::cout << quant_inds[t] - quantizer.get_radius() << " ";
//                        std::cout << pred_inds[t] << " ";
//                    std::cout << data[i] << " , " << quant_inds[i] - 32768 << std::endl;}
                }
                if (i >= 100 && i < 110) {
                    std::cout << std::endl;
                }
            }

//            writefile("pred_int.dat", pred_inds.data(), num_elements);

//            std::cout << *std::min_element(pred_inds.begin(), pred_inds.end()) << ", "
//                      << *std::max_element(pred_inds.begin(), pred_inds.end()) << std::endl;
//            std::cout << "Prediction Error = " << pred_err / num_elements << std::endl;

            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = " << (double) (end.tv_sec - start.tv_sec) +
                                                               (double) (end.tv_nsec - start.tv_nsec) /
                                                               (double) 1000000000 << "s" << std::endl;

            quantizer.postcompress_data();

            uchar *compressed_data;
            compressed_data = new uchar[2 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            write(global_dimensions.data(), N, compressed_data_pos);
            quantizer.save(compressed_data_pos);

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos
            );
            encoder.postprocess_encode();

//            std::cout << *std::min_element(pred_inds.begin(), pred_inds.end()) << std::endl;
//            std::cout << *std::max_element(pred_inds.begin(), pred_inds.end()) << std::endl;

            encoder.preprocess_encode(pred_inds, level_num * 2 + 1);
            encoder.save(compressed_data_pos);
            encoder.encode(pred_inds, compressed_data_pos);
            encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data, compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }


        // compress given the error bound
        uchar *compress(T *data, size_t &compressed_size) {

//            std::vector<int> quant_inds(num_elements);
//            std::vector<int> pred_inds(num_elements);

            std::vector<int> quant_inds(0);
            std::vector<int> pred_inds(0);
            quant_inds.reserve(num_elements);
            pred_inds.reserve(num_elements);

            quantizer.precompress_data();

            struct timespec start, end;

            clock_gettime(CLOCK_REALTIME, &start);

            auto l0 = quantize_to_level(data[0]);

//            pred_inds[0]=l0;
//            quant_inds[0]=quantizer.quantize_and_overwrite(data[0], level(l0));
            pred_inds.push_back(l0);
            quant_inds.push_back(quantizer.quantize_and_overwrite(data[0], level(l0)));

            std::vector<int> levels(num_elements);
            levels[0] = l0;

            size_t n1 = 5;
            size_t n2 = 2371092;
//            for (size_t i = 1; i < n2; i++) {
            for (size_t i = 1; i < num_elements; i++) {
                auto l = quantize_to_level(data[i]);
                levels[i] = l;
                if (i < 96) {
//                    pred_inds[i] = l - levels[i - 1] + level_num;
                    pred_inds.push_back(l - levels[i - 1] + level_num);
                } else {
//                    pred_inds[i] = l - levels[i - 1] + level_num;
//                    pred_inds[i] = l - levels[i - 96] + level_num;
//                    pred_inds.push_back(l - levels[i - 96] + level_num);
                    pred_inds.push_back(l - levels[i - 1] + level_num);
                }
//                quant_inds[i] = quantizer.quantize_and_overwrite(data[i], level(l));
                quant_inds.push_back(quantizer.quantize_and_overwrite(data[i], level(l)));
                l0 = l;
                if (i > 100000 && i < 100400) {
                    std::cout << pred_inds[i] - level_num << " ";
//                    std::cout << quant_inds[i] - quantizer.get_radius() << " ";
//                    std::cout << data[i] << " , " << quant_inds[i] - 32768 << std::endl;
                }
            }
//            for (size_t i = 0; i < n2; i++) {
//                for (size_t t = 1; t < n1; t++) {
//                    size_t idx = t * n2 + i;
//                    auto l = quantize_to_level(data[idx]);
//                    levels[idx] = l;
////                    pred_inds.push_back(l - levels[(t - 1) * n2 + i] + level_num);
////                    pred_inds.push_back(l - levels[idx - 96] + level_num);
////                    pred_inds.push_back(l - levels[idx - 1] + level_num);
////                    quant_inds.push_back(quantizer.quantize_and_overwrite(data[i], level(l)));
//                    quant_inds.push_back(quantizer.quantize_and_overwrite(data[i], data[(t-1)*n2+i]));
////                    pred_inds[idx] = l - levels[(t - 1) * n2 + i] + level_num;
////                    quant_inds[idx] = quantizer.quantize_and_overwrite(data[i], level(l));
//                }
//            }
//            assert(pred_inds.size() == num_elements);
//            assert(quant_inds.size() == num_elements);
//            writefile("pred_int.dat", pred_inds.data(), num_elements);

//            std::cout << *std::min_element(pred_inds.begin(), pred_inds.end()) << ", "
//                      << *std::max_element(pred_inds.begin(), pred_inds.end()) << std::endl;
//            std::cout << "Prediction Error = " << pred_err / num_elements << std::endl;

            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = " << (double) (end.tv_sec - start.tv_sec) +
                                                               (double) (end.tv_nsec - start.tv_nsec) /
                                                               (double) 1000000000 << "s" << std::endl;

            quantizer.postcompress_data();

            uchar *compressed_data;
            compressed_data = new uchar[2 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            write(global_dimensions.data(), N, compressed_data_pos);
            quantizer.save(compressed_data_pos);
            quantizer.print();

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

//            std::cout << *std::min_element(pred_inds.begin(), pred_inds.end()) << std::endl;
//            std::cout << *std::max_element(pred_inds.begin(), pred_inds.end()) << std::endl;

            encoder.preprocess_encode(pred_inds, level_num * 2 + 1);
            encoder.save(compressed_data_pos);
            encoder.encode(pred_inds, compressed_data_pos);
            encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data, compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }


        uchar *compress3(T *data, size_t &compressed_size) {

            std::vector<int> quant_inds(num_elements);
            quantizer.precompress_data();

            struct timespec start, end;

            clock_gettime(CLOCK_REALTIME, &start);

            quant_inds[0] = quantizer.quantize_and_overwrite(data[0], 0);

            for (size_t i = 1; i < num_elements; i++) {
                if (i < 16 || i > 1023) {
                    quant_inds[i] = quantizer.quantize_and_overwrite(data[i], data[i - 1]);
                } else {
//                    quant_inds[i] = quantizer.quantize_and_overwrite(data[i], data[i - 1]);
                    quant_inds[i] = quantizer.quantize_and_overwrite(data[i], data[i - 16]);
                }
//                if (i > 100000 && i < 100400) {
                std::cout << quant_inds[i] - quantizer.get_radius() << " ";
//                    std::cout << data[i] << " , " << quant_inds[i] - 32768 << std::endl;
//                }
            }
            std::cout << std::endl;

//            writefile("pred_int.dat", pred_inds.data(), num_elements);

//            std::cout << *std::min_element(quant_inds.begin(), quant_inds.end()) << ", "
//                      << *std::max_element(quant_inds.begin(), quant_inds.end()) << std::endl;
//            std::cout << "Prediction Error = " << quant_inds / num_elements << std::endl;

            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = " << (double) (end.tv_sec - start.tv_sec) +
                                                               (double) (end.tv_nsec - start.tv_nsec) /
                                                               (double) 1000000000 << "s" << std::endl;

            quantizer.postcompress_data();

            uchar *compressed_data;
            compressed_data = new uchar[2 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            write(global_dimensions.data(), N, compressed_data_pos);
            quantizer.save(compressed_data_pos);
            quantizer.print();

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

//            std::cout << *std::min_element(pred_inds.begin(), pred_inds.end()) << std::endl;
//            std::cout << *std::max_element(pred_inds.begin(), pred_inds.end()) << std::endl;

            uchar *lossless_data = lossless.compress(compressed_data, compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }


        T *decompress(uchar const *lossless_compressed_data, const size_t length) {
            size_t remaining_length = length;
            auto compressed_data = lossless.decompress(lossless_compressed_data, remaining_length);
            uchar const *compressed_data_pos = compressed_data;
            read(global_dimensions.data(), N, compressed_data_pos, remaining_length);
            num_elements = 1;
            for (const auto &d : global_dimensions) {
                num_elements *= d;
                std::cout << d << " ";
            }
            std::cout << std::endl;
//            predictor.load(compressed_data_pos, remaining_length);
            // std::cout << "load predictor done\n";fflush(stdout);
            quantizer.load(compressed_data_pos, remaining_length);
            // std::cout << "load quantizer done\n";fflush(stdout);
            encoder.load(compressed_data_pos, remaining_length);
            // size_t unpred_data_size = 0;
            // read(unpred_data_size, compressed_data_pos, remaining_length);
            // T const * unpred_data_pos = (T const *) compressed_data_pos;
            // compressed_data_pos += unpred_data_size * sizeof(T);
            auto quant_inds = encoder.decode(compressed_data_pos, num_elements);
            encoder.postprocess_decode();

            encoder.load(compressed_data_pos, remaining_length);
            auto pred_inds = encoder.decode(compressed_data_pos, num_elements);
            encoder.postprocess_decode();

            lossless.postdecompress_data(compressed_data);

            auto dec_data = std::make_unique<T[]>(num_elements);

            quantizer.predecompress_data();

            std::cout << "start decompression" << std::endl;
            auto l = pred_inds[0];
            dec_data[0] = quantizer.recover(level(l), quant_inds[0]);
            for (size_t i = 1; i < num_elements; i++) {
                l += pred_inds[i] - level_num;
                dec_data[i] = quantizer.recover(level(l), quant_inds[i]);
            }

            quantizer.postdecompress_data();
            return dec_data.release();
        }


    private:
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
        float level_start = 1;
        float level_offset = 1.8075;
        int level_num = 26;
    };


}
#endif

