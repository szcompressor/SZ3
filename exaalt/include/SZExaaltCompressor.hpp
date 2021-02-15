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


        SZ_Exaalt_Compressor(const Config<T, N> &conf,
                             Quantizer quantizer, Encoder encoder, Lossless lossless) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                global_dimensions(conf.dims), num_elements(conf.num) {
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        void set_level(float level_start_, float level_offset_) {
            this->level_start = level_start_;
            this->level_offset = level_offset_;
        }

        inline int quantize_to_level(T data) {
            return round((data - level_start) / level_offset);
        }

        inline T level(int l) {
            return level_start + l * level_offset;
        }

        // compress given the error bound
        uchar *compress(T *data, size_t &compressed_size) {

            std::vector<int> quant_inds(num_elements);
            std::vector<int> pred_inds(num_elements);
            quantizer.precompress_data();
            struct timespec start, end;

            clock_gettime(CLOCK_REALTIME, &start);

            auto l0 = quantize_to_level(data[0]);
            pred_inds[0] = l0;
            quant_inds[0] = quantizer.quantize_and_overwrite(data[0], level(l0));

            for (size_t i = 1; i < num_elements; i++) {
//                auto d = data[i];
                auto l = quantize_to_level(data[i]);
                pred_inds[i] = l - l0 + pred_offset;
//                if (i % 2 == 0) {
//                    quant_inds.push_back(quantizer.quantize_and_overwrite(data[i], level(1 + quantize_to_level(data[i - 1]))));
//                } else {
//                    quant_inds.push_back(quantizer.quantize_and_overwrite(data[i], data[i - 1]));
//                }
                quant_inds[i] = quantizer.quantize_and_overwrite(data[i], level(l));
                l0 = l;
//                if (i > 10000000 && i < 10000400) {
//                    std::cout << d << " , " << quant_inds[i] - 32768 << std::endl;
//                }
            }

//            writefile("pred_int.dat", pred_inds.data(), num_elements);

//            std::cout << *std::min_element(pred_inds.begin(), pred_inds.end()) << ", "
//                      << *std::max_element(pred_inds.begin(), pred_inds.end()) << std::endl;
//            std::cout << "Prediction Error = " << pred_err / num_elements << std::endl;

            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = "
                      << (double) (end.tv_sec - start.tv_sec) +
                         (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                      << "s" << std::endl;

            quantizer.postcompress_data();
            uchar *compressed_data;
            compressed_data = new uchar[2 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            write(global_dimensions.data(), N, compressed_data_pos);
//            predictor.save(compressed_data_pos);
            quantizer.save(compressed_data_pos);

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

            encoder.preprocess_encode(pred_inds, 100);
            encoder.save(compressed_data_pos);
            encoder.encode(pred_inds, compressed_data_pos);
            encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
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
                l += pred_inds[i] - pred_offset;
                dec_data[i] = quantizer.recover(level(l), quant_inds[i]);
            }

            quantizer.postdecompress_data();
            return dec_data.release();
        }


    private:
        int pred_offset = 26;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
        float level_start = 1;
        float level_offset = 1.8075;
    };


}
#endif

