#ifndef PQL_COMPRESSOR_HPP
#define PQL_COMPRESSOR_HPP

#include "compressor/Compressor.hpp"
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
#include "utils/Timer.hpp"

namespace SZ {
    template<class T, uint N, class Predictor, class Quantizer, class Lossless>
    class PQLCompressor : public concepts::CompressorInterface<T> {
    public:


        PQLCompressor(const Config<T, N> &conf,
                      Predictor predictor, Quantizer quantizer, Lossless lossless) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(conf.eb)),
                predictor(predictor), quantizer(quantizer), lossless(lossless),
                block_size(conf.block_size), stride(conf.stride),
                global_dimensions(conf.dims), num_elements(conf.num) {
            static_assert(std::is_base_of<concepts::PredictorInterface<T, N>, Predictor>::value,
                          "must implement the predictor interface");
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        uchar *compress(T *data, size_t &compressed_size) {

            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions),
                                                                                         stride, 0);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1,
                                                                                         0);
            std::array<size_t, N> intra_block_dims;
            int *quant_inds = new int[num_elements * 2];
            predictor.precompress_data(inter_block_range->begin());
            quantizer.precompress_data();
            size_t quant_count = 0;
            Timer timer;
            timer.start();
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    // std::cout << *block << " " << lp.predict(block) << std::endl;
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1 &&
                                               global_dimensions[i] - cur_index * stride < block_size) ?
                                              global_dimensions[i] - cur_index * stride : block_size;
                    }

                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());
                    concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                    if (!predictor.precompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    predictor_withfallback->precompress_block_commit();
//                    quantizer.precompress_block();
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        quant_inds[quant_count++] = quantizer.quantize_and_overwrite(
                                *element, predictor_withfallback->predict(element));
                    }
                }
            }

            timer.stop("Predition & Quantization");

            predictor.postcompress_data(inter_block_range->begin());
            quantizer.postcompress_data();


            uchar *compressed_data = reinterpret_cast<uchar *>(quant_inds);
            uchar *compressed_data_pos = compressed_data + num_elements * sizeof(int);
            write(block_size, compressed_data_pos);
            predictor.save(compressed_data_pos);
            quantizer.save(compressed_data_pos);

            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
                                                     compressed_size);
//            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }

        T *decompress(uchar const *lossless_compressed_data, const size_t length) {
            size_t remaining_length = length;
            auto compressed_data = lossless.decompress(lossless_compressed_data, remaining_length);
            uchar const *compressed_data_pos = compressed_data;
            compressed_data_pos += num_elements * sizeof(int);
            read(block_size, compressed_data_pos, remaining_length);
            stride = block_size;
            predictor.load(compressed_data_pos, remaining_length);
            quantizer.load(compressed_data_pos, remaining_length);

            int const *quant_inds_pos = reinterpret_cast<const int *> (compressed_data);
            std::array<size_t, N> intra_block_dims;
            auto dec_data = new T[num_elements];
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions),
                                                                                         block_size,
                                                                                         0);

            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1,
                                                                                         0);

            predictor.predecompress_data(inter_block_range->begin());
            quantizer.predecompress_data();

            std::cout << "start decompression" << std::endl;
            Timer timer;
            timer.start();
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; block++) {
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size
                                                                      : block_size;
                    }
                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());

                    concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                    if (!predictor.predecompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        *element = quantizer.recover(predictor_withfallback->predict(element), *(quant_inds_pos++));
                    }
                }
            }
            timer.stop("Predition & Quantization");
            lossless.postdecompress_data(compressed_data);
            predictor.postdecompress_data(inter_block_range->begin());
            quantizer.postdecompress_data();
            return dec_data;
        }


    private:
        Predictor predictor;
        LorenzoPredictor<T, N, 1> fallback_predictor;
        Quantizer quantizer;
        Lossless lossless;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
    };

    template<class T, uint N, class Predictor, class Quantizer, class Lossless>
    PQLCompressor<T, N, Predictor, Quantizer, Lossless>
    make_pql_compressor(const Config<T, N> &conf, Predictor predictor, Quantizer quantizer,
                        Lossless lossless) {
        return PQLCompressor<T, N, Predictor, Quantizer, Lossless>(conf, predictor, quantizer, lossless);
    }
}
#endif
