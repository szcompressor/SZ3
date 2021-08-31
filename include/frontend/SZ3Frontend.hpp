#ifndef SZ3_FRONTEND
#define SZ3_FRONTEND

#include "Frontend.hpp"
#include "def.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "utils/Iterator.hpp"
#include "utils/Config.hpp"
#include "utils/MemoryUtil.hpp"

namespace SZ {


    template<class T, uint N, class Predictor, class Quantizer>
    class SZ3Frontend : public concepts::FrontendInterface<T, N> {
    public:

        SZ3Frontend(const Config <T, N> &conf, Predictor predictor, Quantizer quantizer) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(conf.eb)),
                predictor(predictor),
                quantizer(quantizer),
                block_size(conf.block_size),
                stride(conf.stride),
                global_dimensions(conf.dims),
                num_elements(conf.num) {
        }

        ~SZ3Frontend() = default;

        std::vector<int> compress(T *data) {
            std::vector<int> quant_inds(num_elements);
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions),
                                                                                         stride, 0);

            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1,
                                                                                         0);
            std::array<size_t, N> intra_block_dims;
            predictor.precompress_data(inter_block_range->begin());
            quantizer.precompress_data();
            size_t quant_count = 0;
            struct timespec start, end;
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

            predictor.postcompress_data(inter_block_range->begin());
            quantizer.postcompress_data();
            return quant_inds;
        }

        T *decompress(std::vector<int> &quant_inds) {

            int const *quant_inds_pos = (int const *) quant_inds.data();
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
            predictor.postdecompress_data(inter_block_range->begin());
            quantizer.postdecompress_data();
            return dec_data;
        }

        void save(uchar *&c) {
            write(global_dimensions.data(), N, c);
            write(block_size, c);

            predictor.save(c);
            quantizer.save(c);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            read(global_dimensions.data(), N, c, remaining_length);
            num_elements = 1;
            for (const auto &d : global_dimensions) {
                num_elements *= d;
                std::cout << d << " ";
            }
            std::cout << std::endl;
            read(block_size, c, remaining_length);
            stride = block_size;
            predictor.load(c, remaining_length);
            quantizer.load(c, remaining_length);
        }

        void print() {
//            predictor.print();
//            quantizer.print();
        }

        void clear() {
            predictor.clear();
            fallback_predictor.clear();
            quantizer.clear();
        }

        int get_radius() const { return quantizer.get_radius(); }

        size_t get_num_elements() const { return num_elements; };

    private:
        Predictor predictor;
        LorenzoPredictor<T, N, 1> fallback_predictor;
        Quantizer quantizer;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;

    };

    template<class T, uint N, class Predictor, class Quantizer>
    SZ3Frontend<T, N, Predictor, Quantizer>
    make_sz3_frontend(const Config <T, N> &conf, Predictor predictor, Quantizer quantizer) {
        return SZ3Frontend<T, N, Predictor, Quantizer>(conf, predictor, quantizer);
    }
}

#endif
