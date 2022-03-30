#ifndef SZ3_QOI_FRONTEND
#define SZ3_QOI_FRONTEND
/**
 * This module contains SZ3 Predictor and Quantizer for QOI preservation
 */

#include "Frontend.hpp"
#include "SZ3/def.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ {


    template<class T, uint N, class Predictor, class Quantizer, class Quantizer_EB>
    class SZQoIFrontend : public concepts::FrontendInterface<T, N> {
    public:
        SZQoIFrontend(const Config &conf, Predictor predictor, Quantizer quantizer, Quantizer_EB quantizer_eb, std::shared_ptr<concepts::QoIInterface<T, N>> qoi) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(conf.absErrorBound)),
                predictor(predictor),
                quantizer(quantizer),
                block_size(conf.blockSize),
                num_elements(conf.num),
                quantizer_eb(quantizer_eb),
                qoi(qoi) {
            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
        }

        std::vector<int> compress(T *data) {
            // 0 ~ num_elements - 1: quant_inds_eb
            // num_elements ~ 2*num_elements - 1: quant_inds_data
            std::vector<int> quant_inds(num_elements * 2);

            auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    data, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);

            auto element_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    data, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

            predictor.precompress_data(block_range->begin());
            quantizer.precompress_data();
            size_t quant_count = 0;
            for (auto block = block_range->begin(); block != block_range->end(); ++block) {

                element_range->update_block_range(block, block_size);
                qoi->precompress_block(element_range);

                concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                if (!predictor.precompress_block(element_range)) {
                    predictor_withfallback = &fallback_predictor;
                }
                predictor_withfallback->precompress_block_commit();

                for (auto element = element_range->begin(); element != element_range->end(); ++element) {
                    auto ori_data = *element;
                    // interpret the error bound for current data based on qoi
                    auto eb = qoi->interpret_eb(element);

                    quant_inds[quant_count] = quantizer_eb.quantize_and_overwrite(eb);
                    quant_inds[num_elements + quant_count] = quantizer.quantize_and_overwrite(
                            *element, predictor_withfallback->predict(element), eb);

                    // if(element.get_offset() == 19733396){
                        // auto pred = predictor_withfallback->predict(element);
                        // std::cout << element.get_offset() << "->" << quant_count << ": eb = " << eb << ", quant_inds = " << quant_inds[quant_count] << std::endl;
                        // std::cout << "eb = " << eb << std::endl;
                        // std::cout << quant_inds[quant_count] << " " << quant_inds[num_elements + quant_count] << std::endl;
                        // std::cout << ori_data << " " << *element << std::endl;
                        // std::cout << "pred = " << predictor_withfallback->predict(element) << std::endl;
                    //     auto temp = qoi.check_compliance(ori_data, *element, true);
                    //     std::cout << "check_compliance = " << temp << std::endl;
                    // }
                    // check whether decompressed data is compliant with qoi tolerance
                    if(!qoi->check_compliance(ori_data, *element)){
                        // std::cout << "exceed in " << element.get_offset() << std::endl;
                        // save as unpredictable
                        eb = 0;
                        *element = ori_data;
                        quant_inds[quant_count] = quantizer_eb.quantize_and_overwrite(eb);
                        if(quant_inds[num_elements + quant_count] != 0){
                            // avoid push multiple elements
                            quant_inds[num_elements + quant_count] = quantizer.quantize_and_overwrite(*element, 0, 0);                            
                        }
                    }
                    quant_count ++;
                    // update cumulative tolerance if needed 
                    qoi->update_tolerance(ori_data, *element);
                }
                qoi->postcompress_block();
            }
            predictor.postcompress_data(block_range->begin());
            quantizer.postcompress_data();
            return quant_inds;
        }

        T *decompress(std::vector<int> &quant_inds, T *dec_data) {

            int const *quant_inds_eb_pos = (int const *) quant_inds.data();
            int const *quant_inds_pos = quant_inds_eb_pos + num_elements;
            std::array<size_t, N> intra_block_dims;
            // auto dec_data = new T[num_elements];
            auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    dec_data, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);

            auto element_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    dec_data, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

            predictor.predecompress_data(block_range->begin());
            quantizer.predecompress_data();

            for (auto block = block_range->begin(); block != block_range->end(); ++block) {

                element_range->update_block_range(block, block_size);

                concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                if (!predictor.predecompress_block(element_range)) {
                    predictor_withfallback = &fallback_predictor;
                }
                for (auto element = element_range->begin(); element != element_range->end(); ++element) {
                    auto eb = quantizer_eb.recover(*(quant_inds_eb_pos++));
                    *element = quantizer.recover(predictor_withfallback->predict(element), *(quant_inds_pos++), eb);
                }
            }
            predictor.postdecompress_data(block_range->begin());
            quantizer.postdecompress_data();
            return dec_data;
        }

        void save(uchar *&c) {
            write(global_dimensions.data(), N, c);
            write(block_size, c);

            predictor.save(c);
            quantizer_eb.save(c);
            quantizer.save(c);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            read(global_dimensions.data(), N, c, remaining_length);
            num_elements = 1;
            for (const auto &d: global_dimensions) {
                num_elements *= d;
                //std::cout << d << " ";
            }
            std::cout << std::endl;
            read(block_size, c, remaining_length);
            predictor.load(c, remaining_length);
            quantizer_eb.load(c, remaining_length);
            quantizer.load(c, remaining_length);
        }

        size_t size_est() {
            return quantizer.size_est() + sizeof(T) * num_elements;
        }

        void print() {}

        void clear() {
            predictor.clear();
            fallback_predictor.clear();
            quantizer.clear();
            quantizer_eb.clear();
        }

        int get_radius() const { return quantizer.get_radius(); }

        size_t get_num_elements() const { return num_elements; };

    private:
        Predictor predictor;
        LorenzoPredictor<T, N, 1> fallback_predictor;
        Quantizer_EB quantizer_eb;
        Quantizer quantizer;
        std::shared_ptr<concepts::QoIInterface<T, N>> qoi;
        uint block_size;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
    };

    template<class T, uint N, class Predictor, class Quantizer, class Quantizer_EB>
    SZQoIFrontend<T, N, Predictor, Quantizer, Quantizer_EB>
    make_sz_qoi_frontend(const Config &conf, Predictor predictor, Quantizer quantizer, Quantizer_EB quantizer_eb, std::shared_ptr<concepts::QoIInterface<T, N>> qoi) {
        return SZQoIFrontend<T, N, Predictor, Quantizer, Quantizer_EB>(conf, predictor, quantizer, quantizer_eb, qoi);
    }
}

#endif
