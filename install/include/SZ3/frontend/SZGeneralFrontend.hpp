#ifndef SZ3_FRONTEND
#define SZ3_FRONTEND
/**
 * This module is the implementation of general frontend in SZ3
 */

#include "Frontend.hpp"
#include "SZ3/def.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ {


    template<class T, uint N, class Predictor, class Quantizer>
    class SZGeneralFrontend : public concepts::FrontendInterface<T, N> {
    public:

        SZGeneralFrontend(const Config &conf, Predictor predictor, Quantizer quantizer) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(conf.absErrorBound)),
                predictor(predictor),
                quantizer(quantizer),
                block_size(conf.blockSize),
                num_elements(conf.num) {
            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
            //sprintf(block_info_path, "%s.blockinfo.dat", conf.input_path);
            strcpy(input_path,  conf.input_path);
        }

        ~SZGeneralFrontend() = default;

        std::vector<int> compress(T *data) {
            std::vector<int> quant_inds(num_elements);
            std::cout<< "num-elements " << num_elements <<std::endl;
            std::vector<T> my_pred(num_elements);
            std::vector<int> my_quant(num_elements);


            auto block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    data, std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);

            auto element_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                    data, std::begin(global_dimensions), std::end(global_dimensions), 1, 0);

            predictor.pass_input_path(input_path);
            predictor.precompress_data(block_range->begin());
            quantizer.precompress_data();

            size_t quant_count = 0;
            concepts::PredictorInterface<T, N> *tmp;

            for (auto block = block_range->begin(); block != block_range->end(); ++block) {

                element_range->update_block_range(block, block_size);

                concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;

                if (!predictor.precompress_block(element_range)) {
                    predictor_withfallback = &fallback_predictor;
                }
                predictor_withfallback->precompress_block_commit();
             
                for (auto element = element_range->begin(); element != element_range->end(); ++element) {
                    auto pred = predictor_withfallback->predict(element);
                    my_pred[element.get_offset()]=pred;
                    quant_inds[quant_count++] = quantizer.quantize_and_overwrite(
                            *element, pred);
                    int quant_idx=quant_count-1;
                    my_quant[element.get_offset()]=quant_inds[quant_idx];
                }


            }

            // output the predictor selection in composed predictor. empty for other predictors. 
            
            predictor.write_block_info(input_path);
            char quant_path[1400];
            sprintf(quant_path, "%s.quant_inds", input_path);
            writefile("quant.dat",my_quant.data(),num_elements);
            writefile("pred.dat", my_pred.data(),num_elements);
            std::cout<<"quant size "<< quant_inds.size() << std::endl;

            predictor.postcompress_data(block_range->begin());  // empty function for both lorenzo and regression 
            quantizer.postcompress_data();
            return quant_inds;
        }

        T *decompress(std::vector<int> &quant_inds, T *dec_data) {

            int const *quant_inds_pos = (int const *) quant_inds.data();
            std::array<size_t, N> intra_block_dims;
//            auto dec_data = new T[num_elements];
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
                    *element = quantizer.recover(predictor_withfallback->predict(element), *(quant_inds_pos++));
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
            quantizer.save(c);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            read(global_dimensions.data(), N, c, remaining_length);
            num_elements = 1;
            for (const auto &d: global_dimensions) {
                num_elements *= d;
            }
            read(block_size, c, remaining_length);
            predictor.load(c, remaining_length);
            quantizer.load(c, remaining_length);
        }

        size_t size_est() {
            return quantizer.size_est() + sizeof(T) * num_elements;
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
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
        char input_path[1400];
        

    };

    template<class T, uint N, class Predictor, class Quantizer>
    SZGeneralFrontend<T, N, Predictor, Quantizer>
    make_sz_general_frontend(const Config &conf, Predictor predictor, Quantizer quantizer) {
        return SZGeneralFrontend<T, N, Predictor, Quantizer>(conf, predictor, quantizer);
    }
}

#endif
