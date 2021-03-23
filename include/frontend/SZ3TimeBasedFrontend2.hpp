#ifndef SZ3_TIMEBASED_FRONTEND
#define SZ3_TIMEBASED_FRONTEND

#include <utils/FileUtil.h>
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
    class SZ3TimeBasedFrontend : public concepts::FrontendInterface<T, N> {
    public:

        SZ3TimeBasedFrontend(const Config<T, N> &conf, Predictor predictor, Quantizer quantizer) :
                fallback_predictor(LorenzoPredictor<T, N - 1, 1>(conf.eb)),
                predictor(predictor),
                quantizer(quantizer),
                block_size(conf.block_size),
                stride(conf.stride),
                global_dimensions(conf.dims),
                num_elements(conf.num) {
        }

        ~SZ3TimeBasedFrontend() = default;

        std::vector<int> compress(T *data) {
            std::vector<int> quant_inds(num_elements);
            size_t quant_count = 0;
            size_t num;
            auto data_time0 = readfile<T>("../x0.decompressed", num);
            assert(num == global_dimensions[1]);
            for (size_t j = 0; j < global_dimensions[1]; j++) {
                quant_inds[quant_count++] = quantizer.quantize_and_overwrite(data[j], data_time0[j]);
            }

            for (size_t j = 0; j < global_dimensions[1]; j++) {
                for (int i = 1; i < global_dimensions[0]; i++) {
                    size_t idx = i * global_dimensions[1] + j;
                    size_t idx_prev = (i - 1) * global_dimensions[1] + j;
                    quant_inds[quant_count++] = quantizer.quantize_and_overwrite(data[idx], data[idx_prev]);
                }
            }

            assert(quant_count == num_elements);
            quantizer.postcompress_data();
            return quant_inds;
        }

        T *decompress(std::vector<int> &quant_inds) {

            int const *quant_inds_pos = (int const *) quant_inds.data();
            std::array<size_t, N - 1> intra_block_dims;
            auto dec_data = new T[num_elements];

            size_t num;
            auto data_time0 = readfile<T>("../x0.decompressed", num);
            assert(num == global_dimensions[1]);
            for (size_t j = 0; j < global_dimensions[1]; j++) {
                dec_data[j] = quantizer.recover(data_time0[j], *(quant_inds_pos++));
            }

            for (size_t j = 0; j < global_dimensions[1]; j++) {
                for (int i = 1; i < global_dimensions[0]; i++) {
                    size_t idx = i * global_dimensions[1] + j;
                    size_t idx_prev = (i - 1) * global_dimensions[1] + j;
                    dec_data[idx] = quantizer.recover(dec_data[idx_prev], *(quant_inds_pos++));
                }
            }

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
        LorenzoPredictor<T, N - 1, 1> fallback_predictor;
        Quantizer quantizer;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;

    };

    template<class T, uint N, class Predictor, class Quantizer>
    SZ3TimeBasedFrontend<T, N, Predictor, Quantizer>
    make_sz3_timebased_frontend(const Config<T, N> &conf, Predictor predictor, Quantizer quantizer) {
        return SZ3TimeBasedFrontend<T, N, Predictor, Quantizer>(conf, predictor, quantizer);
    }
}

#endif
