#ifndef _SZ_NO_PREDICTION_DECOMPOSITION_HPP
#define _SZ_NO_PREDICTION_DECOMPOSITION_HPP

#include "Decomposition.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/def.hpp"

namespace SZ3 {
    template<class T, uint N, class Quantizer>
    class NoPredictionDecomposition : public concepts::DecompositionInterface<T, N> {
    public:


        NoPredictionDecomposition(const Config &conf, Quantizer quantizer) : quantizer(quantizer) {
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
        }

        T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) {
            for (size_t i = 0; i < conf.num; i++) {
                dec_data[i] = quantizer.recover(0, quant_inds[i]);
            }
            quantizer.postdecompress_data();
            return dec_data;
        }

        std::vector<int> compress(const Config &conf, T *data) {
            std::vector<int> quant_inds(conf.num);
            for (size_t i = 0; i < conf.num; i++) {
                quant_inds[i] = quantizer.quantize_and_overwrite(data[i], 0);
            }
            quantizer.postcompress_data();
            return quant_inds;
        }

        void save(uchar *&c) {
            quantizer.save(c);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            quantizer.load(c, remaining_length);
        }

    private:
        Quantizer quantizer;
    };

    template<class T, uint N, class Quantizer>
    NoPredictionDecomposition<T, N, Quantizer>
    make_decomposition_noprediction(const Config &conf, Quantizer quantizer) {
        return NoPredictionDecomposition<T, N, Quantizer>(conf, quantizer);
    }

};


#endif

