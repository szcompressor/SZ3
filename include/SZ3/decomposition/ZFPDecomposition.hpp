#ifndef SZ3_NO_PREDICTION_DECOMPOSITION_HPP
#define SZ3_NO_PREDICTION_DECOMPOSITION_HPP

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/zfp/zfpcodec.h"
#include "SZ3/utils/zfp/zfpcodec1.h"
#include "SZ3/utils/zfp/zfpcodec2.h"
#include "SZ3/utils/zfp/zfpcodec3.h"

namespace SZ3 {
template <class T, uint N, class Quantizer>
class NoPredictionDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    NoPredictionDecomposition(const Config &conf, Quantizer quantizer) : quantizer(quantizer) {
        static_assert(std::is_base_of<concepts::QuantizerInterface<T, int>, Quantizer>::value,
                      "must implement the quantizer interface");
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override { return true; }

    std::vector<int> compress(const Config &conf, T *data) override {
        MemoryBitStream stream;
        stream.open(out, outsize);
        int expmin = INT_MIN;
        if (conf.absErrorBound > 0) {
            frexp(conf.absErrorBound, &expmin);
            expmin--;
        }
        // bool dp = (p->type == ZFP_TYPE_DOUBLE);
        const T *p = data;

        ZFP::Codec1<MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
        for (auto x = 0; x < conf.dims[0]; x += 4, p += 4) {
            codec.encode(p, 1, codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u)));
        }
        stream.flush();
        return stream.size();
    }

    void save(uchar *&c) override { quantizer.save(c); }

    void load(const uchar *&c, size_t &remaining_length) override { quantizer.load(c, remaining_length); }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

   private:
    Quantizer quantizer;
};

template <class T, uint N, class Quantizer>
NoPredictionDecomposition<T, N, Quantizer> make_decomposition_noprediction(const Config &conf, Quantizer quantizer) {
    return NoPredictionDecomposition<T, N, Quantizer>(conf, quantizer);
}

}  // namespace SZ3

#endif
