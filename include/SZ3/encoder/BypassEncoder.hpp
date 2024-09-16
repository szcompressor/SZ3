#ifndef _SZ_BYPASS_ENCODER_HPP
#define _SZ_BYPASS_ENCODER_HPP

#include <cassert>
#include <vector>

#include "Encoder.hpp"
#include "SZ3/def.hpp"

namespace SZ3 {

template <class T>
class BypassEncoder : public concepts::EncoderInterface<T> {
   public:
    void preprocess_encode(const std::vector<T> &bins, int stateNum) override {
        assert(stateNum <= 256 && "stateNum should be no more than 256.");
    }

    size_t encode(const std::vector<T> &bins, uchar *&bytes) override {
        for (auto &bin : bins) {
            *bytes++ = uchar(bin);
        }
        return 0;
    }

    void postprocess_encode() override{}

    void preprocess_decode() override{}

    std::vector<T> decode(const uchar *&bytes, size_t targetLength) override {
        std::vector<T> bins(targetLength);
        for (auto &bin : bins) {
            bin = *bytes++;
        }
        return bins;
    }

    void postprocess_decode() override{}

    void save(uchar *&c) override {}

    void load(const uchar *&c, size_t &remaining_length) override{}
};
}  // namespace SZ3
#endif
