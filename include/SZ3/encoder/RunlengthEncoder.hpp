#ifndef SZ3_RUNLENGTH_ENCODER_HPP
#define SZ3_RUNLENGTH_ENCODER_HPP

#include <vector>

#include "Encoder.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

template <class T>
class RunlengthEncoder : public concepts::EncoderInterface<T> {
   public:
    void preprocess_encode(const std::vector<T> &bins, int stateNum) override {}

    size_t encode(const std::vector<T> &bins, uchar *&bytes) override {
        auto bytespos = bytes;
        int max = 0;
        size_t s = 0;
        for (size_t i = 1; i < bins.size(); i++) {
            if (bins[i] != bins[i - 1]) {
                write(bins[i - 1], bytes);
                write(int(i - s), bytes);
                if (int(i - s) > max) {
                    max = int(i - s);
                }
                s = i;
            }
        }
        write(bins[bins.size() - 1], bytes);
        write(int(bins.size() - s), bytes);
        return bytes - bytespos;
    }

    void postprocess_encode() override {}

    void preprocess_decode() override {}

    std::vector<T> decode(const uchar *&bytes, size_t targetLength) override {
        std::vector<T> bins(targetLength, 0);
        T value;
        int cnt;
        for (size_t i = 0; i < bins.size();) {
            read(value, bytes);
            read(cnt, bytes);
            if (i + cnt > bins.size()) {
                throw std::runtime_error("Decoded length exceeds targetLength");
            }
            for (size_t j = i; j < i + cnt; j++) {
                bins[j] = value;
            }
            i += cnt;
        }
        return bins;
    }

    void postprocess_decode() override {}

    void save(uchar *&c) override {}

    void load(const uchar *&c, size_t &remaining_length) override {}
};
}  // namespace SZ3
#endif
