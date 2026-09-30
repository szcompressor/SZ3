#ifndef SZ3_SEGMENTED_ENCODER_HPP
#define SZ3_SEGMENTED_ENCODER_HPP

#include <stdexcept>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * Bins that are a sequence of segments [n, n bins]: each segment gets its own Encoder, so a decomposition that
 * concatenates streams of different statistics gets a code per stream. Decoding returns the bins as they were.
 */
template <class Encoder>
class SegmentedEncoder : public concepts::EncoderInterface<int> {
   public:
    void preprocess_encode(const std::vector<int> &bins, int stateNum) override {
        segments_.clear();
        encoders_.clear();
        for (size_t k = 0; k < bins.size();) {
            const size_t n = size_t(bins[k]);
            if (bins[k] < 0 || n > bins.size() - k - 1) throw std::invalid_argument("SZ3: bins are not segments");
            segments_.emplace_back(bins.begin() + k + 1, bins.begin() + k + 1 + n);
            encoders_.emplace_back();
            if (n) encoders_.back().preprocess_encode(segments_.back(), stateNum);
            k += n + 1;
        }
    }

    size_t encode(const std::vector<int> & /*bins*/, uchar *&bytes) override {
        uchar *const start = bytes;
        for (size_t s = 0; s < segments_.size(); s++)
            if (!segments_[s].empty()) encoders_[s].encode(segments_[s], bytes);
        return bytes - start;
    }

    void postprocess_encode() override { segments_.clear(); }

    void preprocess_decode() override {}

    std::vector<int> decode(const uchar *&bytes, size_t targetLength, size_t &remaining_length) override {
        uint64_t total = encoders_.size();
        for (const uint64_t n : segment_sizes_) total = n > targetLength - total ? targetLength + 1 : total + n;
        if (total != targetLength) throw std::out_of_range("SZ3: segment sizes do not match the bin count");
        std::vector<int> out;
        out.reserve(targetLength);
        for (size_t s = 0; s < encoders_.size(); s++) {
            out.push_back(int(segment_sizes_[s]));
            if (!segment_sizes_[s]) continue;
            const auto v = encoders_[s].decode(bytes, segment_sizes_[s], remaining_length);
            out.insert(out.end(), v.begin(), v.end());
        }
        if (out.size() != targetLength) throw std::out_of_range("SZ3: segment sizes do not match the bin count");
        return out;
    }

    void postprocess_decode() override {}

    void save(uchar *&c) override {
        write(uint32_t(encoders_.size()), c);
        for (size_t s = 0; s < encoders_.size(); s++) {
            uint64_t n = segments_[s].size();  // in 7-bit groups, low first
            for (; n >= 128; n >>= 7) *c++ = uchar(n | 128);
            *c++ = uchar(n);
            if (!segments_[s].empty()) encoders_[s].save(c);
        }
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        uint32_t count = 0;
        read(count, c, remaining_length);
        if (count > remaining_length) throw std::out_of_range("SZ3: segment count exceeds the buffer");
        encoders_.assign(count, Encoder());
        segment_sizes_.resize(count);
        for (uint32_t s = 0; s < count; s++) {
            segment_sizes_[s] = 0;
            for (uint8_t byte = 128, shift = 0; byte >= 128 && shift < 64; shift += 7) {
                read(byte, c, remaining_length);
                segment_sizes_[s] |= uint64_t(byte & 127) << shift;
            }
            if (segment_sizes_[s]) encoders_[s].load(c, remaining_length);
        }
    }

    size_t size_est() override {
        size_t bound = sizeof(uint32_t);
        for (auto &encoder : encoders_) bound += 10 + encoder.size_est();  // 10: a varint of 64 bits
        return bound;
    }

   private:
    std::vector<Encoder> encoders_;
    std::vector<std::vector<int>> segments_;
    std::vector<uint64_t> segment_sizes_;
};

}  // namespace SZ3
#endif
