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
        segs_.clear();
        enc_.clear();
        for (size_t k = 0; k < bins.size();) {
            const size_t n = size_t(bins[k]);
            if (bins[k] < 0 || n > bins.size() - k - 1) throw std::invalid_argument("SZ3: bins are not segments");
            segs_.emplace_back(bins.begin() + k + 1, bins.begin() + k + 1 + n);
            enc_.emplace_back();
            if (n) enc_.back().preprocess_encode(segs_.back(), stateNum);
            k += n + 1;
        }
    }

    size_t encode(const std::vector<int> & /*bins*/, uchar *&bytes) override {
        uchar *const start = bytes;
        for (size_t s = 0; s < segs_.size(); s++)
            if (!segs_[s].empty()) enc_[s].encode(segs_[s], bytes);
        return bytes - start;
    }

    void postprocess_encode() override { segs_.clear(); }

    void preprocess_decode() override {}

    std::vector<int> decode(const uchar *&bytes, size_t targetLength, size_t &remaining_length) override {
        std::vector<int> out;
        out.reserve(targetLength);
        for (size_t s = 0; s < enc_.size(); s++) {
            out.push_back(int(sizes_[s]));
            if (!sizes_[s]) continue;
            const auto v = enc_[s].decode(bytes, sizes_[s], remaining_length);
            out.insert(out.end(), v.begin(), v.end());
        }
        if (out.size() != targetLength) throw std::out_of_range("SZ3: segment sizes do not match the bin count");
        return out;
    }

    void postprocess_decode() override {}

    void save(uchar *&c) override {
        write(uint32_t(enc_.size()), c);
        for (size_t s = 0; s < enc_.size(); s++) {
            write(uint64_t(segs_[s].size()), c);
            if (!segs_[s].empty()) enc_[s].save(c);
        }
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        uint32_t n = 0;
        read(n, c, remaining_length);
        if (n > remaining_length / sizeof(uint64_t)) throw std::out_of_range("SZ3: segment count exceeds the buffer");
        enc_.assign(n, Encoder());
        sizes_.resize(n);
        for (uint32_t s = 0; s < n; s++) {
            read(sizes_[s], c, remaining_length);
            if (sizes_[s]) enc_[s].load(c, remaining_length);
        }
    }

    size_t size_est() override {
        size_t e = sizeof(uint32_t);
        for (auto &enc : enc_) e += sizeof(uint64_t) + enc.size_est();
        return e;
    }

   private:
    std::vector<Encoder> enc_;
    std::vector<std::vector<int>> segs_;
    std::vector<uint64_t> sizes_;
};

}  // namespace SZ3
#endif
