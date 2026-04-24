/**
 * @file BitplaneEncoder.hpp
 * @ingroup Encoder
 * @brief Bit-plane encoder for integer bins.
 *
 * Sweeps MSB→LSB across all bins; for each plane writes one bit per bin
 * packed 8 bins per byte. Sign-vector path is emitted only when at least one
 * input bin is negative (avoiding the overhead for already non-negative
 * inputs such as shifted quantization indices).
 *
 * Pairs naturally with a downstream lossless stage that can exploit the
 * bit-plane structure (e.g., Zstd compresses sparse top planes efficiently).
 */

#ifndef SZ3_BITPLANE_ENCODER_HPP
#define SZ3_BITPLANE_ENCODER_HPP

#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

template <class T>
class BitplaneEncoder : public concepts::EncoderInterface<T> {
    static_assert(std::is_integral<T>::value, "BitplaneEncoder only supports integer types.");

    using UT = typename std::make_unsigned<T>::type;

   public:
    void preprocess_encode(const std::vector<T>& bins, int /*stateNum*/) override {
        n_ = bins.size();
        if (n_ == 0) {
            // Empty input: nothing to encode. Decode produces an empty vector.
            num_planes_ = 0;
            has_signs_ = false;
            return;
        }
        UT max_abs = 0;
        bool has_neg = false;
        for (size_t i = 0; i < n_; i++) {
            const T v = bins[i];
            UT a;
            if (v < 0) {
                has_neg = true;
                // Use unsigned arithmetic to safely negate (avoids -INT_MIN UB).
                a = static_cast<UT>(0) - static_cast<UT>(v);
            } else {
                a = static_cast<UT>(v);
            }
            if (a > max_abs) max_abs = a;
        }
        has_signs_ = has_neg;
        // Reserve a top bit for hosts where T is signed. The shift must stay below
        // the type's signed bit width — max p is 8*sizeof(T)-2 for signed T.
        const int max_p = static_cast<int>(8 * sizeof(T)) - (std::is_signed<T>::value ? 2 : 1);
        num_planes_ = 0;
        for (UT m = max_abs; m > 0; m >>= 1) num_planes_++;
        if (num_planes_ > max_p + 1) num_planes_ = max_p + 1;
    }

    size_t encode(const std::vector<T>& bins, uchar*& bytes) override {
        if (bins.size() != n_) {
            throw std::runtime_error("BitplaneEncoder: encode size mismatch with preprocess_encode.");
        }
        uchar* const start = bytes;
        if (n_ == 0) return 0;
        const size_t plane_bytes = (n_ + 7) / 8;

        if (has_signs_) {
            std::memset(bytes, 0, plane_bytes);
            for (size_t i = 0; i < n_; i++) {
                if (bins[i] < 0) bytes[i >> 3] |= static_cast<uchar>(1u << (i & 7));
            }
            bytes += plane_bytes;
        }

        for (int p = num_planes_ - 1; p >= 0; p--) {
            std::memset(bytes, 0, plane_bytes);
            const UT mask = static_cast<UT>(UT{1} << p);
            for (size_t i = 0; i < n_; i++) {
                UT a;
                const T v = bins[i];
                if (v < 0) {
                    a = static_cast<UT>(0) - static_cast<UT>(v);
                } else {
                    a = static_cast<UT>(v);
                }
                if (a & mask) bytes[i >> 3] |= static_cast<uchar>(1u << (i & 7));
            }
            bytes += plane_bytes;
        }
        return static_cast<size_t>(bytes - start);
    }

    std::vector<T> decode(const uchar*& bytes, size_t targetLength) override {
        if (targetLength != n_) {
            throw std::runtime_error("BitplaneEncoder: decode targetLength does not match saved bin count.");
        }
        std::vector<T> bins(n_, 0);
        if (n_ == 0) return bins;
        const size_t plane_bytes = (n_ + 7) / 8;

        std::vector<uint8_t> signs;
        if (has_signs_) {
            signs.assign(n_, 0);
            for (size_t i = 0; i < n_; i++) {
                signs[i] = static_cast<uint8_t>((bytes[i >> 3] >> (i & 7)) & 1u);
            }
            bytes += plane_bytes;
        }

        std::vector<UT> mags(n_, 0);
        for (int p = num_planes_ - 1; p >= 0; p--) {
            const UT bit = static_cast<UT>(UT{1} << p);
            for (size_t i = 0; i < n_; i++) {
                if ((bytes[i >> 3] >> (i & 7)) & 1u) mags[i] |= bit;
            }
            bytes += plane_bytes;
        }

        for (size_t i = 0; i < n_; i++) {
            if (has_signs_ && signs[i]) {
                bins[i] = static_cast<T>(static_cast<UT>(0) - mags[i]);
            } else {
                bins[i] = static_cast<T>(mags[i]);
            }
        }
        return bins;
    }

    void save(uchar*& c) override {
        write(n_, c);
        write(num_planes_, c);
        const uint8_t flags = has_signs_ ? 1 : 0;
        write(flags, c);
    }

    void load(const uchar*& c, size_t& remaining_length) override {
        read(n_, c, remaining_length);
        read(num_planes_, c, remaining_length);
        uint8_t flags = 0;
        read(flags, c, remaining_length);
        has_signs_ = (flags & 1u) != 0;
    }

    void preprocess_decode() override {}
    void postprocess_encode() override {}
    void postprocess_decode() override {}

    // Size of metadata + encoded payload. SZGenericCompressor uses this to size
    // the intermediate buffer; under-reporting can mask real overruns later.
    size_t size_est() override {
        const size_t header = sizeof(n_) + sizeof(num_planes_) + sizeof(uint8_t);
        if (n_ == 0) return header;
        const size_t plane_bytes = (n_ + 7) / 8;
        const size_t total_planes = static_cast<size_t>(num_planes_) + (has_signs_ ? 1 : 0);
        return header + total_planes * plane_bytes;
    }

   private:
    size_t n_ = 0;
    int num_planes_ = 0;
    bool has_signs_ = false;
};

}  // namespace SZ3

#endif
