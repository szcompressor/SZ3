/**
 * @file BitTruncationQuantizer.hpp
 * @ingroup Quantizer
 * @brief Bit-truncation scalar quantizer: keep the top `keep_bytes` bytes of
 *        each floating-point value's IEEE-754 bit pattern, zero out the rest.
 *
 * Output bin is the truncated bit pattern reinterpreted as an integer, so the
 * round-trip is `recover(quantize(x)) == x_with_low_bits_zeroed`. Pairs
 * naturally with `BypassEncoder` + `Lossless_zstd` (downstream zstd
 * compresses the resulting zero-padded buffer effectively); also composes
 * with `BitplaneEncoder` since the high planes are now mostly identical
 * across the dataset.
 *
 * Notes:
 *  - `pred` is ignored — this is a pure scalar quantizer with no prediction.
 *  - Max abs error is bounded by 2^(e - mantissa_bits_kept) where e is the
 *    floating-point exponent of each value. Unlike a fixed-eb quantizer,
 *    error scales with magnitude.
 *  - For T=float: `keep_bytes ∈ [1, 4]`. For T=double: `keep_bytes ∈ [1, 8]`.
 *    Output bin type is `int64_t` so both fit.
 */

#ifndef SZ3_BIT_TRUNCATION_QUANTIZER_HPP
#define SZ3_BIT_TRUNCATION_QUANTIZER_HPP

#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

template <class T>
class BitTruncationQuantizer : public concepts::QuantizerInterface<T, int64_t> {
    static_assert(std::is_floating_point<T>::value,
                  "BitTruncationQuantizer requires a floating-point input type.");
    static_assert(sizeof(T) <= sizeof(uint64_t),
                  "BitTruncationQuantizer only supports T up to 8 bytes.");

   public:
    explicit BitTruncationQuantizer(int keep_bytes = sizeof(T) / 2) : keep_bytes_(keep_bytes) {
        validate_and_compute_mask();
    }

    ALWAYS_INLINE int64_t quantize_and_overwrite(T& data, T /*pred*/) override {
        uint64_t bits = bits_from(data);
        bits &= mask_;
        data = value_from(bits);
        return static_cast<int64_t>(bits);
    }

    ALWAYS_INLINE T recover(T /*pred*/, int64_t q) override {
        return value_from(static_cast<uint64_t>(q));
    }

    int64_t force_save_unpred(T ori) override {
        // Bit truncation has no out-of-range bucket; just truncate `ori`.
        T tmp = ori;
        return quantize_and_overwrite(tmp, T(0));
    }

    std::pair<int64_t, int64_t> get_out_range() const override {
        return std::make_pair(static_cast<int64_t>(0), std::numeric_limits<int64_t>::max());
    }

    void save(uchar*& c) const override {
        write(uid_, c);
        write(keep_bytes_, c);
    }

    void load(const uchar*& c, size_t& remaining_length) override {
        uchar uid_read = 0;
        read(uid_read, c, remaining_length);
        if (uid_read != uid_) {
            throw std::invalid_argument("BitTruncationQuantizer uid mismatch");
        }
        read(keep_bytes_, c, remaining_length);
        validate_and_compute_mask();
    }

    void print() override {
        printf("[BitTruncationQuantizer<%zuB>] keep_bytes=%d mask=0x%016llx\n", sizeof(T), keep_bytes_,
               static_cast<unsigned long long>(mask_));
    }

    int keep_bytes() const { return keep_bytes_; }

   private:
    void validate_and_compute_mask() {
        if (keep_bytes_ < 1 || keep_bytes_ > static_cast<int>(sizeof(T))) {
            throw std::invalid_argument("BitTruncationQuantizer: keep_bytes out of [1, sizeof(T)].");
        }
        const int drop_bits = (static_cast<int>(sizeof(T)) - keep_bytes_) * 8;
        mask_ = (drop_bits >= 64) ? uint64_t{0} : ~((uint64_t{1} << drop_bits) - 1);
        // Limit mask to the value's actual bit width (T may be 4 bytes).
        if (sizeof(T) < sizeof(uint64_t)) {
            const uint64_t value_mask = (uint64_t{1} << (sizeof(T) * 8)) - 1;
            mask_ &= value_mask;
        }
    }

    static uint64_t bits_from(T data) {
        uint64_t bits = 0;
        std::memcpy(&bits, &data, sizeof(T));
        return bits;
    }

    static T value_from(uint64_t bits) {
        T out;
        std::memcpy(&out, &bits, sizeof(T));
        return out;
    }

    int keep_bytes_;
    uint64_t mask_ = 0;
    static constexpr uchar uid_ = 0b101;  // distinct from Linear (0b10), QuadraticLevel/FixedPoint (0b11)
};

}  // namespace SZ3

#endif
