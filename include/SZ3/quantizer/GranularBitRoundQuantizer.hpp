/**
 * @file GranularBitRoundQuantizer.hpp
 * @ingroup Quantizer
 * @brief Granular BitRound: keep `nsd` significant decimal digits by rounding away the
 *        redundant low mantissa bits of each IEEE-754 value.
 *
 * The number of retained mantissa bits is recomputed per value from that value's magnitude,
 * so with `d10(x) = floor(log10|x|) + 1` every normal finite non-zero input satisfies
 *
 *     |x_hat - x| <= 0.5 * 10^(d10(x) - nsd)
 *
 * Zero, subnormals, infinities and NaN pass through bit-exactly.
 *
 * `pred` is ignored -- this is a pure scalar quantizer. The bin is the rounded IEEE-754 bit
 * pattern reinterpreted as `uint64_t`; unsigned because a negative `double` fills all 64 bits
 * and a signed bin would break the `get_out_range().first == 0` contract.
 */

#ifndef SZ3_GRANULAR_BIT_ROUND_QUANTIZER_HPP
#define SZ3_GRANULAR_BIT_ROUND_QUANTIZER_HPP

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * @brief Granular BitRound quantizer -- keep `nsd` significant decimal digits per value.
 *
 * See the file-level comment for the exact error guarantee.
 *
 * @tparam T Floating-point data type (`float` or `double`).
 */
template <class T>
class GranularBitRoundQuantizer : public concepts::QuantizerInterface<T, uint64_t> {
    static_assert(std::is_floating_point<T>::value, "GranularBitRoundQuantizer requires a floating-point input type.");
    static_assert(sizeof(T) <= sizeof(uint64_t), "GranularBitRoundQuantizer only supports T up to 8 bytes.");

   public:
    /// Number of explicitly stored mantissa bits (23 for float, 52 for double).
    static constexpr int kMantBits = std::numeric_limits<T>::digits - 1;

    /**
     * @brief Construct a Granular BitRound quantizer.
     *
     * @param nsd Number of significant decimal digits to preserve, in [1, 18]. Values
     *            large enough to require more than `kMantBits` mantissa bits degrade
     *            gracefully to a bit-exact (lossless) pass-through: `nsd >= 8` for
     *            `float`, `nsd >= 17` for `double`.
     * @throws std::invalid_argument if `nsd` is outside [1, 18].
     */
    explicit GranularBitRoundQuantizer(int nsd = 3) : nsd_(nsd) { validate(); }

    /// @brief Number of significant decimal digits this quantizer preserves.
    int nsd() const { return nsd_; }

    /**
     * @brief Number of mantissa bits retained for a given value.
     *
     * Exposed for diagnostics and tests; this is the `m` of the file-level comment,
     * clamped to `[0, kMantBits]`. Returns `kMantBits` for values that are passed
     * through unmodified (zero, subnormal, infinity, NaN).
     *
     * @param value The value whose retained-bit count is wanted.
     * @return int Retained mantissa bit count in [0, kMantBits].
     */
    int keep_mantissa_bits(T value) const {
        if (!is_roundable(value)) {
            return kMantBits;
        }
        return keep_bits_for(std::fabs(static_cast<double>(value)));
    }

    /**
     * @brief Guaranteed absolute error bound for a given value, i.e. `0.5 * 10^(d10-nsd)`.
     *
     * Uses the same (slightly conservative) `d10` the quantizer itself uses, so near a decade
     * boundary this can be up to 10x tighter than the bound stated in guarantee (1) -- never
     * looser.
     *
     * @param value The value whose bound is wanted.
     * @return double The right-hand side of guarantee (1); 0 for pass-through values.
     */
    double max_abs_error(T value) const {
        if (!is_roundable(value)) {
            return 0.0;
        }
        const double a = std::fabs(static_cast<double>(value));
        return 0.5 * std::pow(10.0, static_cast<double>(decimal_exponent(a) - nsd_));
    }

    /**
     * @brief Round `data` to `nsd` significant decimal digits and return its bit pattern.
     *
     * @param data Value to quantize; overwritten by the reconstructed (rounded) value.
     * @param pred Ignored -- this quantizer does not use prediction.
     * @return uint64_t The rounded IEEE-754 bit pattern of `data`, zero-extended.
     */
    ALWAYS_INLINE uint64_t quantize_and_overwrite(T& data, T /*pred*/) override {
        const uint64_t bits = round_bits(data);
        data = value_from(bits);
        return bits;
    }

    /**
     * @brief Reconstruct the value from a bin produced by `quantize_and_overwrite()`.
     *
     * @param pred Ignored -- this quantizer does not use prediction.
     * @param quant_index Bin returned by `quantize_and_overwrite()`/`force_save_unpred()`.
     * @return T The reconstructed value.
     */
    ALWAYS_INLINE T recover(T /*pred*/, uint64_t quant_index) override { return value_from(quant_index); }

    /**
     * @brief Store `ori` losslessly: the returned bin is its unrounded bit pattern.
     *
     * @param ori Value to store exactly.
     * @return uint64_t Bin that recovers `ori` bit-exactly.
     */
    uint64_t force_save_unpred(T ori) override { return bits_from(ori); }

    /**
     * @brief Output bin range, `[first, second)`. `first` is 0 as required by SZ3.
     *
     * @return std::pair<uint64_t, uint64_t> `{0, 2^(8*sizeof(T))}`; for `double` the exclusive
     *         bound `2^64` is not representable, so it saturates at `UINT64_MAX`.
     */
    std::pair<uint64_t, uint64_t> get_out_range() const override {
        // 0 means "no bin range"; the encoder derives what it needs from the bins.
        return std::make_pair(uint64_t{0}, uint64_t{0});
    }

    /**
     * @brief Serialize the quantizer state (uid + nsd).
     *
     * @param c Buffer pointer; advanced past the written bytes.
     */
    void save(uchar*& c) const override {
        write(uid_, c);
        write(nsd_, c);
    }

    /**
     * @brief Deserialize the quantizer state written by `save()`.
     *
     * @param c Buffer pointer; advanced past the bytes read.
     * @param remaining_length Remaining readable bytes; decremented accordingly.
     * @throws std::invalid_argument on uid mismatch or an out-of-range `nsd`.
     */
    void load(const uchar*& c, size_t& remaining_length) override {
        uchar uid_read = 0;
        read(uid_read, c, remaining_length);
        if (uid_read != uid_) {
            throw std::invalid_argument("GranularBitRoundQuantizer uid mismatch");
        }
        read(nsd_, c, remaining_length);
        validate();
    }

    /// @brief Print a one-line summary of the quantizer configuration.
    void print() override {
        printf("[GranularBitRoundQuantizer<%zuB>] nsd=%d (mantissa bits available=%d)\n", sizeof(T), nsd_, kMantBits);
    }

   private:
    void validate() const {
        if (nsd_ < 1 || nsd_ > 18) {
            throw std::invalid_argument("GranularBitRoundQuantizer: nsd must be in [1, 18].");
        }
    }

    /// Values that are neither rounded nor modified: zero, subnormal, infinity, NaN.
    /// Note `std::isnormal` is evaluated on `T` itself, so a subnormal `float` is
    /// recognised as subnormal even though it is a normal `double` once promoted.
    static bool is_roundable(T value) { return std::isnormal(value); }

    /**
     * @brief `d10` of the file-level comment: the integer `n` with `10^(n-1) <= a < 10^n`.
     *
     * Biased downward by 1e-12 in log space so that `log10` rounding can only ever make
     * the quantizer keep more bits than required, never fewer.
     */
    static int decimal_exponent(double a) { return static_cast<int>(std::floor(std::log10(a) - 1e-12)) + 1; }

    /// Retained mantissa bit count for `a = |x|`, assumed normal and non-zero.
    int keep_bits_for(double a) const {
        int e = 0;
        std::frexp(a, &e);  // a = f * 2^e, f in [0.5, 1)
        const int d10 = decimal_exponent(a);
        // m >= e - 1 - (d10 - nsd) * log2(10)
        const double m = static_cast<double>(e - 1) - static_cast<double>(d10 - nsd_) * kLog2Of10;
        if (!(m > 0.0)) {
            return 0;
        }
        if (m >= static_cast<double>(kMantBits)) {
            return kMantBits;
        }
        return static_cast<int>(std::ceil(m));
    }

    /// Round `data` to `nsd` significant digits, returning the resulting bit pattern.
    uint64_t round_bits(T data) const {
        const uint64_t bits = bits_from(data);
        if (!is_roundable(data)) {
            return bits;
        }
        const int drop = kMantBits - keep_bits_for(std::fabs(static_cast<double>(data)));
        if (drop <= 0) {
            return bits;
        }
        const uint64_t half = uint64_t{1} << (drop - 1);
        const uint64_t mask = ~((uint64_t{1} << drop) - 1) & kValueMask;
        // Round to nearest, ties to even. Carry out of the mantissa field correctly
        // increments the exponent field, which is the desired IEEE-754 behaviour.
        const uint64_t rounded = (bits + half - 1 + ((bits >> drop) & uint64_t{1})) & mask;
        if (!std::isfinite(value_from(rounded))) {
            // Rounding up overflowed the top binade; keep the value exact instead.
            return bits;
        }
        return rounded;
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

    /// Mask covering the bits actually occupied by a `T` inside the `uint64_t` carrier.
    static constexpr uint64_t kValueMask =
        (sizeof(T) < sizeof(uint64_t)) ? ((uint64_t{1} << (sizeof(T) * 8)) - 1) : ~uint64_t{0};

    /// log2(10), the number of binary digits per decimal digit.
    static constexpr double kLog2Of10 = 3.321928094887362;

    int nsd_;
    // Distinct from Linear (0b10), FixedPoint (0b11), BitTruncation (0b101), Level (0b1010).
    static constexpr uchar uid_ = 0b1000;
};

}  // namespace SZ3

#endif
