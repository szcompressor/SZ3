/**
 * @file FixedPointQuantizer.hpp
 * @ingroup Quantizer
 * @brief Fixed-point (`ldexp`-style) scalar quantizer for floating-point input.
 *
 * Given a target bit-width `num_bits` and a one-shot calibration against the
 * data's max absolute value, the quantizer produces an integer
 *   fp = round( val * 2^(num_bits - 1 - level_exp) )
 * where `level_exp = ceil(log2(max_abs))`. Returned bins are shifted to
 * `[0, 2^num_bits)` so they satisfy SZ3's "first == 0" output-range contract
 * (matching `LinearQuantizer`). `pred` is ignored — this is a plain scalar
 * quantizer with no prediction.
 *
 * Usage:
 *   FixedPointQuantizer<double> q(num_bits);  // num_bits ∈ [2, 62]
 *   q.calibrate(max_abs);                     // set scale once from data range
 *   bin = q.quantize_and_overwrite(val, 0);
 *   ...
 *   val = q.recover(0, bin);
 */

#ifndef SZ3_FIXED_POINT_QUANTIZER_HPP
#define SZ3_FIXED_POINT_QUANTIZER_HPP

#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <type_traits>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

template <class T>
class FixedPointQuantizer : public concepts::QuantizerInterface<T, int64_t> {
    static_assert(std::is_floating_point<T>::value,
                  "FixedPointQuantizer requires a floating-point input type.");

   public:
    explicit FixedPointQuantizer(int num_bits = 32) : num_bits_(num_bits) {
        if (num_bits_ < 2 || num_bits_ > 62) {
            throw std::invalid_argument("FixedPointQuantizer: num_bits must be in [2, 62].");
        }
        offset_ = static_cast<int64_t>(1) << (num_bits_ - 1);
    }

    /**
     * @brief Set the quantization scale from the maximum absolute value in the data.
     *
     * `level_exp = ceil(log2(max_abs))` (0 if max_abs <= 0). The scale is
     * `2^(num_bits - 1 - level_exp)` so that |val| ≤ max_abs ⇒ |fp| < 2^(num_bits-1).
     * Caller MUST invoke this before `quantize_and_overwrite()`.
     */
    void calibrate(double max_abs) {
        if (!(max_abs > 0.0)) {
            level_exp_ = 0;
            scale_ = 1.0;
            return;
        }
        int e = 0;
        std::frexp(max_abs, &e);  // max_abs = mantissa * 2^e, mantissa ∈ [0.5, 1)
        level_exp_ = e;
        scale_ = std::ldexp(1.0, num_bits_ - 1 - level_exp_);
        if (!std::isfinite(scale_) || scale_ <= 0.0) {
            throw std::overflow_error(
                "FixedPointQuantizer::calibrate produced a non-finite or non-positive scale.");
        }
    }

    /// Theoretical max abs reconstruction error for the current calibration.
    double max_abs_error() const {
        if (scale_ <= 0.0) return 0.0;
        return 0.5 / scale_;
    }

    int num_bits() const { return num_bits_; }
    int level_exp() const { return level_exp_; }
    double scale() const { return scale_; }

    ALWAYS_INLINE int64_t quantize_and_overwrite(T& data, T /*pred*/) override {
        if (scale_ == 0.0) {
            throw std::runtime_error(
                "FixedPointQuantizer: must call calibrate() before quantize_and_overwrite().");
        }
        const double scaled = static_cast<double>(data) * scale_;
        int64_t fp = static_cast<int64_t>(std::lround(scaled));
        // Clamp to representable range; protects against rounding up at the boundary.
        const int64_t lo = -offset_;
        const int64_t hi = offset_ - 1;
        if (fp < lo) fp = lo;
        if (fp > hi) fp = hi;
        data = static_cast<T>(static_cast<double>(fp) / scale_);
        return fp + offset_;  // shift to [0, 2^num_bits)
    }

    ALWAYS_INLINE T recover(T /*pred*/, int64_t quant_index) override {
        const int64_t fp = quant_index - offset_;
        return static_cast<T>(static_cast<double>(fp) / scale_);
    }

    int64_t force_save_unpred(T /*ori*/) override {
        // Fixed-point quantizers cover their entire calibrated range with finite precision;
        // there is no out-of-range "unpredictable" path. Returning 0 mirrors LinearQuantizer's
        // sentinel; callers that depend on round-tripping out-of-range samples must clamp first.
        return 0;
    }

    std::pair<int64_t, int64_t> get_out_range() const override {
        return std::make_pair(static_cast<int64_t>(0), 2 * offset_);
    }

    void save(uchar*& c) const override {
        write(uid_, c);
        write(num_bits_, c);
        write(level_exp_, c);
        write(scale_, c);
    }

    void load(const uchar*& c, size_t& remaining_length) override {
        uchar uid_read = 0;
        read(uid_read, c, remaining_length);
        if (uid_read != uid_) {
            throw std::invalid_argument("FixedPointQuantizer uid mismatch");
        }
        read(num_bits_, c, remaining_length);
        read(level_exp_, c, remaining_length);
        read(scale_, c, remaining_length);
        offset_ = static_cast<int64_t>(1) << (num_bits_ - 1);
    }

    void print() override {
        printf("[FixedPointQuantizer] num_bits=%d level_exp=%d scale=%.8g (max_err=%.4g)\n",
               num_bits_, level_exp_, scale_, max_abs_error());
    }

   private:
    int num_bits_;
    int level_exp_ = 0;
    int64_t offset_ = 0;
    double scale_ = 0.0;  // 0 = uncalibrated
    static constexpr uchar uid_ = 0b11;  // distinct from LinearQuantizer (0b10)
};

}  // namespace SZ3

#endif
