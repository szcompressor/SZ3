/**
 * @file ScalarQuantizer.hpp
 * @ingroup Quantizer
 */

#ifndef SZ3_SCALAR_QUANTIZER_HPP
#define SZ3_SCALAR_QUANTIZER_HPP

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * @brief Generic scalar quantizer with one simple reconstruction rule.
 */
template <class Ti, class To = int64_t>
class ScalarQuantizer : public concepts::QuantizerInterface<Ti, To> {
   public:
    /**
     * @brief Reusable scalar quantizer for residual-like values.
     *
     * Designed for decomposition modules that need a small, predictable mapping
     * between floating values and integer bins with configurable reconstruction.
     */
    ScalarQuantizer(double step = 1.0, double one_bin_reconstruct = 1.0, double tail_offset = 0.0)
        : step_(step), one_bin_reconstruct_(one_bin_reconstruct), tail_offset_(tail_offset) {
        static_assert(std::is_arithmetic<Ti>::value, "ScalarQuantizer requires arithmetic Ti.");
        static_assert(std::is_integral<To>::value, "ScalarQuantizer requires integral To.");
        validate_positive_step(step_);
        inv_step_ = 1.0 / step_;
    }

    ALWAYS_INLINE To quantize_and_overwrite(Ti &data, Ti pred) override {
        const long long q =
            std::llrint((static_cast<double>(data) - static_cast<double>(pred)) * inv_step_);
        data = static_cast<Ti>(static_cast<double>(pred) + dequantize_delta(q));
        return static_cast<To>(q);
    }

    ALWAYS_INLINE Ti recover(Ti pred, To quant_index) override {
        const long long q = static_cast<long long>(quant_index);
        return static_cast<Ti>(static_cast<double>(pred) + dequantize_delta(q));
    }

    To force_save_unpred(Ti ori) override {
        return static_cast<To>(std::llrint(static_cast<double>(ori) * inv_step_));
    }

    void save(uchar *&c) const override {
        write(uid(), c);
        write(step_, c);
        write(one_bin_reconstruct_, c);
        write(tail_offset_, c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        uchar uid_read = 0;
        read(uid_read, c, remaining_length);
        if (uid_read != uid()) {
            throw std::invalid_argument("ScalarQuantizer uid mismatch");
        }
        read(step_, c, remaining_length);
        read(one_bin_reconstruct_, c, remaining_length);
        read(tail_offset_, c, remaining_length);
        validate_positive_step(step_);
        inv_step_ = 1.0 / step_;
    }

    std::pair<To, To> get_out_range() const override {
        return std::make_pair(std::numeric_limits<To>::lowest(), std::numeric_limits<To>::max());
    }

    void print() override {
        printf("[ScalarQuantizer] step=%.8G, one_bin_reconstruct=%.8G, tail_offset=%.8G\n", step_,
               one_bin_reconstruct_, tail_offset_);
    }

   private:
    static void validate_positive_step(double step) {
        if (!(step > 0.0)) {
            throw std::invalid_argument("ScalarQuantizer requires positive step.");
        }
    }

    double dequantize_delta(long long q) const {
        if (q == 0) {
            return 0.0;
        }

        const double sign = q < 0 ? -1.0 : 1.0;
        const uint64_t magnitude = static_cast<uint64_t>(std::llabs(q));
        const double reconstructed_magnitude =
            (magnitude == 1) ? one_bin_reconstruct_ : (static_cast<double>(magnitude) + tail_offset_);
        return sign * reconstructed_magnitude * step_;
    }

    static uchar uid() { return 0b100; }

    double step_ = 1.0;
    double inv_step_ = 1.0;
    double one_bin_reconstruct_ = 1.0;
    double tail_offset_ = 0.0;
};

}  // namespace SZ3

#endif
