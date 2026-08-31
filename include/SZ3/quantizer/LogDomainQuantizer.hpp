/**
 * @file LogDomainQuantizer.hpp
 * @ingroup Quantizer
 */

#ifndef SZ3_LOG_DOMAIN_QUANTIZER_HPP
#define SZ3_LOG_DOMAIN_QUANTIZER_HPP

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * @brief Quantizer that works on `log(|x|)`, yielding a pointwise *relative* error bound.
 *
 * Quantizes `ln|x|` on a uniform grid of half-width `h` instead of quantizing the data domain.
 * Since `|x_hat - x| / |x| = |exp(L_hat - L) - 1| <= exp(h) - 1` for `|L_hat - L| <= h`, choosing
 * `h = log1p(rel_eb)` bounds the relative error by `rel_eb`. The reconstruction magnitudes are
 * the geometric ladder `mag(m) = min_abs * exp(log_step * m)`, `log_step = 2 * h`,
 * `m = 0 .. radius-1`.
 *
 * `pred` is ignored by both `quantize_and_overwrite()` and `recover()`: this encodes the value
 * itself, not a residual, so it is not a drop-in replacement for `LinearQuantizer` inside a
 * residual-based decomposition.
 *
 * @par Bin layout
 * - `0` -- unpredictable; the original value is stored verbatim.
 * - `1` -- exact zero.
 * - `2 .. radius+1` -- positive values, magnitude level `bin - 2`.
 * - `radius+2 .. 2*radius+1` -- negative values, magnitude level `bin - radius - 2`.
 *
 * @par Error guarantee
 * - `x == 0` (including `-0`): `x_hat == 0` exactly, `-0` normalized to `+0`.
 * - unpredictable path: `x_hat == x` bit-identically.
 * - otherwise: `|x_hat - x| <= rel_eb * |x|` and `sign(x_hat) == sign(x)`.
 *
 * The achieved relative error is re-checked in double precision against the value actually
 * written back; anything that fails is stored verbatim. That covers non-finite values, magnitudes
 * below `min_abs / (1 + rel_eb)` or above the top rung, and values lost to rounding into `Ti`.
 *
 * @par Guard band
 * The grid is built for `rel_eb - guard`, `guard = 4 * epsilon<Ti>`, because rounding the
 * reconstruction into `Ti` alone costs up to `epsilon/2` of relative error. This is also why
 * `rel_eb` must exceed `4 * epsilon`.
 *
 * @par Choosing `min_abs` and `radius`
 * Together they fix the representable magnitude window, spanning a factor of
 * `exp(log_step * (radius - 1))`; data outside it falls back to verbatim storage. See
 * `calibrate()` for anchoring the window's top at a known data maximum. `radius` is clamped at
 * construction so the top rung stays finite in `Ti`.
 *
 * @par Reproducibility
 * Magnitudes are recomputed with `std::exp` rather than serialized, so bit-identical
 * reconstruction assumes both hosts agree on `std::exp` to the last ulp.
 *
 * @tparam Ti Data type
 */
template <class Ti>
class LogDomainQuantizer : public concepts::QuantizerInterface<Ti, int> {
   public:
    /**
     * @brief Construct a new log-domain quantizer.
     *
     * @param rel_eb Pointwise relative error bound. Must be > `4 * epsilon<Ti>` and < 1.
     * @param min_abs Smallest representable magnitude (the ladder's bottom rung). Must be > 0.
     * @param r Number of magnitude levels per sign. Must be >= 1; clamped down if the requested
     *          ladder would overflow `Ti`.
     * @throw std::invalid_argument if any argument is out of range.
     */
    explicit LogDomainQuantizer(double rel_eb, double min_abs = 1e-20, int r = 32768)
        : rel_eb_(rel_eb), min_abs_(min_abs), radius_(r) {
        static_assert(std::is_floating_point<Ti>::value, "LogDomainQuantizer requires a floating-point type.");
        validate(rel_eb_, min_abs_, radius_);
        configure();
        clamp_radius();
    }

    /**
     * @brief Default-construct with `rel_eb = 1e-3`, `min_abs = 1e-20`, `radius = 32768`.
     */
    LogDomainQuantizer() : LogDomainQuantizer(1e-3) {}

    /**
     * @brief Get the pointwise relative error bound.
     * @return double The bound actually guaranteed to callers.
     */
    double get_eb() const { return rel_eb_; }

    /**
     * @brief Set a new pointwise relative error bound and rebuild the ladder.
     * @param rel_eb New relative bound. Must be > `4 * epsilon<Ti>` and < 1.
     * @throw std::invalid_argument if `rel_eb` is out of range.
     */
    void set_eb(double rel_eb) {
        validate(rel_eb, min_abs_, radius_);
        rel_eb_ = rel_eb;
        configure();
        clamp_radius();
    }

    /**
     * @brief Set the ladder's bottom rung.
     * @param min_abs Smallest representable magnitude. Must be > 0.
     * @throw std::invalid_argument if `min_abs <= 0`.
     */
    void set_min_abs(double min_abs) {
        validate(rel_eb_, min_abs, radius_);
        min_abs_ = min_abs;
        configure();
        clamp_radius();
    }

    /**
     * @brief Anchor the representable window's top at a known data maximum.
     *
     * Sets `min_abs` so the top rung sits at or just above `max_abs`, keeping `radius` fixed. This
     * pushes the window as low as the level budget allows, minimizing how many small magnitudes
     * spill onto the unpredictable path.
     *
     * @param max_abs Largest magnitude present in the data. Must be > 0 and finite.
     * @throw std::invalid_argument if `max_abs` is not a positive finite number.
     */
    void calibrate(double max_abs) {
        if (!(max_abs > 0) || !std::isfinite(max_abs)) {
            throw std::invalid_argument("LogDomainQuantizer::calibrate requires a positive finite max_abs.");
        }
        set_min_abs(max_abs * std::exp(-log_step_ * (radius_ - 1)));
    }

    /**
     * @brief Get the number of magnitude levels per sign (possibly clamped from the requested value).
     * @return int The effective radius.
     */
    int get_radius() const { return radius_; }

    /**
     * @brief Get the ladder's bottom rung.
     * @return double The smallest reconstruction magnitude.
     */
    double min_magnitude() const { return min_abs_; }

    /**
     * @brief Get the ladder's top rung.
     * @return double The largest reconstruction magnitude.
     */
    double max_magnitude() const { return magnitude(radius_ - 1); }

    /**
     * @brief Uniform grid spacing in the log domain.
     * @return double `2 * log1p(rel_eb - guard)`.
     */
    double log_step() const { return log_step_; }

    /**
     * @brief Bin range produced by this quantizer.
     *
     * Bin `0` is the unpredictable marker and bin `1` is the exact zero, so the alphabet is
     * `2 * radius + 2` symbols. The upper bound is exclusive.
     *
     * @return std::pair<int, int> `{0, 2 * radius + 2}`.
     */
    std::pair<int, int> get_out_range() const override { return std::make_pair(0, 2 * radius_ + 2); }

    /**
     * @brief Quantize `data` in the log domain and overwrite it with the reconstruction.
     *
     * @param data Data point (in/out). Overwritten with the reconstructed value, or left untouched
     *             (and pushed to the unpredictable buffer) when the relative bound cannot be met.
     * @param pred Ignored -- see the class documentation.
     * @return int Bin index in `[0, 2 * radius + 1]`; `0` means unpredictable, `1` means exact zero.
     */
    ALWAYS_INLINE int quantize_and_overwrite(Ti &data, Ti pred) override {
        (void)pred;
        const double value = static_cast<double>(data);
        if (value == 0.0) {
            data = static_cast<Ti>(0);  // normalizes -0 to +0
            return kZeroBin;
        }
        if (std::isfinite(value)) {
            const double abs_value = std::fabs(value);
            const double level_real = (std::log(abs_value) - log_min_) / log_step_;
            // A false comparison also catches NaN. The window is [-0.5, radius-0.5] after rounding;
            // testing against [-1, radius] keeps the llrint argument safely in range.
            if (level_real >= -1.0 && level_real <= static_cast<double>(radius_)) {
                const auto level = static_cast<int>(std::llrint(level_real));
                if (level >= 0 && level < radius_) {
                    const bool negative = value < 0.0;
                    const Ti reconstructed = reconstruct(level, negative);
                    // Verify on the value we would actually write back, in double precision.
                    const double err = std::fabs(static_cast<double>(reconstructed) - value) / abs_value;
                    if (err <= rel_eb_) {
                        data = reconstructed;
                        return bin_of(level, negative);
                    }
                }
            }
        }
        unpred_.push_back(data);
        return kUnpredBin;
    }

    /**
     * @brief Reconstruct a data point from its bin index.
     *
     * @param pred Ignored -- see the class documentation.
     * @param quant_index Bin index as returned by `quantize_and_overwrite()`.
     * @return Ti The reconstructed value.
     */
    ALWAYS_INLINE Ti recover(Ti pred, int quant_index) override {
        (void)pred;
        if (quant_index == kUnpredBin) {
            return recover_unpred();
        }
        if (quant_index == kZeroBin) {
            return static_cast<Ti>(0);
        }
        const int offset = quant_index - kFirstPositiveBin;
        const bool negative = offset >= radius_;
        return reconstruct(negative ? offset - radius_ : offset, negative);
    }

    /**
     * @brief Pop the next verbatim value from the unpredictable buffer.
     * @return Ti The stored original value.
     */
    ALWAYS_INLINE Ti recover_unpred() {
        if (index_ >= unpred_.size()) {
            throw std::runtime_error("LogDomainQuantizer: more unpredictable bins than stored values");
        }
        return unpred_[index_++];
    }

    /**
     * @brief Force a value onto the unpredictable path.
     * @param ori Original value, stored verbatim.
     * @return int Always `0`.
     */
    ALWAYS_INLINE int force_save_unpred(Ti ori) override {
        unpred_.push_back(ori);
        return kUnpredBin;
    }

    /**
     * @brief Estimated serialized size of the unpredictable buffer, in bytes.
     * @return size_t Byte count.
     */
    size_t size_est() const { return unpred_.size() * sizeof(Ti); }

    /**
     * @brief Number of values currently on the unpredictable path.
     * @return size_t Element count.
     */
    size_t unpred_count() const { return unpred_.size(); }

    /**
     * @brief Serialize the quantizer state.
     * @param c Buffer pointer, advanced past the written bytes.
     */
    void save(uchar *&c) const override {
        write(uid_, c);
        write(rel_eb_, c);
        write(min_abs_, c);
        write(radius_, c);
        size_t unpred_size = unpred_.size();
        write(unpred_size, c);
        if (unpred_size > 0) {
            write(unpred_.data(), unpred_.size(), c);
        }
    }

    /**
     * @brief Deserialize the quantizer state and rebuild the derived ladder parameters.
     *
     * The stored `radius` is the already-clamped effective value and is used as-is.
     *
     * @param c Buffer pointer, advanced past the read bytes.
     * @param remaining_length Remaining buffer length, decremented as bytes are consumed.
     * @throw std::invalid_argument on uid mismatch or an invalid stored configuration.
     */
    void load(const uchar *&c, size_t &remaining_length) override {
        uchar uid_read = 0;
        read(uid_read, c, remaining_length);
        if (uid_read != uid_) {
            throw std::invalid_argument("LogDomainQuantizer uid mismatch");
        }
        read(rel_eb_, c, remaining_length);
        read(min_abs_, c, remaining_length);
        read(radius_, c, remaining_length);
        validate(rel_eb_, min_abs_, radius_);
        configure();
        size_t unpred_size = 0;
        read(unpred_size, c, remaining_length);
        unpred_.clear();
        if (unpred_size > 0) {
            unpred_.resize(unpred_size);
            read(unpred_.data(), unpred_size, c, remaining_length);
        }
        index_ = 0;
    }

    /**
     * @brief Print a one-line summary of the quantizer configuration.
     */
    void print() override {
        printf("[LogDomainQuantizer] rel_eb = %.8G, window = [%.8G, %.8G], radius = %d, unpred = %zu\n", rel_eb_,
               min_abs_, max_magnitude(), radius_, unpred_.size());
    }

   private:
    static constexpr int kUnpredBin = 0;
    static constexpr int kZeroBin = 1;
    static constexpr int kFirstPositiveBin = 2;

    /// Rounding the reconstruction into `Ti` costs up to epsilon/2 of relative error; reserve 4x.
    static constexpr double guard_band() { return 4.0 * static_cast<double>(std::numeric_limits<Ti>::epsilon()); }

    static void validate(double rel_eb, double min_abs, int r) {
        if (!(rel_eb > guard_band())) {
            throw std::invalid_argument(
                "LogDomainQuantizer: rel_eb must exceed 4 * epsilon of the data type; a relative bound finer than "
                "the storage type's own resolution cannot be honored.");
        }
        if (!(rel_eb < 1.0)) {
            throw std::invalid_argument("LogDomainQuantizer: rel_eb must be < 1.");
        }
        if (!(min_abs > 0) || !std::isfinite(min_abs)) {
            throw std::invalid_argument("LogDomainQuantizer: min_abs must be a positive finite number.");
        }
        if (r < 1) {
            throw std::invalid_argument("LogDomainQuantizer: radius must be >= 1.");
        }
    }

    /// Recompute the derived log-domain grid parameters from `rel_eb_` and `min_abs_`.
    void configure() {
        log_step_ = 2.0 * std::log1p(rel_eb_ - guard_band());
        log_min_ = std::log(min_abs_);
    }

    /// Shrink `radius_` if the requested ladder's top rung would overflow `Ti`.
    void clamp_radius() {
        const double max_log = std::log(static_cast<double>(std::numeric_limits<Ti>::max()));
        const auto affordable = static_cast<long long>(std::floor((max_log - log_min_) / log_step_)) + 1;
        if (affordable < 1) {
            throw std::invalid_argument("LogDomainQuantizer: min_abs is too large to represent in this data type.");
        }
        if (static_cast<long long>(radius_) > affordable) {
            radius_ = static_cast<int>(affordable);
        }
    }

    double magnitude(int level) const { return min_abs_ * std::exp(log_step_ * level); }

    /// Single source of truth for reconstruction, shared by compression and decompression.
    ALWAYS_INLINE Ti reconstruct(int level, bool negative) const {
        const double mag = magnitude(level);
        return static_cast<Ti>(negative ? -mag : mag);
    }

    /// Inverse of the offset decoding in `recover()`: positives occupy the first `radius_` bins.
    ALWAYS_INLINE int bin_of(int level, bool negative) const {
        return kFirstPositiveBin + level + (negative ? radius_ : 0);
    }

    std::vector<Ti> unpred_;
    size_t index_ = 0;  // used in decompression only
    static constexpr uchar uid_ = 0b111;

    double rel_eb_ = 1e-3;
    double min_abs_ = 1e-20;
    int radius_ = 32768;
    double log_step_ = 0.0;  // derived: uniform grid spacing in the log domain
    double log_min_ = 0.0;   // derived: log(min_abs_)
};

}  // namespace SZ3
#endif
