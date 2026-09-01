/**
 * @file OutlierQuantizer.hpp
 * @ingroup Quantizer
 */

#ifndef SZ3_OUTLIER_QUANTIZER_HPP
#define SZ3_OUTLIER_QUANTIZER_HPP

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/quantizer/ScalarQuantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * @brief Sparse "second pass" quantizer for values a lossy stage failed to represent.
 *
 * Plays the role `LinearQuantizer::unpred` does: values whose reconstruction error exceeds the
 * tolerance are recorded as a quantized *correction* in the module's own serialized state rather
 * than in the integer stream handed to the encoder.
 *
 * Per element, given an approximation `pred` of a true value `data`:
 *  - `|data - pred| <= tolerance` -> bin 0, nothing recorded, `data` snaps to `pred`.
 *  - otherwise -> the residual is quantized with an inner
 *    `ScalarQuantizer(tolerance, 1.1, -0.25)` and a non-zero bin is returned.
 *
 * `|residual| > tolerance` implies `|residual| / tolerance > 1`, so the inner quantizer never
 * rounds an outlier to bin 0 and "bin == 0" unambiguously means "no correction". The `1.1` /
 * `-0.25` tweaks bias a one-bin residual to `1.1 * tolerance` and an `m`-bin residual to
 * `(m - 0.25) * tolerance`, keeping the corrected value inside the tolerance band.
 *
 * @tparam Ti Value type being corrected (floating point)
 * @tparam To Integer bin type
 */
template <class Ti, class To = int64_t>
class OutlierQuantizer : public concepts::QuantizerInterface<Ti, To> {
   public:
    /// One recorded correction: position in the field, and its quantized residual bin.
    struct Correction {
        size_t pos = 0;
        To bin = To{0};
    };

    /**
     * @brief Construct with the tolerance that defines an outlier.
     * @param tolerance Absolute error bound; residuals with `|r| > tolerance` become corrections
     */
    explicit OutlierQuantizer(double tolerance = 1.0)
        : tolerance_(tolerance), quantizer_(tolerance, kOneBinReconstruct, kTailOffset) {
        static_assert(std::is_arithmetic<Ti>::value, "OutlierQuantizer requires arithmetic Ti.");
        static_assert(std::is_integral<To>::value, "OutlierQuantizer requires integral To.");
        validate_positive_tolerance(tolerance_);
    }

    /// Tolerance currently in effect.
    double get_tolerance() const { return tolerance_; }

    /**
     * @brief Quantize one residual, snapping `data` to the value the decoder will rebuild.
     * @param data Value to be corrected; overwritten with the reconstruction
     * @param pred Approximation produced by the upstream lossy stage
     * @return Correction bin, or 0 when the approximation is already within tolerance
     */
    ALWAYS_INLINE To quantize_and_overwrite(Ti &data, Ti pred) override {
        const double residual = static_cast<double>(data) - static_cast<double>(pred);
        if (!(std::abs(residual) > tolerance_)) {
            data = pred;
            return To{0};
        }

        Ti scratch = static_cast<Ti>(residual);
        const To bin = quantizer_.quantize_and_overwrite(scratch, static_cast<Ti>(0));
        data = static_cast<Ti>(static_cast<double>(pred) + static_cast<double>(scratch));
        return bin;
    }

    /**
     * @brief Rebuild a value from an approximation and a correction bin.
     */
    ALWAYS_INLINE Ti recover(Ti pred, To quant_index) override {
        if (quant_index == To{0}) {
            return pred;
        }
        return static_cast<Ti>(static_cast<double>(pred) +
                               static_cast<double>(quantizer_.recover(static_cast<Ti>(0), quant_index)));
    }

    To force_save_unpred(Ti ori) override { return quantizer_.force_save_unpred(ori); }

    /**
     * @brief Scan a whole field and record every element the approximation misses.
     *
     * Appends to the internal correction list (call `clear()` first to restart).
     *
     * @param original Reference values
     * @param approximation Values reconstructed by the upstream lossy stage
     * @return Number of corrections recorded by this call
     */
    size_t collect(const std::vector<Ti> &original, const std::vector<Ti> &approximation) {
        if (original.size() != approximation.size()) {
            throw std::runtime_error("OutlierQuantizer length mismatch between original and approximation.");
        }

        const size_t before = corrections_.size();
        corrections_.reserve(corrections_.size() + original.size() / 20);
        for (size_t i = 0; i < original.size(); i++) {
            Ti value = original[i];
            const To bin = quantize_and_overwrite(value, approximation[i]);
            if (bin != To{0}) {
                corrections_.push_back(Correction{i, bin});
            }
        }
        return corrections_.size() - before;
    }

    /**
     * @brief Apply every recorded correction in place.
     * @param values Field in the same domain `collect()` ran on
     */
    void apply(std::vector<Ti> &values) {
        for (size_t i = 0; i < corrections_.size(); i++) {
            const Correction &correction = corrections_[i];
            if (correction.pos >= values.size()) {
                throw std::runtime_error("OutlierQuantizer correction index out of bounds.");
            }
            values[correction.pos] = recover(values[correction.pos], correction.bin);
        }
    }

    /// Recorded corrections, in ascending position order when produced by `collect()`.
    const std::vector<Correction> &get_corrections() const { return corrections_; }

    /// Replace the recorded corrections wholesale.
    void set_corrections(std::vector<Correction> corrections) { corrections_ = std::move(corrections); }

    /// Number of recorded corrections.
    size_t size() const { return corrections_.size(); }

    /// Drop all recorded corrections (the tolerance is kept).
    void clear() { corrections_.clear(); }

    /**
     * @brief Expand the sparse correction list into a dense, mostly-zero bin array.
     *
     * Provided for interoperability with pipelines (such as `SPERRFusedDecomposition`) that
     * ship outlier corrections as a full-length integer stream.
     */
    std::vector<To> to_dense(size_t len) const {
        std::vector<To> dense(len, To{0});
        for (size_t i = 0; i < corrections_.size(); i++) {
            if (corrections_[i].pos >= len) {
                throw std::runtime_error("OutlierQuantizer correction index out of bounds.");
            }
            dense[corrections_[i].pos] = corrections_[i].bin;
        }
        return dense;
    }

    /// Inverse of `to_dense()`: adopt the non-zero entries of a dense bin array.
    void from_dense(const std::vector<To> &dense) {
        corrections_.clear();
        for (size_t i = 0; i < dense.size(); i++) {
            if (dense[i] != To{0}) {
                corrections_.push_back(Correction{i, dense[i]});
            }
        }
    }

    /**
     * @brief Serialize tolerance and correction list.
     * @param c Buffer pointer; advanced past the written bytes.
     */
    void save(uchar *&c) const override {
        write(uid(), c);
        write(tolerance_, c);
        const size_t count = corrections_.size();
        write(count, c);
        for (size_t i = 0; i < count; i++) {
            write(corrections_[i].pos, c);
            write(corrections_[i].bin, c);
        }
    }

    /**
     * @brief Deserialize tolerance and correction list.
     * @param c Buffer pointer; advanced past the consumed bytes.
     * @param remaining_length Remaining buffer length; decremented.
     */
    void load(const uchar *&c, size_t &remaining_length) override {
        uchar uid_read = 0;
        read(uid_read, c, remaining_length);
        if (uid_read != uid()) {
            throw std::invalid_argument("OutlierQuantizer uid mismatch");
        }
        read(tolerance_, c, remaining_length);
        validate_positive_tolerance(tolerance_);
        quantizer_ = ScalarQuantizer<Ti, To>(tolerance_, kOneBinReconstruct, kTailOffset);

        size_t count = 0;
        read(count, c, remaining_length);
        if (count * (sizeof(size_t) + sizeof(To)) > remaining_length) {
            throw std::runtime_error("OutlierQuantizer correction list exceeds the remaining buffer.");
        }
        corrections_.resize(count);
        for (size_t i = 0; i < count; i++) {
            read(corrections_[i].pos, c, remaining_length);
            read(corrections_[i].bin, c, remaining_length);
        }
    }

    /// Serialized size of `save()`.
    size_t size_est() const {
        return sizeof(uchar) + sizeof(double) + sizeof(size_t) + corrections_.size() * (sizeof(size_t) + sizeof(To));
    }

    std::pair<To, To> get_out_range() const override {
        return std::make_pair(std::numeric_limits<To>::lowest(), std::numeric_limits<To>::max());
    }

    void print() override {
        printf("[OutlierQuantizer] tolerance=%.8G, corrections=%zu\n", tolerance_, corrections_.size());
    }

   private:
    static void validate_positive_tolerance(double tolerance) {
        if (!(tolerance > 0.0)) {
            throw std::invalid_argument("OutlierQuantizer requires a positive tolerance.");
        }
    }

    static uchar uid() { return 0b1011; }

    /// SPERR's reconstruction tweaks for outlier corrections.
    static constexpr double kOneBinReconstruct = 1.1;
    static constexpr double kTailOffset = -0.25;

    double tolerance_ = 1.0;
    ScalarQuantizer<Ti, To> quantizer_;
    std::vector<Correction> corrections_;
};

}  // namespace SZ3

#endif  // SZ3_OUTLIER_QUANTIZER_HPP
