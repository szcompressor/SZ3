/**
 * @file LevelQuantizer.hpp
 * @ingroup Quantizer
 */

#ifndef SZ3_LEVEL_QUANTIZER_HPP
#define SZ3_LEVEL_QUANTIZER_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * @brief Shape of the non-uniform reconstruction-level curve used by `LevelQuantizer`.
 *
 * Selected at construction (runtime, not a template parameter) and written into the serialized
 * state, so a stream is self-describing and `load()` can reject a mismatched configuration instead
 * of silently rebuilding the wrong table. The underlying values are part of the compressed format:
 * append new curves, never renumber existing ones.
 */
enum class LevelCurve : uchar {
    /// `pos[i] = eb + 2 * eb * (radius - 1) * (i / (radius - 1))^2` -- spacing grows quadratically.
    Quadratic = 0,
    /// `pos[i] = eb * (2 * radius - 1)^(i / (radius - 1))` -- geometric, i.e. log-linear in `i`.
    Log = 1,
};

/**
 * @brief Lookup-table quantizer with non-uniformly spaced reconstruction levels.
 *
 * Precomputes `2 * radius` reconstruction levels for the prediction residual, symmetric around
 * zero, and picks the nearest by binary search. `LevelCurve` selects the spacing at construction:
 *
 *  - `LevelCurve::Quadratic` -- `pos[i] = eb + 2 * eb * (radius - 1) * (i / (radius - 1))^2`
 *  - `LevelCurve::Log`       -- `pos[i] = eb * (2 * radius - 1)^(i / (radius - 1))`
 *
 * Both start at `eb` and reach the same nominal top level, differing only in how densely they
 * pack levels near zero. Neither curve bounds the error on its own, so the table is post-clamped
 * so consecutive positive levels are never more than `2 * eb` apart; this turns each curve's tail
 * into a uniform `2 * eb` ladder and lowers the largest representable residual.
 *
 * @par Error guarantee
 * `|reconstructed - original| <= eb` for every element, both curves. Bin `0` is reserved for
 * unpredictable values: whenever the best level would miss the bound -- residual outside the
 * table, `data`/`pred` not finite, or `pred + level` rounding away in `T` -- the original is
 * stored verbatim and reconstructs exactly. The check is made on the value actually written back.
 *
 * @par No zero level
 * The smallest positive level is `eb`, so a residual of exactly `0` reconstructs to `+/-eb`.
 * Perfectly predicted points pay full error. Changing this would change the compressed format.
 *
 * @par Reproducibility
 * The table is recomputed from `(curve, eb, radius)` on load rather than serialized, and the log
 * curve uses `std::exp`, so bit-identical reconstruction assumes both hosts agree on `std::exp`
 * to the last ulp. A disagreement shifts a level by at most one ulp. The quadratic curve uses
 * only arithmetic and is unaffected.
 *
 * @tparam T Data type
 */
template <class T>
class LevelQuantizer : public concepts::QuantizerInterface<T, int> {
   public:
    /**
     * @brief Construct a new quantizer.
     *
     * @param eb Absolute error bound. Must be > 0. The smallest positive level equals `eb`.
     * @param r Quantization radius; the table holds `2 * r` levels. Must be >= 2.
     * @param curve Shape of the level curve. Defaults to `LevelCurve::Quadratic`.
     * @throw std::invalid_argument if `eb <= 0`, `r < 2`, or `curve` is not a known enumerator.
     */
    LevelQuantizer(double eb, int r = 32768, LevelCurve curve = LevelCurve::Quadratic)
        : error_bound(eb), radius(r), curve_(curve) {
        validate(error_bound, radius, curve_);
        build_quant_table();
    }

    /**
     * @brief Default-construct with `eb = 1e-4`, `radius = 32768`, `LevelCurve::Quadratic`.
     */
    LevelQuantizer() { build_quant_table(); }

    /**
     * @brief Get the absolute error bound.
     * @return double The current error bound.
     */
    double get_eb() const { return error_bound; }

    /**
     * @brief Set a new absolute error bound and rebuild the level table.
     * @param eb New absolute error bound. Must be > 0.
     * @throw std::invalid_argument if `eb <= 0`.
     */
    void set_eb(double eb) {
        validate(eb, radius, curve_);
        error_bound = eb;
        build_quant_table();
    }

    /**
     * @brief Get the quantization radius.
     * @return int The number of positive levels (the table holds `2 * radius`).
     */
    int get_radius() const { return radius; }

    /**
     * @brief Get the level curve this quantizer was configured with.
     * @return LevelCurve The active curve.
     */
    LevelCurve get_curve() const { return curve_; }

    /**
     * @brief Largest residual magnitude the table can represent after post-clamping.
     * @return double The top positive reconstruction level.
     */
    double max_level() const { return static_cast<double>(quant_table.back()); }

    /**
     * @brief Bin range produced by this quantizer.
     *
     * Bin `0` is the unpredictable marker; bins `1 .. 2 * radius` index the level table. The upper
     * bound is exclusive, matching the alphabet size expected by the encoders.
     *
     * @return std::pair<int, int> `{0, 2 * radius + 1}`.
     */
    std::pair<int, int> get_out_range() const override { return std::make_pair(0, 2 * radius + 1); }

    /**
     * @brief Build the post-clamped level table for the active curve.
     *
     * Called by the constructors, `set_eb()` and `load()`. Fills `quant_table` with `2 * radius`
     * non-decreasing values `[-pos[radius-1] .. -pos[0], pos[0] .. pos[radius-1]]`. The sequence is
     * strictly increasing except where `T`'s precision collapses adjacent levels, which the binary
     * search tolerates.
     */
    void build_quant_table() {
        quant_table.resize(static_cast<size_t>(2) * radius);
        std::vector<double> pos_levels(radius);
        pos_levels[0] = error_bound;  // Smallest positive level, identical for both curves.

        const double inv = 1.0 / static_cast<double>(radius - 1);
        if (curve_ == LevelCurve::Log) {
            // Geometric spacing: log(pos_levels[i]) is linear in i.
            const double log_span = std::log(static_cast<double>(2 * radius - 1));
            for (int i = 1; i < radius; i++) {
                pos_levels[i] = error_bound * std::exp(log_span * (static_cast<double>(i) * inv));
            }
        } else {
            // Quadratic spacing. At i = radius-1 this reaches eb * (2 * radius - 1), the same
            // nominal top level as the log curve.
            const double span = 2 * error_bound * (radius - 1);
            for (int i = 1; i < radius; i++) {
                const double r = static_cast<double>(i) * inv;
                pos_levels[i] = error_bound + span * r * r;
            }
        }

        // Post-clamp: consecutive levels never farther apart than 2 * eb, which is what preserves
        // the absolute error bound. Forward pass, so the clamp is cumulative.
        const double max_gap = 2 * error_bound;
        for (int i = 1; i < radius; i++) {
            if (pos_levels[i] - pos_levels[i - 1] > max_gap) {
                pos_levels[i] = pos_levels[i - 1] + max_gap;
            }
        }

        for (int i = 0; i < radius; i++) {
            quant_table[radius + i] = static_cast<T>(pos_levels[i]);
            quant_table[radius - 1 - i] = static_cast<T>(-pos_levels[i]);
        }
    }

    /**
     * @brief Quantize the residual `data - pred` and overwrite `data` with the reconstruction.
     *
     * @param data Data point (in/out). Overwritten with the reconstructed value, or left untouched
     *             (and pushed to the unpredictable buffer) when the bound cannot be met.
     * @param pred Predicted value.
     * @return int Bin index in `[0, 2 * radius]`; `0` means unpredictable.
     */
    ALWAYS_INLINE int quantize_and_overwrite(T &data, T pred) override {
        const T diff = data - pred;
        // Rejects NaN/Inf data, NaN/Inf pred, and residual overflow. Also keeps NaN out of the
        // binary search, whose comparator would otherwise not be a strict weak ordering.
        if (std::isfinite(static_cast<double>(diff))) {
            const int level = search_best_level(diff);
            const T reconstructed = pred + quant_table[level];
            // Verify on the value we would actually write back, not on the idealized level.
            if (std::fabs(static_cast<double>(reconstructed) - static_cast<double>(data)) <= error_bound) {
                data = reconstructed;
                return level + 1;
            }
        }
        unpred.push_back(data);
        return 0;
    }

    /**
     * @brief Reconstruct a data point from its bin index.
     *
     * @param pred Predicted value (ignored for bin `0`).
     * @param quant_index Bin index as returned by `quantize_and_overwrite()`.
     * @return T The reconstructed value.
     */
    ALWAYS_INLINE T recover(T pred, int quant_index) override {
        if (quant_index) {
            return recover_pred(pred, quant_index);
        }
        return recover_unpred();
    }

    /**
     * @brief Reconstruct a predictable value. Bin index must be non-zero.
     * @param pred Predicted value.
     * @param quant_index Bin index in `[1, 2 * radius]`.
     * @return T The reconstructed value.
     */
    ALWAYS_INLINE T recover_pred(T pred, int quant_index) { return pred + quant_table[quant_index - 1]; }

    /**
     * @brief Pop the next verbatim value from the unpredictable buffer.
     * @return T The stored original value.
     */
    ALWAYS_INLINE T recover_unpred() {
        if (index >= unpred.size()) {
            throw std::runtime_error("LevelQuantizer: more unpredictable bins than stored values");
        }
        return unpred[index++];
    }

    /**
     * @brief Force a value onto the unpredictable path.
     * @param ori Original value, stored verbatim.
     * @return int Always `0`.
     */
    ALWAYS_INLINE int force_save_unpred(T ori) override {
        unpred.push_back(ori);
        return 0;
    }

    /**
     * @brief Estimated serialized size of the unpredictable buffer, in bytes.
     * @return size_t Byte count.
     */
    size_t size_est() const { return unpred.size() * sizeof(T); }

    /**
     * @brief Number of values currently on the unpredictable path.
     * @return size_t Element count.
     */
    size_t unpred_count() const { return unpred.size(); }

    /**
     * @brief Serialize the quantizer state (uid, curve, eb, radius, unpredictable buffer).
     * @param c Buffer pointer, advanced past the written bytes.
     */
    void save(unsigned char *&c) const override {
        write(uid_, c);
        write(static_cast<uchar>(curve_), c);
        write(this->error_bound, c);
        write(this->radius, c);
        size_t unpred_size = unpred.size();
        write(unpred_size, c);
        if (unpred_size > 0) {
            write(unpred.data(), unpred.size(), c);
        }
    }

    /**
     * @brief Deserialize the quantizer state and rebuild the level table.
     *
     * `eb` and `radius` are taken from the stream, but the curve is *validated* against the curve
     * this instance was constructed with rather than adopted silently: the curve determines the
     * meaning of every bin, so a stream written with the other curve is a configuration error, not
     * something to paper over. Construct the reader with the matching `LevelCurve`.
     *
     * @param c Buffer pointer, advanced past the read bytes.
     * @param remaining_length Remaining buffer length, decremented as bytes are consumed.
     * @throw std::invalid_argument on uid mismatch, an unknown or mismatched curve, or an invalid
     *        stored configuration.
     */
    void load(const unsigned char *&c, size_t &remaining_length) override {
        uchar uid_read;
        read(uid_read, c, remaining_length);
        if (uid_read != uid_) {
            throw std::invalid_argument("LevelQuantizer uid mismatch");
        }
        uchar curve_read;
        read(curve_read, c, remaining_length);
        if (!is_known_curve(curve_read)) {
            throw std::invalid_argument("LevelQuantizer: unknown level curve in stream.");
        }
        if (static_cast<LevelCurve>(curve_read) != curve_) {
            throw std::invalid_argument("LevelQuantizer level curve mismatch: stream and quantizer disagree.");
        }
        read(this->error_bound, c, remaining_length);
        read(this->radius, c, remaining_length);
        validate(this->error_bound, this->radius, curve_);
        size_t unpred_size = 0;
        read(unpred_size, c, remaining_length);
        unpred.clear();
        if (unpred_size > 0) {
            unpred.resize(unpred_size);
            read(unpred.data(), unpred_size, c, remaining_length);
        }
        index = 0;
        build_quant_table();
    }

    /**
     * @brief Print a one-line summary of the quantizer configuration.
     */
    void print() override {
        printf("[LevelQuantizer] curve = %s, error_bound = %.8G, radius = %d, max_level = %.8G, unpred = %zu\n",
               curve_name(curve_), error_bound, radius, max_level(), unpred.size());
    }

    /**
     * @brief Human-readable name of a level curve, for diagnostics.
     * @param curve The curve to name.
     * @return const char* A static string, `"?"` for an unknown enumerator.
     */
    static const char *curve_name(LevelCurve curve) {
        switch (curve) {
            case LevelCurve::Quadratic:
                return "quadratic";
            case LevelCurve::Log:
                return "log";
        }
        return "?";
    }

   private:
    /// True if `v` is one of the `LevelCurve` enumerators this build knows about.
    static bool is_known_curve(uchar v) {
        return v == static_cast<uchar>(LevelCurve::Quadratic) || v == static_cast<uchar>(LevelCurve::Log);
    }

    static void validate(double eb, int r, LevelCurve curve) {
        if (!(eb > 0)) {
            throw std::invalid_argument("LevelQuantizer requires a positive error bound.");
        }
        if (r < 2) {
            throw std::invalid_argument("LevelQuantizer requires radius >= 2.");
        }
        if (!is_known_curve(static_cast<uchar>(curve))) {
            throw std::invalid_argument("LevelQuantizer: unknown level curve.");
        }
    }

    /**
     * @brief Index of the table level nearest to `diff`.
     *
     * `quant_table` is sorted, so a binary search is valid and keeps the cost at O(log radius)
     * instead of the O(radius) linear scan.
     */
    int search_best_level(T diff) const {
        const auto begin = quant_table.begin();
        const auto it = std::lower_bound(begin, quant_table.end(), diff);
        if (it == quant_table.end()) {
            return static_cast<int>(quant_table.size()) - 1;
        }
        if (it == begin) {
            return 0;
        }
        const int hi = static_cast<int>(it - begin);
        const int lo = hi - 1;
        return (quant_table[hi] - diff < diff - quant_table[lo]) ? hi : lo;
    }

    std::vector<T> unpred;
    size_t index = 0;  // used in decompression only
    static constexpr uchar uid_ = 0b1010;

    double error_bound = 1e-4;
    int radius = 32768;  // number of positive levels; the table holds 2 * radius entries
    LevelCurve curve_ = LevelCurve::Quadratic;
    std::vector<T> quant_table;  // sorted, symmetric around zero
};

}  // namespace SZ3
#endif
