/**
 * @file ClusterQuantizer.hpp
 * @ingroup Quantizer
 * @brief Cluster-based (codebook) quantizer over a uniform lattice of levels, for data that
 *        clusters around a few discrete values rather than forming a smooth field.
 *
 * Building the codebook needs the data distribution, which `precompress_data()` cannot see,
 * so the two halves are separate: `derive_cluster_levels()` samples the data and returns a
 * `ClusterLevels`; `ClusterQuantizer` is constructed from one. `save()`/`load()` carry the
 * levels, so the decompressor never re-derives them.
 *
 * With `level(l) = level_start + l * level_offset`, an input within `eb` of its nearest level
 * is coded as bin `l+1`; otherwise it takes bin `0` and is stored verbatim. So for every input
 * `|x_hat - x| <= eb`, and bin 0 recovers bit-exactly. `eb` defaults to `level_offset / 2`.
 *
 * `pred` is ignored -- the codebook is absolute. Bin `0` is the unpredictable sentinel, as in
 * `LinearQuantizer`, so `recover()` must be called in encode order.
 */

#ifndef SZ3_CLUSTER_QUANTIZER_HPP
#define SZ3_CLUSTER_QUANTIZER_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/KmeansUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * @brief A uniform lattice codebook: `level(l) = level_start + l * level_offset`, `l in [0, level_num)`.
 *
 * This is exactly the triple `KmeansUtil::get_cluster()` produces. `level_num == 0` is that
 * function's "no usable clustering found" signal and is propagated here as `!valid()`.
 */
struct ClusterLevels {
    float level_start = 0.0f;   ///< Value of level 0.
    float level_offset = 0.0f;  ///< Spacing between consecutive levels (> 0 when valid).
    int level_num = 0;          ///< Number of levels; 0 means "no usable clustering".

    /// @brief True when the triple describes a usable codebook.
    bool valid() const { return level_num > 0 && level_offset > 0.0f && std::isfinite(level_start); }
};

/**
 * @brief Derive a lattice codebook from data using the 1-D k-means in `KmeansUtil.hpp`.
 *
 * Samples `data`, hands the sample to `KmeansUtil::get_cluster()`, and then re-phases the
 * returned lattice so its index range covers the full `[min, max]` of `data` (the lattice
 * *pitch and phase*, which is what k-means actually determines, are preserved; only the index
 * origin and count move). Non-finite entries are ignored.
 *
 * Notable implementation choices, both forced by `KmeansUtil.hpp` being used as-is:
 *  - The sample is materialised here, as a contiguous `std::vector<float>`, and passed to
 *    `get_cluster()` with `num == sample_num` so that its *internal* sampling branch is not
 *    taken. That branch writes past the end of a merely-`reserve()`d vector and can index one
 *    past the input, so it is deliberately bypassed. Sampling here is a deterministic uniform
 *    stride, which also makes compression reproducible run to run.
 *  - `get_cluster()` only instantiates for `float` (its centroid buffer is hard-coded to
 *    `std::vector<float>`), so a `double` input is sampled down to `float`. The returned
 *    levels are `float` regardless, matching `ClusterLevels`.
 *
 * @tparam T Input data type (`float`, `double`, ...); must be convertible to `float`.
 * @param data Pointer to the data to analyse.
 * @param num Number of elements in `data`.
 * @param sample_num Number of samples to cluster. 0 selects the MDZ heuristic
 *                   `clamp(num / 10, min(5000, num), 20000)`.
 * @param max_levels Reject codebooks that would need more than this many levels to span the
 *                   data; such a "clustering" is not a codebook and is reported as invalid.
 * @return ClusterLevels A valid codebook, or a `!valid()` one when the data does not cluster.
 */
template <class T>
ClusterLevels derive_cluster_levels(const T *data, size_t num, size_t sample_num = 0, int max_levels = 65536) {
    ClusterLevels levels;
    if (data == nullptr || num == 0 || max_levels < 1) {
        return levels;
    }

    // ---- full-range min/max over finite entries -------------------------------------
    double min_val = std::numeric_limits<double>::infinity();
    double max_val = -std::numeric_limits<double>::infinity();
    size_t finite_num = 0;
    for (size_t i = 0; i < num; i++) {
        const double v = static_cast<double>(data[i]);
        if (!std::isfinite(v)) continue;
        finite_num++;
        if (v < min_val) min_val = v;
        if (v > max_val) max_val = v;
    }
    if (finite_num < 2 || min_val >= max_val) {
        return levels;  // constant or degenerate data: nothing to cluster
    }

    // ---- deterministic uniform-stride sample ----------------------------------------
    if (sample_num == 0) {
        sample_num = std::max(std::min(static_cast<size_t>(5000), num), std::min(num / 10, static_cast<size_t>(20000)));
    }
    sample_num = std::min(sample_num, num);
    std::vector<float> sample;
    sample.reserve(sample_num);
    for (size_t i = 0; i < sample_num; i++) {
        const size_t idx = (num == sample_num) ? i : (i * num) / sample_num;
        const float v = static_cast<float>(data[idx]);
        if (std::isfinite(v)) {
            sample.push_back(v);
        }
    }
    if (sample.size() < 8) {
        return levels;  // too few points for the k-means elbow heuristic to mean anything
    }

    // ---- optimal 1-D k-means (KmeansUtil.hpp, used as-is) ----------------------------
    float level_start = 0.0f;
    float level_offset = 0.0f;
    int level_num = 0;
    // `get_cluster()` leaves start/offset untouched when it finds no clusters, hence the
    // explicit initialisation above; passing num == sample_num takes its safe branch.
    get_cluster(sample.data(), sample.size(), level_start, level_offset, level_num, sample.size());
    if (level_num <= 0 || !std::isfinite(level_start) || !std::isfinite(level_offset) || level_offset <= 0.0f) {
        return levels;
    }

    // ---- re-phase so index 0..level_num-1 spans the whole data range -----------------
    const double offset = static_cast<double>(level_offset);
    double start = static_cast<double>(level_start);
    const double lo_units = (min_val - start) / offset;
    if (!std::isfinite(lo_units) || std::fabs(lo_units) > static_cast<double>(max_levels)) {
        return levels;
    }
    if (lo_units < 0.0) {
        start += std::floor(lo_units) * offset;
    }
    const double span_units = (max_val - start) / offset;
    if (!std::isfinite(span_units) || span_units < 0.0) {
        return levels;
    }
    const double needed = std::ceil(span_units) + 1.0;
    if (needed > static_cast<double>(max_levels)) {
        return levels;  // lattice too fine relative to the data range to be a codebook
    }

    levels.level_start = static_cast<float>(start);
    levels.level_offset = level_offset;
    levels.level_num = static_cast<int>(needed);
    return levels;
}

/**
 * @brief Cluster/codebook quantizer over a uniform lattice of reconstruction levels.
 *
 * See the file-level comment for the exact error guarantee. Build the levels with
 * `derive_cluster_levels()` and hand them to the constructor.
 *
 * @tparam T Data type (floating point).
 */
template <class T>
class ClusterQuantizer : public concepts::QuantizerInterface<T, int> {
   public:
    /// @brief Construct an empty quantizer; only useful as a target for `load()`.
    ClusterQuantizer() = default;

    /**
     * @brief Construct from an explicit lattice.
     *
     * @param level_start Value of level 0.
     * @param level_offset Spacing between levels; must be finite and > 0.
     * @param level_num Number of levels; must be > 0.
     * @param eb Absolute error bound. Values <= 0 select the default `level_offset / 2`.
     * @throws std::invalid_argument if the lattice is not usable.
     */
    ClusterQuantizer(float level_start, float level_offset, int level_num, double eb = 0.0)
        : level_start_(level_start), level_offset_(level_offset), level_num_(level_num), eb_(eb) {
        validate();
    }

    /**
     * @brief Construct from levels produced by `derive_cluster_levels()`.
     *
     * @param levels Lattice codebook; must satisfy `levels.valid()`.
     * @param eb Absolute error bound. Values <= 0 select the default `level_offset / 2`.
     * @throws std::invalid_argument if `levels` is not valid.
     */
    explicit ClusterQuantizer(const ClusterLevels &levels, double eb = 0.0)
        : ClusterQuantizer(levels.level_start, levels.level_offset, levels.level_num, eb) {}

    /// @brief Absolute error bound actually enforced (the default is `level_offset / 2`).
    double get_eb() const { return eb_; }

    /// @brief Number of levels in the codebook.
    int get_level_num() const { return level_num_; }

    /// @brief Value of level 0.
    float get_level_start() const { return level_start_; }

    /// @brief Spacing between consecutive levels.
    float get_level_offset() const { return level_offset_; }

    /// @brief Number of samples that took the exact (unpredictable) path so far.
    size_t get_unpred_num() const { return unpred_.size(); }

    /**
     * @brief Reconstruction value of a level index.
     *
     * @param l Level index in `[0, level_num)`.
     * @return T The lattice point `level_start + l * level_offset`.
     */
    ALWAYS_INLINE T level(int l) const {
        return static_cast<T>(static_cast<double>(level_start_) + static_cast<double>(l) * level_offset_);
    }

    /**
     * @brief Map `data` onto the codebook, overwriting it with the reconstructed value.
     *
     * Falls back to storing `data` verbatim (bin 0) whenever the nearest in-range level is
     * farther away than the error bound, which keeps the guarantee unconditional. `data` is
     * left untouched on that path, exactly as `LinearQuantizer` does.
     *
     * @param data Value to quantize; overwritten by the reconstructed value on the codebook path.
     * @param pred Ignored -- the codebook is absolute, this quantizer does not use prediction.
     * @return int Bin: `level_index + 1`, or 0 for the exact/unpredictable path.
     */
    ALWAYS_INLINE int quantize_and_overwrite(T &data, T /*pred*/) override {
        // NaN fails every comparison below and therefore takes the exact path.
        const double x = static_cast<double>(data);
        double units = (x - static_cast<double>(level_start_)) / static_cast<double>(level_offset_);
        int l;
        if (units > 0.0) {
            units = std::round(units);
            l = (units > static_cast<double>(level_num_ - 1)) ? level_num_ - 1 : static_cast<int>(units);
        } else {
            l = 0;  // covers units <= 0 and NaN
        }
        const T rec = level(l);
        if (std::fabs(static_cast<double>(rec) - x) <= eb_) {
            data = rec;
            return l + 1;
        }
        unpred_.push_back(data);
        return 0;
    }

    /**
     * @brief Reconstruct a value from its bin.
     *
     * @param pred Ignored -- the codebook is absolute, this quantizer does not use prediction.
     * @param quant_index Bin from `quantize_and_overwrite()`; 0 pops the next unpredictable value.
     * @return T The reconstructed value.
     */
    ALWAYS_INLINE T recover(T /*pred*/, int quant_index) override {
        if (quant_index > 0) {
            return level(quant_index - 1);
        }
        return recover_unpred();
    }

    /// @brief Pop the next value from the unpredictable list (bin 0 path).
    ALWAYS_INLINE T recover_unpred() {
        if (index_ >= unpred_.size()) {
            throw std::runtime_error("ClusterQuantizer: more unpredictable bins than stored values");
        }
        return unpred_[index_++];
    }

    /**
     * @brief Store `ori` exactly, bypassing the codebook.
     *
     * @param ori Value to store verbatim.
     * @return int Always 0, the unpredictable sentinel.
     */
    int force_save_unpred(T ori) override {
        unpred_.push_back(ori);
        return 0;
    }

    /**
     * @brief Output bin range, `[first, second)`. `first` is 0 as required by SZ3.
     *
     * @return std::pair<int, int> `{0, level_num + 1}` -- bin 0 plus one bin per level.
     */
    std::pair<int, int> get_out_range() const override { return std::make_pair(0, level_num_ + 1); }

    /// @brief Estimated serialized size of the unpredictable payload, in bytes.
    size_t size_est() const { return unpred_.size() * sizeof(T); }

    /**
     * @brief Serialize the codebook, the error bound and the unpredictable list.
     *
     * @param c Buffer pointer; advanced past the written bytes.
     */
    void save(uchar *&c) const override {
        write(uid_, c);
        write(level_start_, c);
        write(level_offset_, c);
        write(level_num_, c);
        write(eb_, c);
        size_t unpred_size = unpred_.size();
        write(unpred_size, c);
        if (unpred_size > 0) {
            write(unpred_.data(), unpred_size, c);
        }
    }

    /**
     * @brief Deserialize the state written by `save()`.
     *
     * @param c Buffer pointer; advanced past the bytes read.
     * @param remaining_length Remaining readable bytes; decremented accordingly.
     * @throws std::invalid_argument on uid mismatch or an unusable lattice.
     */
    void load(const uchar *&c, size_t &remaining_length) override {
        uchar uid_read = 0;
        read(uid_read, c, remaining_length);
        if (uid_read != uid_) {
            throw std::invalid_argument("ClusterQuantizer uid mismatch");
        }
        read(level_start_, c, remaining_length);
        read(level_offset_, c, remaining_length);
        read(level_num_, c, remaining_length);
        read(eb_, c, remaining_length);
        size_t unpred_size = 0;
        read(unpred_size, c, remaining_length);
        unpred_.clear();
        if (unpred_size > 0) {
            if (unpred_size > remaining_length / sizeof(T)) {
                throw std::invalid_argument("ClusterQuantizer: unpredictable list exceeds buffer");
            }
            unpred_.resize(unpred_size);
            read(unpred_.data(), unpred_size, c, remaining_length);
        }
        index_ = 0;
        validate_loaded();
    }

    /// @brief Rewind the unpredictable-list cursor before a decompression pass.
    void predecompress_data() override { index_ = 0; }

    /// @brief Print a one-line summary of the codebook.
    void print() override {
        printf("[ClusterQuantizer] start=%.8G offset=%.8G levels=%d eb=%.8G unpred=%zu\n", level_start_, level_offset_,
               level_num_, eb_, unpred_.size());
    }

   private:
    void validate() {
        if (!(level_offset_ > 0.0f) || !std::isfinite(level_offset_) || !std::isfinite(level_start_)) {
            throw std::invalid_argument("ClusterQuantizer: level_offset must be finite and positive.");
        }
        if (level_num_ <= 0) {
            throw std::invalid_argument("ClusterQuantizer: level_num must be positive.");
        }
        if (!(eb_ > 0.0)) {
            eb_ = 0.5 * static_cast<double>(level_offset_);
        }
    }

    /// Same checks as `validate()` but without the eb defaulting -- eb comes from the stream.
    void validate_loaded() const {
        if (!(level_offset_ > 0.0f) || !std::isfinite(level_offset_) || !std::isfinite(level_start_) ||
            level_num_ <= 0 || !(eb_ > 0.0)) {
            throw std::invalid_argument("ClusterQuantizer: corrupted codebook in stream.");
        }
    }

    std::vector<T> unpred_;
    size_t index_ = 0;  // used in decompression only

    float level_start_ = 0.0f;
    float level_offset_ = 0.0f;
    int level_num_ = 0;
    double eb_ = 0.0;
    // Distinct from Linear (0b10), FixedPoint (0b11), BitTruncation (0b101), Level (0b1010).
    static constexpr uchar uid_ = 0b1001;
};

}  // namespace SZ3

#endif
