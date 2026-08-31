#ifndef SZ3_MULTI_LEVEL_DECOMPOSITION_HPP
#define SZ3_MULTI_LEVEL_DECOMPOSITION_HPP

/**
 * @file MultiLevelDecomposition.hpp
 * @ingroup Decomposition
 * @brief Multi-resolution transform + per-level error budget + per-level
 *        quantization, assembled from three independent modules.
 *
 * This is the composable counterpart of `MGARDFusedDecomposition`, which fuses the
 * same three concerns into one class. Here each is a separate, substitutable
 * piece:
 *
 *  - the **transform** (`Transform`, default `MGARDTransform<T, N>`) supplies
 *    the forward/inverse basis change plus its level geometry;
 *  - the **schedule** (`geometric_level_ebs()` in `MultiLevelErrorBound.hpp`)
 *    splits one absolute error bound across the levels;
 *  - the **quantizer** (`Quantizer`, default `LinearQuantizer<T>`) is built once
 *    per level by a caller-supplied factory, so the multigrid transform can be
 *    paired with any `concepts::QuantizerInterface` implementation.
 *
 * Wired with its defaults, this reproduces `MGARDFusedDecomposition` bit for bit.
 *
 * Example — MGARD multigrid driven by a non-default quantizer:
 * @code
 *   using Q = SZ3::LevelQuantizer<float>;
 *   SZ3::MultiLevelDecomposition<float, 3, Q> decomp(
 *       eb, [](double level_eb) { return Q(level_eb, 32768, SZ3::LevelCurve::Log); });
 *   auto bins = decomp.compress(conf, data);
 * @endcode
 *
 * The quantizer factory must produce quantizers with a uniform output range
 * whose lower end is 0, as `SZGenericCompressor` requires; this is checked at
 * construction.
 *
 * @note Choosing a quantizer: the coarsest level holds *nodal values*, not
 * residuals, so its magnitudes are O(data) while its error bound is the
 * smallest of the schedule. A quantizer with a bounded representable range and
 * no unpredictable-value escape hatch therefore silently clamps level-0
 * coefficients once the hierarchy is deep enough, and the error bound is lost.
 * Cap `target_level` for such quantizers. Quantizers with an unpredictable path
 * -- `LinearQuantizer`, `LevelQuantizer` -- store what they cannot represent
 * verbatim and hold their bound at any depth, at some cost in size.
 */

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/preprocessor/MGARDTransform.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/MultiLevelErrorBound.hpp"
#include "SZ3/utils/MultiLevelQuantization.hpp"

namespace SZ3 {

/// Detects the optional (non-virtual) `size_est()` some quantizers expose.
template <class Q, class = void>
struct quantizer_has_size_est : std::false_type {};
template <class Q>
struct quantizer_has_size_est<Q, std::void_t<decltype(std::declval<Q &>().size_est())>> : std::true_type {};

/**
 * @brief Decomposition over a coarse-to-fine transform with one quantizer per level.
 *
 * @tparam T Floating-point data type
 * @tparam N Data dimension (1, 2 or 3)
 * @tparam Quantizer Any `concepts::QuantizerInterface<T, To>` implementation
 * @tparam Transform Multi-resolution transform (see `MGARDTransform` for the required surface)
 * @tparam To Bin type, deduced from `Quantizer`
 */
template <class T, uint N, class Quantizer = LinearQuantizer<T>, class Transform = MGARDTransform<T, N>,
          class To = quantizer_bin_t<T, Quantizer>>
class MultiLevelDecomposition : public concepts::DecompositionInterface<T, To, N> {
   public:
    /// Builds the quantizer for one level from that level's absolute error bound.
    using QuantizerFactory = std::function<Quantizer(double level_eb)>;

    /// Pass as `target_level` to let the transform pick its own level count.
    static constexpr size_t kAutoTargetLevel = static_cast<size_t>(-1);

    /**
     * @brief Construct the composition.
     *
     * @param abs_error_bound Total absolute error bound (must be > 0)
     * @param make_quantizer Factory invoked once per level with that level's error bound
     * @param target_level Level count, or `kAutoTargetLevel` for the transform's default
     * @param level_eb_growth Per-level growth factor of the error schedule
     * @param level_eb_amplification Basis amplification constant of the error schedule
     */
    MultiLevelDecomposition(double abs_error_bound, QuantizerFactory make_quantizer,
                            size_t target_level = kAutoTargetLevel,
                            double level_eb_growth = Transform::level_eb_growth(),
                            double level_eb_amplification = Transform::level_eb_amplification())
        : eb_(abs_error_bound),
          make_quantizer_(std::move(make_quantizer)),
          requested_level_(target_level),
          growth_(level_eb_growth),
          amplification_(level_eb_amplification) {
        if (!(eb_ > 0.0)) {
            throw std::invalid_argument("MultiLevelDecomposition: error bound must be positive.");
        }
        if (!make_quantizer_) {
            throw std::invalid_argument("MultiLevelDecomposition: quantizer factory must be callable.");
        }
        out_range_ = make_quantizer_(eb_).get_out_range();
        if (out_range_.first != To(0)) {
            throw std::invalid_argument("MultiLevelDecomposition: the quantizer's output range must start from 0.");
        }
    }

    std::vector<To> compress(const Config &conf, T *data) override {
        dims_ = resolve_dims(conf);
        target_level_ = resolve_target_level(dims_);

        transform_.preprocess(data, dims_, target_level_);

        const auto levels = Transform::level_dims(dims_, target_level_);
        const auto level_ebs = geometric_level_ebs(eb_, target_level_, growth_, amplification_);

        quantizers_.clear();
        quantizers_.reserve(level_ebs.size());
        for (double level_eb : level_ebs) {
            quantizers_.push_back(make_quantizer_(level_eb));
        }

        std::vector<To> bins;
        bins.reserve(conf.num);
        multilevel_quantize<N>(data, dims_, levels, quantizers_, bins);
        return bins;
    }

    T *decompress(const Config &conf, std::vector<To> &bins, T *dec_data) override {
        dims_ = resolve_dims(conf);
        const auto levels = Transform::level_dims(dims_, target_level_);

        // The level slabs partition the buffer, so this fill is belt-and-braces:
        // it keeps a truncated level set from leaking uninitialised memory.
        std::fill(dec_data, dec_data + conf.num, T(0));

        const size_t consumed = multilevel_recover<N>(dec_data, dims_, levels, quantizers_, bins);
        if (consumed != bins.size()) {
            throw std::runtime_error("MultiLevelDecomposition: bin count mismatch in decompress.");
        }

        transform_.postprocess(dec_data, dims_, target_level_);
        return dec_data;
    }

    void save(uchar *&c) override {
        write(target_level_, c);
        write(eb_, c);
        const size_t num_quantizers = quantizers_.size();
        write(num_quantizers, c);
        for (auto &q : quantizers_) {
            q.save(c);
        }
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        read(target_level_, c, remaining_length);
        read(eb_, c, remaining_length);
        size_t num_quantizers = 0;
        read(num_quantizers, c, remaining_length);
        quantizers_.assign(num_quantizers, make_quantizer_(eb_));
        for (auto &q : quantizers_) {
            q.load(c, remaining_length);
        }
        if (!quantizers_.empty()) {
            out_range_ = quantizers_.front().get_out_range();
        }
    }

    std::pair<To, To> get_out_range() override { return out_range_; }

    size_t size_est() override {
        size_t bytes = sizeof(target_level_) + sizeof(eb_) + sizeof(size_t) + 64 * (quantizers_.size() + 1);
        if constexpr (quantizer_has_size_est<Quantizer>::value) {
            for (auto &q : quantizers_) {
                bytes += q.size_est();
            }
        } else {
            // No introspection available: assume a quantizer stores at most one
            // value per input element (LinearQuantizer's worst case).
            bytes += num_elements() * sizeof(T);
        }
        return bytes;
    }

    void print() override {
        printf("[MultiLevelDecomposition] eb=%.8G target_level=%zu quantizers=%zu\n", eb_, target_level_,
               quantizers_.size());
        for (auto &q : quantizers_) {
            q.print();
        }
    }

   private:
    static std::vector<size_t> resolve_dims(const Config &conf) {
        if (conf.N != N) {
            throw std::invalid_argument("MultiLevelDecomposition: dimensionality mismatch.");
        }
        if (conf.dims.size() != N) {
            throw std::invalid_argument("MultiLevelDecomposition: dims vector size != N.");
        }
        return std::vector<size_t>(conf.dims.begin(), conf.dims.end());
    }

    size_t resolve_target_level(const std::vector<size_t> &dims) const {
        const size_t level =
            (requested_level_ == kAutoTargetLevel) ? Transform::default_target_level(dims) : requested_level_;
        return std::min(level, Transform::max_level(dims));
    }

    size_t num_elements() const {
        size_t n = 1;
        for (size_t d : dims_) {
            n *= d;
        }
        return n;
    }

    double eb_;
    QuantizerFactory make_quantizer_;
    size_t requested_level_;
    double growth_;
    double amplification_;

    Transform transform_;
    size_t target_level_ = 0;
    std::vector<size_t> dims_;
    std::vector<Quantizer> quantizers_;
    std::pair<To, To> out_range_{To(0), To(0)};
};

/**
 * @brief MGARD multigrid as a composable decomposition: the default MGARD wiring.
 *
 * `MGARDFusedDecomposition` is the fused equivalent, kept for the `ALGO_MGARD` pipeline.
 */
template <class T, uint N, class Quantizer = LinearQuantizer<T>>
using MGARDDecomposition = MultiLevelDecomposition<T, N, Quantizer, MGARDTransform<T, N>>;

}  // namespace SZ3

#endif
