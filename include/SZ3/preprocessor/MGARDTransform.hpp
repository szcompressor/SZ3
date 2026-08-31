#ifndef SZ3_MGARD_TRANSFORM_HPP
#define SZ3_MGARD_TRANSFORM_HPP

/**
 * @file MGARDTransform.hpp
 * @ingroup Preprocessor
 * @brief In-place MGARD multigrid transform, with **no** quantization policy attached.
 *
 * The transform half of `MGARDFusedDecomposition`, split out so the multigrid basis change can be
 * reused on its own. `preprocess()` runs the forward (analysis) pass and `postprocess()` the
 * inverse (synthesis) pass, both in place on a row-major buffer of shape `dims`.
 *
 * After the forward pass the buffer holds `target_level + 1` coefficient slabs laid out as
 * nested nodal sub-blocks: `level_dims()[0]` is the coarsest nodal grid, and level `l` (for
 * `l >= 1`) owns everything inside `level_dims()[l]` that is not already inside
 * `level_dims()[l - 1]`. Those slabs partition the buffer exactly.
 *
 * Level-count policy lives here as well:
 *   - `max_level()` -- the deepest level the bundled MGARDx kernels accept.
 *   - `default_target_level()` -- the MGARDx-AC default, one level shallower and capped at 8.
 *
 * `kLevelEbGrowth` and `kLevelEbAmplification` are exposed as static members for
 * `geometric_level_ebs()`.
 *
 * Constraints: floating-point T (float/double); 1D / 2D / 3D.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/preprocessor/PreProcessor.hpp"
#include "SZ3/utils/thirdparty/mgard/MGARDHeaderOnly.hpp"

namespace SZ3 {

/**
 * @brief Forward/inverse MGARD multigrid transform as a standalone preprocessor.
 *
 * @tparam T Floating-point data type
 * @tparam N Data dimension (1, 2 or 3)
 */
template <class T, uint N>
class MGARDTransform : public concepts::PreprocessorInterface<T, N> {
    static_assert(std::is_floating_point<T>::value, "MGARDTransform requires a floating-point data type.");
    static_assert(N >= 1 && N <= 3, "MGARDTransform supports 1D, 2D, or 3D data only.");

   public:
    /**
     * @brief Operator norm of the piecewise-linear synthesis basis (`sqrt(8)` in 3D).
     *
     * Per-level error is amplified by at most this factor for every synthesis
     * step, so a per-level error budget that grows by this factor from coarse
     * to fine spends the total budget evenly. Consumed by `geometric_level_ebs()`.
     */
    static double level_eb_growth() { return std::sqrt(8.0); }

    /// Geometric constant `C2 = 1 + 3*sqrt(3)/4` of the MGARD error analysis.
    static double level_eb_amplification() { return 1.0 + 3.0 * std::sqrt(3.0) / 4.0; }

    /**
     * @brief Deepest level the bundled MGARDx kernels will actually run.
     *
     * `Decomposer::decompose()` silently clamps to this value but
     * `Recomposer::recompose()` does not, so callers must clamp themselves to
     * keep the two passes in agreement. `preprocess()`/`postprocess()` below do.
     */
    static size_t max_level(const std::vector<size_t> &dims) {
        const size_t min_dim = *std::min_element(dims.begin(), dims.end());
        if (min_dim <= 1) {
            return 0;
        }
        return static_cast<size_t>(std::log2(static_cast<double>(min_dim)));
    }

    /**
     * @brief Default level count: `floor(log2(min_dim)) - 1`, capped at 8.
     *
     * Mirrors MGARDx-AC's `max_level` policy — one level shallower than what the
     * kernels allow, because the very coarsest grid carries too few nodes to pay
     * for its own quantizer.
     */
    static size_t default_target_level(const std::vector<size_t> &dims) {
        const size_t min_dim = *std::min_element(dims.begin(), dims.end());
        if (min_dim <= 2) {
            return 0;
        }
        size_t level = static_cast<size_t>(std::floor(std::log2(static_cast<double>(min_dim)))) - 1;
        if (level > 8) {
            level = 8;  // AC default cap
        }
        return level;
    }

    /**
     * @brief Nested nodal sub-grid shape for every level.
     *
     * `level_dims[target_level]` is `dims` itself and each coarser entry halves
     * every axis as `(n >> 1) + 1`. Matches MGARDx-AC `init_levels`.
     */
    static std::vector<std::vector<size_t>> level_dims(const std::vector<size_t> &dims, size_t target_level) {
        std::vector<std::vector<size_t>> levels(target_level + 1, std::vector<size_t>(dims.size()));
        for (size_t i = 0; i < dims.size(); i++) {
            size_t n = dims[i];
            for (size_t j = 0; j <= target_level; j++) {
                levels[target_level - j][i] = n;
                n = (n >> 1) + 1;
            }
        }
        return levels;
    }

    /**
     * @brief Forward multigrid transform, in place.
     *
     * @param data Row-major buffer of shape `dims` (modified in place)
     * @param dims Data shape, `dims.size() == N`
     * @param target_level Number of multigrid levels; clamped to `max_level(dims)`
     */
    void preprocess(T *data, const std::vector<size_t> &dims, size_t target_level) {
        validate(data, dims);
        MGARD::Decomposer<T> decomposer;
        decomposer.decompose(data, dims, std::min(target_level, max_level(dims)), /*hierarchical=*/false);
    }

    /**
     * @brief Inverse multigrid transform, in place. Exact inverse of `preprocess()`.
     *
     * @param data Row-major buffer of shape `dims` (modified in place)
     * @param dims Data shape, `dims.size() == N`
     * @param target_level Must match the value handed to `preprocess()`
     */
    void postprocess(T *data, const std::vector<size_t> &dims, size_t target_level) {
        validate(data, dims);
        MGARD::Recomposer<T> recomposer;
        recomposer.recompose(data, dims, std::min(target_level, max_level(dims)), /*hierarchical=*/false);
    }

   private:
    static void validate(const T *data, const std::vector<size_t> &dims) {
        if (data == nullptr) {
            throw std::invalid_argument("MGARDTransform: null data pointer.");
        }
        if (dims.size() != N) {
            throw std::invalid_argument("MGARDTransform: dims vector size != N.");
        }
        if (std::any_of(dims.begin(), dims.end(), [](size_t d) { return d == 0; })) {
            throw std::invalid_argument("MGARDTransform: zero-sized dimension.");
        }
    }
};

}  // namespace SZ3

#endif
