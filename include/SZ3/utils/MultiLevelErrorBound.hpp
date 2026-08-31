#ifndef SZ3_MULTI_LEVEL_ERROR_BOUND_HPP
#define SZ3_MULTI_LEVEL_ERROR_BOUND_HPP

/**
 * @file MultiLevelErrorBound.hpp
 * @ingroup Utils
 * @brief Splitting one absolute error budget across the levels of a multi-resolution transform.
 *
 * Coarse coefficients are amplified by every subsequent synthesis step, so they need the tightest
 * bound; fine coefficients are barely amplified and can be quantized loosely, which is what makes
 * the bin stream cheap to entropy-code. This header holds that schedule as a standalone policy,
 * independent of the transform and of the quantizer that consumes it. The transform supplies the
 * two constants (see `MGARDTransform::level_eb_growth()` and
 * `MGARDTransform::level_eb_amplification()`).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "SZ3/def.hpp"

namespace SZ3 {

/**
 * @brief Geometric per-level error-bound schedule.
 *
 * Solves the geometric budget so that the per-level bounds sum to
 * `eb / amplification`, growing by `growth` per level from coarsest to finest:
 *
 *     level_eb[0] = (1 - g) / (1 - g^(L+1)) * eb / A
 *     level_eb[l] = level_eb[l - 1] * g                       (l = 1..L)
 *
 * so that `sum_l level_eb[l] == eb / A`. With `g` set to the synthesis operator
 * norm and `A` to the basis constant, the worst-case amplified sum of the
 * per-level errors is bounded by `eb`.
 *
 * @param eb           Total absolute error budget (must be > 0)
 * @param target_level Deepest level index `L`; the result has `L + 1` entries
 * @param growth       Per-level growth factor `g` (must be > 0)
 * @param amplification Basis constant `A` (must be > 0)
 * @return Per-level absolute error bounds, index 0 = coarsest
 */
inline std::vector<double> geometric_level_ebs(double eb, size_t target_level, double growth, double amplification) {
    if (!(eb > 0.0)) {
        throw std::invalid_argument("geometric_level_ebs: error bound must be positive.");
    }
    if (!(growth > 0.0)) {
        throw std::invalid_argument("geometric_level_ebs: growth factor must be positive.");
    }
    if (!(amplification > 0.0)) {
        throw std::invalid_argument("geometric_level_ebs: amplification must be positive.");
    }

    std::vector<double> ebs(target_level + 1);
    if (growth == 1.0) {
        // Degenerate case: the geometric series collapses to an even split.
        const double share = eb / (amplification * static_cast<double>(target_level + 1));
        std::fill(ebs.begin(), ebs.end(), share);
        return ebs;
    }
    const double head = (1.0 - growth) / (1.0 - std::pow(growth, static_cast<double>(target_level + 1)));
    ebs[0] = head * eb / amplification;
    for (size_t l = 1; l <= target_level; l++) {
        ebs[l] = ebs[l - 1] * growth;
    }
    return ebs;
}

}  // namespace SZ3

#endif
