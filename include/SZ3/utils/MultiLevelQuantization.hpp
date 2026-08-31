#ifndef SZ3_MULTI_LEVEL_QUANTIZATION_HPP
#define SZ3_MULTI_LEVEL_QUANTIZATION_HPP

/**
 * @file MultiLevelQuantization.hpp
 * @ingroup Utils
 * @brief Walking the nested coefficient slabs of a multi-resolution transform
 *        and driving one quantizer per level.
 *
 * A coarse-to-fine transform (multigrid, dyadic wavelet, ...) leaves its output
 * in a *nested nodal* layout: level `L` is the whole array, and each coarser
 * level occupies the sub-block anchored at the origin whose extent is given by
 * `level_dims[l]`. The coefficients *owned* by level `l` are therefore the
 * elements inside `level_dims[l]` that are not already inside `level_dims[l-1]`
 * (level 0 owns its whole sub-block). Those slabs partition the array exactly.
 *
 * `for_each_level_coefficient()` is the raw walk over one such slab.
 * `multilevel_quantize()` / `multilevel_recover()` layer per-level quantization
 * on top of it, against any `concepts::QuantizerInterface` implementation — the
 * transform, the error-budget policy and the quantizer stay independent.
 *
 * Iteration is outermost-dimension-first (ZYX) so the produced bin stream has a
 * stable, transform-independent order.
 *
 * All quantization is done against a zero prediction: these are transform
 * coefficients, not prediction residuals.
 */

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "SZ3/def.hpp"

namespace SZ3 {

/**
 * @brief Visit the flat offsets of every coefficient owned by one level.
 *
 * @tparam N Data dimension (1, 2 or 3)
 * @param dims Full data shape (supplies the row-major strides)
 * @param coarse_dims Shape of the enclosed coarser sub-block; pass all-zeros
 *        (or an empty vector) for the coarsest level, which owns everything
 *        inside `fine_dims`
 * @param fine_dims Shape of this level's sub-block
 * @param fn Callable invoked as `fn(size_t offset)` for each owned element
 */
template <uint N, class F>
void for_each_level_coefficient(const std::vector<size_t> &dims, const std::vector<size_t> &coarse_dims,
                                const std::vector<size_t> &fine_dims, F &&fn) {
    if (dims.size() != N || fine_dims.size() != N) {
        throw std::invalid_argument("for_each_level_coefficient: dims/fine_dims size != N.");
    }
    const auto coarse = [&coarse_dims](size_t axis) -> size_t {
        return axis < coarse_dims.size() ? coarse_dims[axis] : 0;
    };

    if constexpr (N == 1) {
        const size_t c0 = coarse(0);
        for (size_t i = 0; i < fine_dims[0]; i++) {
            if (i < c0) continue;
            fn(i);
        }
    } else if constexpr (N == 2) {
        const size_t c0 = coarse(0);
        const size_t c1 = coarse(1);
        for (size_t i = 0; i < fine_dims[0]; i++) {
            for (size_t j = 0; j < fine_dims[1]; j++) {
                if (i < c0 && j < c1) continue;
                fn(i * dims[1] + j);
            }
        }
    } else {
        const size_t c0 = coarse(0);
        const size_t c1 = coarse(1);
        const size_t c2 = coarse(2);
        const size_t stride0 = dims[1] * dims[2];
        for (size_t i = 0; i < fine_dims[0]; i++) {
            for (size_t j = 0; j < fine_dims[1]; j++) {
                for (size_t k = 0; k < fine_dims[2]; k++) {
                    if (i < c0 && j < c1 && k < c2) continue;
                    fn(i * stride0 + j * dims[2] + k);
                }
            }
        }
    }
}

/**
 * @brief Quantize every level of a transformed buffer with its own quantizer.
 *
 * `quantizers[l]` handles level `l` and is expected to already carry that
 * level's error bound (see `geometric_level_ebs()`). Bins are *appended* to
 * `bins`, coarsest level first. `data` is overwritten in place with the
 * reconstructed coefficients, as `QuantizerInterface` prescribes.
 *
 * Each quantizer's `precompress_data()` / `postcompress_data()` hooks are run
 * around its own level.
 *
 * @tparam N Data dimension (1, 2 or 3)
 * @param data Transformed buffer of shape `dims`, modified in place
 * @param dims Full data shape
 * @param level_dims Per-level sub-block shapes, coarsest first (see the transform)
 * @param quantizers One quantizer per level; `quantizers.size() == level_dims.size()`
 * @param bins Output bin stream, appended to
 */
template <uint N, class T, class Quantizer, class To>
void multilevel_quantize(T *data, const std::vector<size_t> &dims, const std::vector<std::vector<size_t>> &level_dims,
                         std::vector<Quantizer> &quantizers, std::vector<To> &bins) {
    if (quantizers.size() != level_dims.size()) {
        throw std::invalid_argument("multilevel_quantize: one quantizer per level is required.");
    }
    const std::vector<size_t> no_coarse(dims.size(), 0);
    for (size_t l = 0; l < level_dims.size(); l++) {
        Quantizer &quantizer = quantizers[l];
        quantizer.precompress_data();
        for_each_level_coefficient<N>(
            dims, (l == 0) ? no_coarse : level_dims[l - 1], level_dims[l],
            [&](size_t offset) { bins.push_back(quantizer.quantize_and_overwrite(data[offset], T(0))); });
        quantizer.postcompress_data();
    }
}

/**
 * @brief Inverse of `multilevel_quantize()`: rebuild the transformed buffer.
 *
 * Reads `bins` starting at `bin_offset` and writes the reconstructed
 * coefficients into `dec_data`. Because the level slabs partition the array,
 * every element of `dec_data` is written exactly once.
 *
 * @return The offset one past the last bin consumed
 */
template <uint N, class T, class Quantizer, class To>
size_t multilevel_recover(T *dec_data, const std::vector<size_t> &dims,
                          const std::vector<std::vector<size_t>> &level_dims, std::vector<Quantizer> &quantizers,
                          const std::vector<To> &bins, size_t bin_offset = 0) {
    if (quantizers.size() != level_dims.size()) {
        throw std::invalid_argument("multilevel_recover: one quantizer per level is required.");
    }
    const std::vector<size_t> no_coarse(dims.size(), 0);
    for (size_t l = 0; l < level_dims.size(); l++) {
        Quantizer &quantizer = quantizers[l];
        quantizer.predecompress_data();
        for_each_level_coefficient<N>(dims, (l == 0) ? no_coarse : level_dims[l - 1], level_dims[l],
                                      [&](size_t offset) {
                                          if (bin_offset >= bins.size()) {
                                              throw std::runtime_error("multilevel_recover: bin stream underrun.");
                                          }
                                          dec_data[offset] = quantizer.recover(T(0), bins[bin_offset++]);
                                      });
        quantizer.postdecompress_data();
    }
    return bin_offset;
}

}  // namespace SZ3

#endif
