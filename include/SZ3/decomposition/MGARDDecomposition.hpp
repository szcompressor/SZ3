#ifndef SZ3_MGARD_PERLEVEL_DECOMPOSITION_HPP
#define SZ3_MGARD_PERLEVEL_DECOMPOSITION_HPP

/**
 * @file MGARDDecomposition.hpp
 * @ingroup Decomposition
 * @brief MGARD multigrid + **per-level** LinearQuantizer (one quantizer per
 *        multigrid level, with geometrically-scaled error bound).
 *
 * Mirrors the MGARDx-AC `puremg` strategy: after the in-place forward
 * multigrid transform, the L+1 coefficient slabs are quantized with
 * per-level error bounds that grow geometrically from the coarsest to the
 * finest level. The coarse-level coefficients (few of them but with the
 * largest magnitudes and the largest amplification factors during inverse
 * synthesis) get tighter quantization; the fine-level coefficients (many
 * of them, with the smallest amplification) get looser quantization. This
 * matches the multigrid energy budget and produces a quant-index stream
 * that's much cheaper to Huffman-code than a single global quantizer's.
 *
 * Theory: the MGARD inverse synthesis amplifies per-level error by at most
 * `c^L * C2`, where `c = sqrt(8) ≈ 2.83` is the 3D operator norm of the
 * piecewise-linear synthesis basis, and `C2 = 1 + 3*sqrt(3)/4 ≈ 2.30` is a
 * geometric constant. Solving the geometric-series budget so that the sum
 * of per-level contributions equals `eb` gives:
 *
 *     level_eb_l = (1 - c) / (1 - c^(L+1)) * eb / C2 * c^l   (l = 0..L)
 *
 * The output is a single concatenated `std::vector<int>` of all per-level
 * quantization indices (level 0 first, then level 1, etc.). The order
 * walks the same ZYX iteration as MGARDx-AC `quantize_level`. Pairs
 * naturally with `HuffmanEncoder<int>`.
 *
 * Constraints: floating-point T (float/double); 1D / 2D / 3D.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/thirdparty/mgard/MGARDHeaderOnly.hpp"

namespace SZ3 {

template <class T, uint N>
class MGARDDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    explicit MGARDDecomposition(double abs_error_bound, int radius = 32768)
        : eb_(abs_error_bound), radius_(radius) {
        static_assert(std::is_floating_point<T>::value,
                      "MGARDDecomposition requires a floating-point data type.");
        static_assert(N >= 1 && N <= 3,
                      "MGARDDecomposition supports 1D, 2D, or 3D data only.");
        if (!(eb_ > 0.0)) {
            throw std::invalid_argument("MGARDDecomposition: error bound must be positive.");
        }
    }

    std::vector<int> compress(const Config& conf, T* data) override {
        const std::vector<size_t> dims = resolve_dims(conf);
        target_level_ = compute_target_level(dims);
        dims_ = dims;

        // Forward multigrid transform (in place).
        MGARD::Decomposer<T> decomposer;
        decomposer.decompose(data, dims, target_level_, /*hierarchical=*/false);

        // Build per-level dim slabs and per-level eb.
        const auto level_dims = init_levels(dims, target_level_);
        const auto level_ebs = build_level_ebs(eb_, target_level_);

        // Quantize per-level. Reserve one quantizer per level, save them all.
        quantizers_.clear();
        quantizers_.reserve(target_level_ + 1);
        std::vector<int> bins;
        bins.reserve(conf.num);

        // Level 0: ALL of level_dims[0] (coarsest nodal subgrid).
        std::vector<size_t> dummy_dims(dims.size(), 0);
        quantize_level_walk(data, dims, dummy_dims, level_dims[0], level_ebs[0], bins);

        // Levels 1..L: difference region between level_dims[l-1] and level_dims[l].
        for (size_t l = 1; l <= target_level_; l++) {
            quantize_level_walk(data, dims, level_dims[l - 1], level_dims[l], level_ebs[l], bins);
        }
        return bins;
    }

    T* decompress(const Config& conf, std::vector<int>& bins, T* dec_data) override {
        const std::vector<size_t> dims = resolve_dims(conf);
        if (dims_ != dims) {
            // dims wasn't restored from save() (compute_target_level relies on
            // dims_), but resolve_dims already returns conf.dims; fall through.
            dims_ = dims;
        }
        const auto level_dims = init_levels(dims, target_level_);

        // Zero-init the buffer (per-level walks only write the coefficient region).
        std::fill(dec_data, dec_data + conf.num, T(0));

        size_t bin_off = 0;
        std::vector<size_t> dummy_dims(dims.size(), 0);
        bin_off = recover_level_walk(dec_data, dims, dummy_dims, level_dims[0], 0, bins, bin_off);
        for (size_t l = 1; l <= target_level_; l++) {
            bin_off = recover_level_walk(dec_data, dims, level_dims[l - 1], level_dims[l], l, bins, bin_off);
        }
        if (bin_off != bins.size()) {
            throw std::runtime_error("MGARDDecomposition: bin count mismatch in decompress.");
        }
        for (auto& q : quantizers_) q.postdecompress_data();

        MGARD::Recomposer<T> recomposer;
        recomposer.recompose(dec_data, dims, target_level_, /*hierarchical=*/false);
        return dec_data;
    }

    void save(uchar*& c) override {
        write(target_level_, c);
        write(eb_, c);
        write(radius_, c);
        size_t nq = quantizers_.size();
        write(nq, c);
        for (auto& q : quantizers_) q.save(c);
    }

    void load(const uchar*& c, size_t& remaining_length) override {
        read(target_level_, c, remaining_length);
        read(eb_, c, remaining_length);
        read(radius_, c, remaining_length);
        size_t nq = 0;
        read(nq, c, remaining_length);
        quantizers_.assign(nq, LinearQuantizer<T>(eb_, radius_));
        for (auto& q : quantizers_) q.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return std::make_pair(0, radius_ * 2); }

    size_t size_est() override {
        // Header + per-quantizer estimate (each quantizer carries unpredictable list).
        size_t s = sizeof(target_level_) + sizeof(eb_) + sizeof(radius_) + sizeof(size_t);
        for (auto& q : quantizers_) s += q.size_est() + 64;
        return s + 64;
    }

   private:
    static std::vector<size_t> resolve_dims(const Config& conf) {
        if (conf.N != N) {
            throw std::invalid_argument("MGARDDecomposition: dimensionality mismatch.");
        }
        if (conf.dims.size() != N) {
            throw std::invalid_argument("MGARDDecomposition: dims vector size != N.");
        }
        return std::vector<size_t>(conf.dims.begin(), conf.dims.end());
    }

    // floor(log2(min_dim)) - 1 to mirror MGARDx-AC's `max_level` cap.
    static size_t compute_target_level(const std::vector<size_t>& dims) {
        const size_t min_dim = *std::min_element(dims.begin(), dims.end());
        if (min_dim <= 2) return 0;
        size_t lvl = static_cast<size_t>(std::floor(std::log2(static_cast<double>(min_dim)))) - 1;
        if (lvl > 8) lvl = 8;  // AC default cap
        return lvl;
    }

    // Match MGARDx-AC `init_levels`: level_dims[L] = full dims, level_dims[l-1] = (dims_l>>1)+1.
    static std::vector<std::vector<size_t>> init_levels(const std::vector<size_t>& dims, size_t target_level) {
        std::vector<std::vector<size_t>> level_dims(target_level + 1, std::vector<size_t>(dims.size()));
        for (size_t i = 0; i < dims.size(); i++) {
            size_t n = dims[i];
            for (size_t j = 0; j <= target_level; j++) {
                level_dims[target_level - j][i] = n;
                n = (n >> 1) + 1;
            }
        }
        return level_dims;
    }

    // Geometric per-level eb schedule from MGARDx-AC.
    //   level_eb[0] = (1 - c) / (1 - c^(L+1)) * eb / C2
    //   level_eb[l] = level_eb[l-1] * c
    static std::vector<double> build_level_ebs(double eb, size_t target_level) {
        const double c = std::sqrt(8.0);
        const double C2 = 1.0 + 3.0 * std::sqrt(3.0) / 4.0;
        const double cc = (1.0 - c) / (1.0 - std::pow(c, static_cast<double>(target_level + 1)));
        std::vector<double> ebs(target_level + 1);
        ebs[0] = cc * eb / C2;
        for (size_t l = 1; l <= target_level; l++) ebs[l] = ebs[l - 1] * c;
        return ebs;
    }

    // Walk the coefficient region of level l (everything in fine_dims that is
    // NOT in coarse_dims) and quantize. Mirrors `quantize_level` in MGARDx-AC.
    void quantize_level_walk(T* data, const std::vector<size_t>& dims,
                             const std::vector<size_t>& coarse_dims,
                             const std::vector<size_t>& fine_dims, double level_eb,
                             std::vector<int>& bins) {
        quantizers_.emplace_back(level_eb, radius_);
        auto& quantizer = quantizers_.back();
        if (N == 1) {
            const size_t s0 = (coarse_dims.size() > 0) ? coarse_dims[0] : 0;
            for (size_t i = 0; i < fine_dims[0]; i++) {
                if (i < s0) continue;
                bins.push_back(quantizer.quantize_and_overwrite(data[i], T(0)));
            }
        } else if (N == 2) {
            const size_t c0 = (coarse_dims.size() > 0) ? coarse_dims[0] : 0;
            const size_t c1 = (coarse_dims.size() > 1) ? coarse_dims[1] : 0;
            for (size_t i = 0; i < fine_dims[0]; i++) {
                for (size_t j = 0; j < fine_dims[1]; j++) {
                    if (i < c0 && j < c1) continue;
                    bins.push_back(quantizer.quantize_and_overwrite(data[i * dims[1] + j], T(0)));
                }
            }
        } else {  // N == 3
            const size_t c0 = (coarse_dims.size() > 0) ? coarse_dims[0] : 0;
            const size_t c1 = (coarse_dims.size() > 1) ? coarse_dims[1] : 0;
            const size_t c2 = (coarse_dims.size() > 2) ? coarse_dims[2] : 0;
            const size_t stride0 = dims[1] * dims[2];
            for (size_t i = 0; i < fine_dims[0]; i++) {
                for (size_t j = 0; j < fine_dims[1]; j++) {
                    for (size_t k = 0; k < fine_dims[2]; k++) {
                        if (i < c0 && j < c1 && k < c2) continue;
                        bins.push_back(
                            quantizer.quantize_and_overwrite(data[i * stride0 + j * dims[2] + k], T(0)));
                    }
                }
            }
        }
        quantizer.postcompress_data();
    }

    // Inverse walk: recover decompressed values from bins[] for the level-l region.
    size_t recover_level_walk(T* dec_data, const std::vector<size_t>& dims,
                              const std::vector<size_t>& coarse_dims,
                              const std::vector<size_t>& fine_dims, size_t level_idx,
                              std::vector<int>& bins, size_t bin_off) {
        if (level_idx >= quantizers_.size()) {
            throw std::runtime_error("MGARDDecomposition: missing quantizer for level.");
        }
        auto& quantizer = quantizers_[level_idx];
        if (N == 1) {
            const size_t s0 = (coarse_dims.size() > 0) ? coarse_dims[0] : 0;
            for (size_t i = 0; i < fine_dims[0]; i++) {
                if (i < s0) continue;
                dec_data[i] = quantizer.recover(T(0), bins[bin_off++]);
            }
        } else if (N == 2) {
            const size_t c0 = (coarse_dims.size() > 0) ? coarse_dims[0] : 0;
            const size_t c1 = (coarse_dims.size() > 1) ? coarse_dims[1] : 0;
            for (size_t i = 0; i < fine_dims[0]; i++) {
                for (size_t j = 0; j < fine_dims[1]; j++) {
                    if (i < c0 && j < c1) continue;
                    dec_data[i * dims[1] + j] = quantizer.recover(T(0), bins[bin_off++]);
                }
            }
        } else {  // N == 3
            const size_t c0 = (coarse_dims.size() > 0) ? coarse_dims[0] : 0;
            const size_t c1 = (coarse_dims.size() > 1) ? coarse_dims[1] : 0;
            const size_t c2 = (coarse_dims.size() > 2) ? coarse_dims[2] : 0;
            const size_t stride0 = dims[1] * dims[2];
            for (size_t i = 0; i < fine_dims[0]; i++) {
                for (size_t j = 0; j < fine_dims[1]; j++) {
                    for (size_t k = 0; k < fine_dims[2]; k++) {
                        if (i < c0 && j < c1 && k < c2) continue;
                        dec_data[i * stride0 + j * dims[2] + k] =
                            quantizer.recover(T(0), bins[bin_off++]);
                    }
                }
            }
        }
        return bin_off;
    }

    double eb_;
    int radius_;
    size_t target_level_ = 0;
    std::vector<size_t> dims_;
    std::vector<LinearQuantizer<T>> quantizers_;
};

}  // namespace SZ3

#endif
