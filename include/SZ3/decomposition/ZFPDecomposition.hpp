/**
 * @file ZFPDecomposition.hpp
 * @ingroup Decomposition
 */

#ifndef SZ3_ZFP_DECOMPOSITION_HPP
#define SZ3_ZFP_DECOMPOSITION_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/thirdparty/zfp/zfp_codec.hpp"

namespace SZ3 {

/**
 * @brief zfp CODEC 5's transform stage: block-floating-point alignment, the decorrelating
 *        transform, and the reordering into negabinary coefficients.
 *
 * Emits the per-block exponents followed by the coefficients, which ZFPEncoder consumes.
 * Carries no error bound of its own: the bound enters downstream, where each block's
 * precision is derived from its exponent. Pair it with ZFPEncoder -- the two are one codec
 * split in half and neither works with anything else.
 */
template <class T, class To, uint N>
class ZFPDecomposition : public concepts::DecompositionInterface<T, To, N> {
    static_assert(N >= 1 && N <= 3, "zfp modules are instantiated for 1D, 2D and 3D");
    static_assert(sizeof(To) == sizeof(T), "the bin type must be the integer of the scalar's width");

    using ops = ZFP::ops<T, N>;
    static constexpr int bs = ops::block_size;

    /// Number of 4^N blocks, and the strides zfp walks the array with.
    static size_t geometry(const Config &conf, std::array<size_t, N> &nb, std::array<ptrdiff_t, N> &stride) {
        // Config::setDims drops extents of 1, so conf.dims.size() can be smaller than N.
        // SZ_compress dispatches on conf.N, but a direct caller can hand over a mismatch.
        if (conf.dims.size() != N) {
            throw std::invalid_argument("SZ3 zfp: the configuration's dimensionality does not match the module's");
        }
        // SZ3 is row-major (conf.dims[N-1] is contiguous); zfp numbers axes from the fastest.
        // Reverse the mapping, as SZAlgoSPERR does. Getting this wrong still round-trips -- it
        // just stops the 4^N blocks being neighbourhoods, which costs ratio.
        size_t total = 1;
        ptrdiff_t s = 1;
        for (uint i = 0; i < N; i++) {
            const uint axis = N - 1 - i;  // zfp axis i is SZ3 dimension N-1-i
            nb[i] = (conf.dims[axis] + 3) / 4;
            stride[i] = s;
            s *= static_cast<ptrdiff_t>(conf.dims[axis]);
            total *= nb[i];
        }
        return total;
    }

    /// Walks the array in zfp's own block order, handing each block to `op`.
    template <class Op>
    static void for_each_block(const Config &conf, Op op) {
        std::array<size_t, N> nb;
        std::array<ptrdiff_t, N> stride;
        geometry(conf, nb, stride);
        std::array<size_t, N> b{};
        size_t n_blocks = 1;
        for (uint i = 0; i < N; i++) n_blocks *= nb[i];
        for (size_t k = 0; k < n_blocks; k++) {
            ptrdiff_t offset = 0;
            std::array<size_t, N> extent;
            for (uint i = 0; i < N; i++) {
                const size_t x = b[i] * 4;
                offset += stride[i] * static_cast<ptrdiff_t>(x);
                extent[i] = std::min<size_t>(conf.dims[N - 1 - i] - x, 4);
            }
            op(offset, stride.data(), extent.data());
            // odometer, slowest dimension last -- the order zfp's own loops produce
            for (uint i = 0; i < N; i++) {
                if (++b[i] < nb[i]) break;
                b[i] = 0;
            }
        }
    }

   public:
    std::vector<To> compress(const Config &conf, T *data) override {
        std::array<size_t, N> nb;
        std::array<ptrdiff_t, N> stride;
        const size_t n_blocks = geometry(conf, nb, stride);
        std::vector<To> out(1 + n_blocks + n_blocks * bs);
        out[0] = static_cast<To>(n_blocks);
        To *emax_pos = &out[1];
        auto *coeff_pos = reinterpret_cast<typename ops::uint_type *>(&out[1 + n_blocks]);
        for_each_block(conf, [&](ptrdiff_t off, const ptrdiff_t *s, const size_t *n) {
            *emax_pos++ = static_cast<To>(ops::fwd(data + off, s, n, coeff_pos));
            coeff_pos += bs;
        });
        return out;
    }

    T *decompress(const Config &conf, std::vector<To> &transformed, T *dec_data) override {
        std::array<size_t, N> nb;
        std::array<ptrdiff_t, N> stride;
        const size_t n_blocks = geometry(conf, nb, stride);
        if (transformed.size() < 1 + n_blocks + n_blocks * bs ||
            static_cast<size_t>(transformed[0]) != n_blocks) {
            throw std::out_of_range("SZ3 zfp: coefficient stream does not match the configuration");
        }
        const To *emax_pos = &transformed[1];
        const auto *coeff_pos = reinterpret_cast<const typename ops::uint_type *>(&transformed[1 + n_blocks]);
        for_each_block(conf, [&](ptrdiff_t off, const ptrdiff_t *s, const size_t *n) {
            ops::inv(dec_data + off, s, n, static_cast<int>(*emax_pos++), coeff_pos);
            coeff_pos += bs;
        });
        return dec_data;
    }

    void save(uchar *&c) override {}

    void load(const uchar *&c, size_t &remaining_length) override {}

    std::pair<To, To> get_out_range() override { return std::make_pair(0, 0); }
};

template <class T, class To, uint N>
ZFPDecomposition<T, To, N> make_decomposition_zfp() {
    return ZFPDecomposition<T, To, N>();
}

}  // namespace SZ3
#endif
