/**
 * @file PaSTRIDecomposition.hpp
 * @ingroup Decomposition
 */

#ifndef SZ3_PASTRI_DECOMPOSITION_HPP
#define SZ3_PASTRI_DECOMPOSITION_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * @brief PaSTRI: pattern-scaling decomposition for two-electron repulsion integrals (ERI).
 *
 * PaSTRI is an error-bounded lossy decomposition designed for the two-electron repulsion
 * integrals produced by quantum chemistry codes such as GAMESS. Reference:
 * A. M. Gok et al., "PaSTRI: Error-Bounded Lossy Compression for Two-Electron Integrals in
 * Quantum Chemistry", IEEE Cluster 2018.
 *
 * ERI data is laid out as a sequence of equally sized blocks; within a block the sub-blocks are
 * near-scalar-multiples of one another. PaSTRI exploits exactly that: for each block it picks one
 * sub-block as a *pattern*, derives one *scale* per sub-block, and stores a quantized
 * *error correction* (EC) term per point.
 *
 * Per block, with @c binSize = 2*eb :
 *   1. Find the global extremum @c data[extIdx] of the block (largest magnitude).
 *   2. The pattern is the sub-block that contains it: @c patternIdx = (extIdx/sbSize)*sbSize .
 *      This is the paper's "ratio of extremums" pattern metric.
 *   3. @c patternQ[i] = round(data[patternIdx+i] / binSize) .
 *   4. @c patternBits = bitsNeeded(|patternExt|/binSize + 1) + 1 , @c scaleBits = patternBits ,
 *      @c scalesBinSize = 1 / (2^(scaleBits-1) - 1) . Because the pattern is taken from the
 *      extremum's sub-block, every scale is guaranteed to fall in [-1,1], so the whole integer
 *      range can be spent on that interval.
 *   5. One scale per sub-block, from the single point at the same local offset as the global
 *      extremum (@c localExtIdx = extIdx % sbSize):
 *      @c scalesQ[sb] = round((data[sb*sbSize+localExtIdx] / patternExt) / scalesBinSize) .
 *   6. @c ECQ[i] = round((scalesQ[sb]*patternQ[i]*scalesBinSize*binSize - data[i]) / binSize) .
 *
 * Reconstruction is @c scalesQ[sb]*patternQ[i]*scalesBinSize*binSize - ECQ[i]*binSize . Note the
 * sign convention: ECQ is *predicted minus actual*, so reconstruction subtracts it.
 *
 * **The caller must supply the basis-function configuration.** The sub-block geometry comes from
 * the four orbital / basis-function types of the integral shell quartet, not from the data or from
 * @c Config::dims :
 * @code
 *   idxRange[i] = (bf[i]+1)*(bf[i]+2)/2       for i in 0..3
 *   sbSize      = idxRange[2]*idxRange[3]
 *   sbNum       = idxRange[0]*idxRange[1]
 *   bSize       = sbSize*sbNum
 * @endcode
 * They are passed to the constructor (typical values are in [0,3]) and serialized by @c save().
 *
 * **Tail policy.** The paper only forms patterns on full-sized blocks. When @c conf.num is not an
 * exact multiple of @c bSize, this module compresses the leading @c floor(conf.num/bSize) full
 * blocks with the PaSTRI pattern/scale model and falls back to plain quantization against a zero
 * prediction (@c pred = 0) for the trailing @c conf.num % bSize values, which also covers
 * @c conf.num < bSize. The tail is typically incompressible; for real work supply data whose
 * length is a multiple of @c bSize.
 *
 * **Output type.** The emitted bins are radius-shifted EC terms in [0, 2*radius], so @c int is
 * ample (radius defaults to 32768). The wide intermediates -- @c patternQ and @c scalesQ, which
 * can need up to ~62 bits -- are @c int64_t but live in the @c save() metadata, never in the bin
 * stream. Bin 0 is reserved for "unpredictable": a point whose EC term falls outside the radius,
 * or whose reconstruction would violate the error bound, is stored verbatim instead, which makes
 * the error bound unconditional.
 *
 * @c compress() does not modify the input array.
 *
 * @tparam T Original data type (float or double)
 * @tparam N Original data dimension (unused; PaSTRI treats the input as a flat block sequence)
 */
template <class T, uint N>
class PaSTRIDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    /**
     * @brief Construct a PaSTRI decomposition
     *
     * @param eb Absolute error bound (must be > 0)
     * @param bf The four basis-function (orbital) types of the shell quartet, typically in [0,3]
     * @param radius Quantization radius; bins are in [0, 2*radius] and bin 0 means "unpredictable"
     */
    PaSTRIDecomposition(double eb, const std::array<int, 4> &bf, int radius = 32768)
        : error_bound(eb), bin_size(2 * eb), bf(bf), radius(radius) {
        if (!(eb > 0)) {
            throw std::invalid_argument("PaSTRIDecomposition: error bound must be positive");
        }
        if (radius < 2) {
            throw std::invalid_argument("PaSTRIDecomposition: radius must be at least 2");
        }
        compute_geometry();
    }

    ~PaSTRIDecomposition() override = default;

    std::vector<int> compress(const Config &conf, T *data) override {
        unpred.clear();
        unpred_index = 0;

        num_full_blocks = conf.num / block_size;
        tail_len = conf.num - num_full_blocks * block_size;

        pattern_q.assign(num_full_blocks * sb_size, 0);
        scales_q.assign(num_full_blocks * sb_num, 0);
        scale_bits.assign(num_full_blocks, 0);

        std::vector<int> quant_inds(conf.num);

        for (size_t b = 0; b < num_full_blocks; b++) {
            const T *blk = data + b * block_size;
            int *bins = quant_inds.data() + b * block_size;
            int64_t *pat = pattern_q.data() + b * sb_size;
            int64_t *sca = scales_q.data() + b * sb_num;

            // 1. Global extremum of the block.
            double abs_ext = 0;
            size_t ext_idx = 0;
            for (size_t i = 0; i < block_size; i++) {
                double a = std::fabs(static_cast<double>(blk[i]));
                if (a > abs_ext) {
                    abs_ext = a;
                    ext_idx = i;
                }
            }
            // Note: the reference C implementation seeds ext_idx with -1 and reads data[-1] on an
            // all-zero block. Seeding with 0 keeps that case in bounds and yields patternExt == 0,
            // which the scale loop below already guards.

            // 2. The pattern is the sub-block holding the extremum.
            const size_t pattern_idx = (ext_idx / sb_size) * sb_size;
            const double pattern_ext = static_cast<double>(blk[ext_idx]);

            // 3. Quantize the pattern.
            for (size_t i = 0; i < sb_size; i++) {
                pat[i] = round_half_away(static_cast<double>(blk[pattern_idx + i]) / bin_size);
            }

            // 4. Bit widths. abs_ext/bin_size + 1 >= 1, so bits_needed() >= 1 and sb >= 2.
            const int pattern_bits = bits_needed(abs_ext / bin_size + 1) + 1;
            const int sb_bits = std::min(std::max(pattern_bits, 2), kMaxScaleBits);
            scale_bits[b] = static_cast<uint8_t>(sb_bits);
            const double scales_bin_size = scales_bin_size_of(sb_bits);

            // 5. One scale per sub-block, taken from the local offset of the global extremum.
            const size_t local_ext_idx = ext_idx % sb_size;
            const bool pattern_ext_zero = (pattern_ext == 0);
            for (size_t s = 0; s < sb_num; s++) {
                const double scale =
                    pattern_ext_zero ? 0.0 : static_cast<double>(blk[s * sb_size + local_ext_idx]) / pattern_ext;
                sca[s] = round_half_away(scale / scales_bin_size);
            }

            // 6. Per-point error correction.
            const double ps_bin_size = scales_bin_size * bin_size;
            for (size_t s = 0; s < sb_num; s++) {
                for (size_t i = 0; i < sb_size; i++) {
                    const size_t idx = s * sb_size + i;
                    bins[idx] = quantize_ec(blk[idx], predict(sca[s], pat[i], ps_bin_size));
                }
            }
        }

        // Tail: no pattern is formed, quantize against a zero prediction.
        for (size_t i = num_full_blocks * block_size; i < conf.num; i++) {
            quant_inds[i] = quantize_ec(data[i], 0.0);
        }

        return quant_inds;
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        const size_t expected_blocks = conf.num / block_size;
        if (expected_blocks != num_full_blocks || quant_inds.size() < conf.num) {
            throw std::runtime_error("PaSTRIDecomposition::decompress: input size does not match the saved state");
        }
        unpred_index = 0;

        for (size_t b = 0; b < num_full_blocks; b++) {
            T *out = dec_data + b * block_size;
            const int *bins = quant_inds.data() + b * block_size;
            const int64_t *pat = pattern_q.data() + b * sb_size;
            const int64_t *sca = scales_q.data() + b * sb_num;

            const double ps_bin_size = scales_bin_size_of(scale_bits[b]) * bin_size;
            for (size_t s = 0; s < sb_num; s++) {
                for (size_t i = 0; i < sb_size; i++) {
                    const size_t idx = s * sb_size + i;
                    out[idx] = recover_ec(predict(sca[s], pat[i], ps_bin_size), bins[idx]);
                }
            }
        }

        for (size_t i = num_full_blocks * block_size; i < conf.num; i++) {
            dec_data[i] = recover_ec(0.0, quant_inds[i]);
        }
        return dec_data;
    }

    void save(uchar *&c) override {
        write(kUid, c);
        write(bf.data(), bf.size(), c);
        write(error_bound, c);
        write(radius, c);
        write(static_cast<uint64_t>(num_full_blocks), c);
        write(static_cast<uint64_t>(tail_len), c);
        if (num_full_blocks > 0) {
            write(scale_bits.data(), scale_bits.size(), c);
        }
        write_varints(pattern_q, c);
        write_varints(scales_q, c);
        write(static_cast<uint64_t>(unpred.size()), c);
        if (!unpred.empty()) {
            write(unpred.data(), unpred.size(), c);
        }
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        uchar uid_read = 0;
        read(uid_read, c, remaining_length);
        if (uid_read != kUid) {
            throw std::invalid_argument("PaSTRIDecomposition uid mismatch");
        }
        read(bf.data(), bf.size(), c, remaining_length);
        read(error_bound, c, remaining_length);
        read(radius, c, remaining_length);
        if (!(error_bound > 0) || radius < 2) {
            throw std::runtime_error("PaSTRIDecomposition::load: corrupted header");
        }
        bin_size = 2 * error_bound;
        compute_geometry();

        uint64_t nb = 0, tl = 0;
        read(nb, c, remaining_length);
        read(tl, c, remaining_length);
        num_full_blocks = static_cast<size_t>(nb);
        tail_len = static_cast<size_t>(tl);

        require(remaining_length >= num_full_blocks, "truncated scale bit widths");
        scale_bits.resize(num_full_blocks);
        if (num_full_blocks > 0) {
            read(scale_bits.data(), scale_bits.size(), c, remaining_length);
            for (uint8_t sb : scale_bits) {
                if (sb < 2 || sb > kMaxScaleBits) {
                    throw std::runtime_error("PaSTRIDecomposition::load: corrupted scale bit width");
                }
            }
        }

        pattern_q.resize(num_full_blocks * sb_size);
        scales_q.resize(num_full_blocks * sb_num);
        read_varints(pattern_q, c, remaining_length);
        read_varints(scales_q, c, remaining_length);

        uint64_t unpred_size = 0;
        read(unpred_size, c, remaining_length);
        require(unpred_size <= remaining_length / sizeof(T), "truncated unpredictable value list");
        unpred.resize(static_cast<size_t>(unpred_size));
        if (unpred_size > 0) {
            read(unpred.data(), unpred.size(), c, remaining_length);
        }
        unpred_index = 0;
    }

    size_t size_est() override {
        // header + per-block scale width + worst-case 10-byte varints + unpredictable values
        return 64 + num_full_blocks * (1 + kMaxVarintBytes * (sb_size + sb_num)) + unpred.size() * sizeof(T);
    }

    std::pair<int, int> get_out_range() override { return std::make_pair(0, radius * 2); }

    void print() override {
        printf("[PaSTRIDecomposition] eb = %.8G, bf = {%d,%d,%d,%d}, sbSize = %zu, sbNum = %zu, bSize = %zu\n",
               error_bound, bf[0], bf[1], bf[2], bf[3], sb_size, sb_num, block_size);
        printf("                      blocks = %zu, tail = %zu, radius = %d, unpred = %zu\n", num_full_blocks, tail_len,
               radius, unpred.size());
    }

    /// Number of values in one PaSTRI block (sbSize * sbNum).
    size_t get_block_size() const { return block_size; }

    /// Number of values in one sub-block (idxRange[2] * idxRange[3]).
    size_t get_sub_block_size() const { return sb_size; }

    /// Number of sub-blocks in one block (idxRange[0] * idxRange[1]).
    size_t get_sub_block_num() const { return sb_num; }

   private:
    static constexpr uchar kUid = 0xA7;
    /// Cap on scaleBits so that (1ULL << (scaleBits-1)) - 1 stays well defined.
    static constexpr int kMaxScaleBits = 62;
    static constexpr size_t kMaxVarintBytes = 10;

    void compute_geometry() {
        std::array<size_t, 4> idx_range{};
        for (size_t i = 0; i < 4; i++) {
            if (bf[i] < 0 || bf[i] > 30) {
                throw std::invalid_argument("PaSTRIDecomposition: basis function type must be in [0,30]");
            }
            const size_t v = static_cast<size_t>(bf[i]);
            idx_range[i] = (v + 1) * (v + 2) / 2;
        }
        sb_size = idx_range[2] * idx_range[3];
        sb_num = idx_range[0] * idx_range[1];
        block_size = sb_size * sb_num;
    }

    /// Round half away from zero, matching the reference implementation's quantizer.
    static int64_t round_half_away(double x) { return static_cast<int64_t>(std::llround(x)); }

    /**
     * @brief Minimum number of bits needed to represent x, i.e. floor(log2(x)) + 1 for x >= 1.
     *
     * Matches @c bitsNeeded_double() in the reference implementation over its used domain (x >= 1).
     * Returns 0 for x < 1, where the reference's exponent arithmetic would go negative.
     */
    static int bits_needed(double x) {
        if (!(x >= 1.0)) {
            return 0;
        }
        int exponent = 0;
        std::frexp(x, &exponent);  // x = m * 2^exponent with m in [0.5, 1)
        return exponent;
    }

    static double scales_bin_size_of(int sb_bits) {
        return 1.0 / static_cast<double>((static_cast<uint64_t>(1) << (sb_bits - 1)) - 1);
    }

    /**
     * @brief The PaSTRI prediction for one point: scale * pattern, in physical units.
     *
     * Computed in double rather than as an int64 product; @c scalesQ and @c patternQ can each need
     * up to ~62 bits, so the integer product in the reference implementation can overflow. Both
     * @c compress() and @c decompress() go through this one function, so the two agree bit for bit.
     */
    static double predict(int64_t scale_q, int64_t pattern_q_i, double ps_bin_size) {
        return static_cast<double>(scale_q) * static_cast<double>(pattern_q_i) * ps_bin_size;
    }

    /**
     * @brief The single reconstruction expression, shared by compress and decompress.
     *
     * Keeping it in one place means the two directions cannot drift by an ULP through differing
     * floating-point contraction of @c pred - ecq*binSize into an FMA.
     */
    double reconstruct(double pred, int64_t ecq) const { return pred - static_cast<double>(ecq) * bin_size; }

    static void require(bool ok, const char *what) {
        if (!ok) {
            throw std::runtime_error(std::string("PaSTRIDecomposition::load: ") + what);
        }
    }

    /**
     * @brief Quantize one error-correction term and shift it into [0, 2*radius].
     *
     * ECQ follows the paper's sign convention (predicted minus actual). The returned bin is
     * @c radius - ECQ, so bin 0 stays free to mean "unpredictable"; out-of-range terms and any
     * point whose reconstruction would exceed the error bound are pushed to @c unpred instead,
     * which makes the bound unconditional.
     */
    int quantize_ec(T value, double pred) {
        const double diff = pred - static_cast<double>(value);
        const double ratio = diff / bin_size;
        // The +1 margin keeps |ECQ| < radius after rounding, so the shifted bin never hits 0.
        // NaN and infinity fail this test and fall through to the unpredictable path.
        if (std::fabs(ratio) + 1 < static_cast<double>(radius)) {
            const int64_t ecq = round_half_away(ratio);
            const T rec = static_cast<T>(reconstruct(pred, ecq));
            if (std::fabs(static_cast<double>(rec) - static_cast<double>(value)) <= error_bound) {
                return static_cast<int>(static_cast<int64_t>(radius) - ecq);
            }
        }
        unpred.push_back(value);
        return 0;
    }

    T recover_ec(double pred, int bin) {
        if (bin == 0) {
            if (unpred_index >= unpred.size()) {
                throw std::runtime_error("PaSTRIDecomposition: unpredictable value list exhausted");
            }
            return unpred[unpred_index++];
        }
        const int64_t ecq = static_cast<int64_t>(radius) - static_cast<int64_t>(bin);
        return static_cast<T>(reconstruct(pred, ecq));
    }

    /// Zigzag + LEB128 varint. The pattern and scale magnitudes are typically far below 2^62.
    static void write_varints(const std::vector<int64_t> &v, uchar *&c) {
        for (int64_t x : v) {
            uint64_t u = (static_cast<uint64_t>(x) << 1) ^ static_cast<uint64_t>(x >> 63);
            while (u >= 0x80) {
                *c++ = static_cast<uchar>((u & 0x7F) | 0x80);
                u >>= 7;
            }
            *c++ = static_cast<uchar>(u);
        }
    }

    static void read_varints(std::vector<int64_t> &v, const uchar *&c, size_t &remaining_length) {
        for (auto &x : v) {
            uint64_t u = 0;
            int shift = 0;
            while (true) {
                if (remaining_length == 0) {
                    throw std::runtime_error("PaSTRIDecomposition::load: buffer underflow");
                }
                const uchar byte = *c++;
                remaining_length--;
                u |= static_cast<uint64_t>(byte & 0x7F) << shift;
                if ((byte & 0x80) == 0) {
                    break;
                }
                shift += 7;
                if (shift > 63) {
                    throw std::runtime_error("PaSTRIDecomposition::load: malformed varint");
                }
            }
            x = static_cast<int64_t>((u >> 1) ^ (~(u & 1) + 1));
        }
    }

    double error_bound;
    double bin_size;  // 2 * error_bound
    std::array<int, 4> bf;
    int radius;

    // Geometry derived from bf.
    size_t sb_size = 0;
    size_t sb_num = 0;
    size_t block_size = 0;

    // Per-block model, serialized by save().
    std::vector<int64_t> pattern_q;  // num_full_blocks * sb_size
    std::vector<int64_t> scales_q;   // num_full_blocks * sb_num
    std::vector<uint8_t> scale_bits;
    size_t num_full_blocks = 0;
    size_t tail_len = 0;

    std::vector<T> unpred;
    size_t unpred_index = 0;  // decompression only
};

/**
 * @brief Factory for PaSTRIDecomposition, taking the error bound from the Config
 *
 * @param conf Compression configuration (uses @c conf.absErrorBound)
 * @param bf The four basis-function (orbital) types of the shell quartet
 * @param radius Quantization radius
 */
template <class T, uint N>
PaSTRIDecomposition<T, N> make_decomposition_pastri(const Config &conf, const std::array<int, 4> &bf,
                                                    int radius = 32768) {
    return PaSTRIDecomposition<T, N>(conf.absErrorBound, bf, radius);
}

}  // namespace SZ3

#endif
