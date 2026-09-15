/**
 * @file ZFPEncoder.hpp
 * @ingroup Encoder
 */

#ifndef SZ3_ZFP_ENCODER_HPP
#define SZ3_ZFP_ENCODER_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Encoder.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/thirdparty/zfp/zfp_codec.hpp"

namespace SZ3 {

/**
 * @brief zfp CODEC 5's embedded coder: writes each block's exponent and then its
 *        coefficients under the group-testing bitplane coder.
 *
 * This is where zfp's error control lives. The precision granted to a block follows from
 * that block's exponent and from the absolute error bound handed to the constructor;
 * ZFPDecomposition upstream of it carries no bound at all.
 *
 * Consumes the layout ZFPDecomposition produces and nothing else. Pair it with
 * Lossless_bypass -- the output is already packed.
 */
template <class Int, uint N>
class ZFPEncoder : public concepts::EncoderInterface<Int> {
    static_assert(N >= 1 && N <= 3, "zfp modules are instantiated for 1D, 2D and 3D");
    static_assert(sizeof(Int) == 4 || sizeof(Int) == 8, "bins are the integer of the scalar's width");

    /// The bin width identifies the scalar the coefficients came from.
    using T = std::conditional_t<sizeof(Int) == 4, float, double>;
    using ops = ZFP::ops<T, N>;
    static constexpr int bs = ops::block_size;

   public:
    explicit ZFPEncoder(const Config &conf) : minexp(ZFP::minexp_for(conf.absErrorBound)) {}

    size_t encode(const std::vector<Int> &data, uchar *&bytes) override {
        if (data.empty()) return 0;
        const size_t n_blocks = static_cast<size_t>(data[0]);
        if (data.size() != 1 + n_blocks + n_blocks * bs) {
            throw std::out_of_range("SZ3 zfp: coefficient stream length does not match its block count");
        }
        // zfp's bitstream does no bounds checking. Record the block count and the coded size so
        // decode() can cross-check them against the length its caller passes.
        uchar *const start = bytes;
        uchar *payload = start + header_bytes;
        ZFP::bitstream *s = ZFP::stream_open(payload, capacity_for(n_blocks));
        const Int *emax_pos = &data[1];
        const auto *coeff_pos = reinterpret_cast<const typename ops::uint_type *>(&data[1 + n_blocks]);
        for (size_t i = 0; i < n_blocks; i++) {
            ops::enc(s, static_cast<int>(*emax_pos++), coeff_pos, minexp);
            coeff_pos += bs;
        }
        ZFP::stream_flush(s);
        const size_t coded = ZFP::stream_size(s);
        ZFP::stream_close(s);
        uchar *hpos = start;
        write(static_cast<uint64_t>(n_blocks), hpos);
        write(static_cast<uint64_t>(coded), hpos);
        const size_t total = header_bytes + coded;
        bytes += total;
        return total;
    }

    std::vector<Int> decode(const uchar *&bytes, size_t targetLength) override {
        // targetLength counts ZFPDecomposition's layout: a block count, then an exponent and
        // `bs` coefficients per block. Anything else did not come from encode().
        if (targetLength < 1 || (targetLength - 1) % (1 + bs) != 0) {
            throw std::out_of_range("SZ3 zfp: coefficient count is not a whole number of blocks");
        }
        const size_t n_blocks = (targetLength - 1) / (1 + bs);
        const uchar *hpos = bytes;
        uint64_t recorded_blocks = 0, coded = 0;
        read(recorded_blocks, hpos);
        read(coded, hpos);
        // The stream says how many blocks it holds; the caller says how many it expects. They are
        // written by different layers, so requiring them to agree rejects a stream whose count was
        // edited in one place.
        if (recorded_blocks != n_blocks) {
            throw std::out_of_range("SZ3 zfp: block count disagrees with the coefficient count");
        }
        if (coded > capacity_for(n_blocks)) {
            throw std::out_of_range("SZ3 zfp: coded length exceeds what this block count can produce");
        }
        std::vector<Int> out(1 + n_blocks + n_blocks * bs);
        out[0] = static_cast<Int>(n_blocks);
        // zfp's bitstream does not stop at the end it is handed, and a corrupted block asks for
        // more bits than the encoder wrote -- one flipped byte is enough. Read from a padded copy,
        // so a run-on stays inside this buffer, and reject the stream once the reader passes the
        // length the encoder recorded.
        std::vector<uchar> payload(coded + block_bytes + sizeof(uint64_t), 0);
        std::memcpy(payload.data(), bytes + header_bytes, coded);
        ZFP::bitstream *s = ZFP::stream_open(payload.data(), payload.size());
        Int *emax_pos = &out[1];
        auto *coeff_pos = reinterpret_cast<typename ops::uint_type *>(&out[1 + n_blocks]);
        for (size_t i = 0; i < n_blocks; i++) {
            int emax = 0;
            ops::dec(s, &emax, coeff_pos, minexp);
            if (ZFP::stream_rtell(s) > 8 * coded) {
                ZFP::stream_close(s);
                throw std::out_of_range("SZ3 zfp: a block reads past the coded length");
            }
            *emax_pos++ = static_cast<Int>(emax);
            coeff_pos += bs;
        }
        ZFP::stream_close(s);
        bytes += header_bytes + coded;
        return out;
    }

    /// Worst case for the variable-rate path: every coefficient at full precision, plus the
    /// exponent and the nonzero-block flag. Recorded by preprocess_encode().
    size_t size_est() override { return capacity_for(n_blocks_hint); }

    void preprocess_encode(const std::vector<Int> &bins, int stateNum) override {
        n_blocks_hint = bins.empty() ? 0 : static_cast<size_t>(bins[0]);
    }
    void postprocess_encode() override {}
    void preprocess_decode() override {}
    void postprocess_decode() override {}
    /// Nothing to serialise: the bound reaches decode() through the Config, which the API has
    /// already resolved to EB_ABS by the time it is written.
    void save(uchar *&c) override {}
    void load(const uchar *&c, size_t &remaining_length) override {}

   private:
    /// Block count and coded length, both little-endian like the rest of the stream.
    static constexpr size_t header_bytes = 2 * sizeof(uint64_t);

    /// Most a block can occupy: 1 flag + exponent + one sign-and-value bit per coefficient bitplane.
    static constexpr size_t block_bytes = (1 + 8 * sizeof(T) + bs * (8 * sizeof(Int) + 1)) / 8 + 1;

    static size_t capacity_for(size_t n_blocks) { return header_bytes + n_blocks * block_bytes + 64; }

    int minexp;
    size_t n_blocks_hint = 0;
};

}  // namespace SZ3
#endif
