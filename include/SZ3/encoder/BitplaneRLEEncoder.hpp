/**
 * @file BitplaneRLEEncoder.hpp
 * @ingroup Encoder
 * @brief Bit-plane encoder with per-plane run-length coding and a raw-pack
 *        fallback.
 *
 * For each MSB→LSB plane the encoder computes both the RLE cost (alternating
 * 0/1 run lengths as `uint32_t`) and the raw bit-pack cost (8 bins per byte),
 * then writes a 1-byte mode flag (0 = RAW, 1 = RLE) followed by the smaller
 * payload. Sparse planes (typical of high-magnitude bits) compress strongly
 * with RLE; dense planes (typical of low-magnitude bits) fall back to raw
 * pack so they don't pay an RLE overhead penalty.
 *
 * Sign vector (only emitted when the input contains negatives) is always
 * raw-packed.
 *
 * Drop-in alternative to `BitplaneEncoder<T>` with the same interface.
 */

#ifndef SZ3_BITPLANE_RLE_ENCODER_HPP
#define SZ3_BITPLANE_RLE_ENCODER_HPP

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

template <class T>
class BitplaneRLEEncoder : public concepts::EncoderInterface<T> {
    static_assert(std::is_integral<T>::value, "BitplaneRLEEncoder only supports integer types.");

    using UT = typename std::make_unsigned<T>::type;
    enum PlaneMode : uint8_t { MODE_RAW = 0, MODE_RLE = 1 };

   public:
    void preprocess_encode(const std::vector<T>& bins, int /*stateNum*/) override {
        n_ = bins.size();
        if (n_ == 0) {
            num_planes_ = 0;
            has_signs_ = false;
            return;
        }
        UT max_abs = 0;
        bool has_neg = false;
        for (size_t i = 0; i < n_; i++) {
            const T v = bins[i];
            UT a;
            if (v < 0) {
                has_neg = true;
                a = static_cast<UT>(0) - static_cast<UT>(v);
            } else {
                a = static_cast<UT>(v);
            }
            if (a > max_abs) max_abs = a;
        }
        has_signs_ = has_neg;
        num_planes_ = 0;
        for (UT m = max_abs; m > 0; m >>= 1) num_planes_++;
        const int max_p = static_cast<int>(8 * sizeof(T)) - (std::is_signed<T>::value ? 2 : 1);
        if (num_planes_ > max_p + 1) num_planes_ = max_p + 1;
    }

    size_t encode(const std::vector<T>& bins, uchar*& bytes) override {
        if (bins.size() != n_) {
            throw std::runtime_error("BitplaneRLEEncoder: encode size mismatch with preprocess_encode.");
        }
        uchar* const start = bytes;
        if (n_ == 0) return 0;

        const size_t plane_bytes = (n_ + 7) / 8;

        // Build magnitudes once so plane sweeps are pure ints.
        std::vector<UT> mags(n_);
        for (size_t i = 0; i < n_; i++) {
            const T v = bins[i];
            mags[i] = (v < 0) ? static_cast<UT>(static_cast<UT>(0) - static_cast<UT>(v))
                              : static_cast<UT>(v);
        }

        if (has_signs_) {
            std::memset(bytes, 0, plane_bytes);
            for (size_t i = 0; i < n_; i++) {
                if (bins[i] < 0) bytes[i >> 3] |= static_cast<uchar>(1u << (i & 7));
            }
            bytes += plane_bytes;
        }

        std::vector<uint32_t> runs;
        runs.reserve(64);

        const bool log_planes = std::getenv("BITPLANE_LOG") != nullptr;
        for (int p = num_planes_ - 1; p >= 0; p--) {
            const UT mask = static_cast<UT>(static_cast<UT>(1) << p);
            uchar* plane_start = bytes;

            // Build RLE runs (alternating false/true counts; if first bit is true,
            // a leading 0 placeholder represents the empty initial false run).
            runs.clear();
            bool last_bit = (mags[0] & mask) != 0;
            if (last_bit) runs.push_back(0);
            uint32_t count = 1;
            for (size_t i = 1; i < n_; i++) {
                const bool b = (mags[i] & mask) != 0;
                if (b == last_bit) {
                    count++;
                } else {
                    runs.push_back(count);
                    count = 1;
                    last_bit = b;
                }
            }
            runs.push_back(count);

            // Pick the cheaper representation. RLE pays 1 (mode) + 4 (num_runs) + 4*runs.
            const size_t rle_bytes = 1 + sizeof(uint32_t) + runs.size() * sizeof(uint32_t);
            const size_t raw_bytes = 1 + plane_bytes;

            if (rle_bytes < raw_bytes) {
                *bytes++ = static_cast<uchar>(MODE_RLE);
                const uint32_t num_runs = static_cast<uint32_t>(runs.size());
                std::memcpy(bytes, &num_runs, sizeof(num_runs));
                bytes += sizeof(num_runs);
                std::memcpy(bytes, runs.data(), runs.size() * sizeof(uint32_t));
                bytes += runs.size() * sizeof(uint32_t);
            } else {
                *bytes++ = static_cast<uchar>(MODE_RAW);
                std::memset(bytes, 0, plane_bytes);
                for (size_t i = 0; i < n_; i++) {
                    if (mags[i] & mask) bytes[i >> 3] |= static_cast<uchar>(1u << (i & 7));
                }
                bytes += plane_bytes;
            }
            if (log_planes) {
                size_t plane_size = static_cast<size_t>(bytes - plane_start);
                const char* mode_str = (rle_bytes < raw_bytes) ? "RLE" : "RAW";
                fprintf(stderr, "BITPLANE_LOG plane=%2d mode=%s bytes=%zu runs=%zu\n",
                        p, mode_str, plane_size, runs.size());
            }
        }
        return static_cast<size_t>(bytes - start);
    }

    std::vector<T> decode(const uchar*& bytes, size_t targetLength) override {
        if (targetLength != n_) {
            throw std::runtime_error("BitplaneRLEEncoder: decode targetLength does not match saved bin count.");
        }
        std::vector<T> bins(n_, 0);
        if (n_ == 0) return bins;
        const size_t plane_bytes = (n_ + 7) / 8;

        std::vector<uint8_t> signs;
        if (has_signs_) {
            signs.assign(n_, 0);
            for (size_t i = 0; i < n_; i++) {
                signs[i] = static_cast<uint8_t>((bytes[i >> 3] >> (i & 7)) & 1u);
            }
            bytes += plane_bytes;
        }

        std::vector<UT> mags(n_, 0);
        for (int p = num_planes_ - 1; p >= 0; p--) {
            const uchar mode = *bytes++;
            const UT bit_val = static_cast<UT>(static_cast<UT>(1) << p);
            if (mode == MODE_RLE) {
                uint32_t num_runs = 0;
                std::memcpy(&num_runs, bytes, sizeof(num_runs));
                bytes += sizeof(num_runs);
                bool current_bit = false;
                size_t pos = 0;
                for (uint32_t r = 0; r < num_runs; r++) {
                    uint32_t len;
                    std::memcpy(&len, bytes, sizeof(len));
                    bytes += sizeof(len);
                    if (current_bit) {
                        for (uint32_t k = 0; k < len; k++) {
                            if (pos < n_) mags[pos] |= bit_val;
                            pos++;
                        }
                    } else {
                        pos += len;
                    }
                    current_bit = !current_bit;
                }
                if (pos != n_) {
                    throw std::runtime_error("BitplaneRLEEncoder: RLE run-length total != n.");
                }
            } else if (mode == MODE_RAW) {
                for (size_t i = 0; i < n_; i++) {
                    if ((bytes[i >> 3] >> (i & 7)) & 1u) mags[i] |= bit_val;
                }
                bytes += plane_bytes;
            } else {
                throw std::runtime_error("BitplaneRLEEncoder: unknown plane mode flag.");
            }
        }

        for (size_t i = 0; i < n_; i++) {
            if (has_signs_ && signs[i]) {
                bins[i] = static_cast<T>(static_cast<UT>(0) - mags[i]);
            } else {
                bins[i] = static_cast<T>(mags[i]);
            }
        }
        return bins;
    }

    void save(uchar*& c) override {
        write(n_, c);
        write(num_planes_, c);
        const uint8_t flags = has_signs_ ? 1 : 0;
        write(flags, c);
    }

    void load(const uchar*& c, size_t& remaining_length) override {
        read(n_, c, remaining_length);
        read(num_planes_, c, remaining_length);
        uint8_t flags = 0;
        read(flags, c, remaining_length);
        has_signs_ = (flags & 1u) != 0;
    }

    void preprocess_decode() override {}
    void postprocess_encode() override {}
    void postprocess_decode() override {}

    // Worst-case bound: assume every plane chose RAW (RLE never used).
    // Covers SZGenericCompressor's buffer sizing without overruns.
    size_t size_est() override {
        const size_t header = sizeof(n_) + sizeof(num_planes_) + sizeof(uint8_t);
        if (n_ == 0) return header;
        const size_t plane_bytes = (n_ + 7) / 8;
        const size_t per_plane_max = 1 /*mode*/ + plane_bytes;  // RAW upper bound
        const size_t total_planes = static_cast<size_t>(num_planes_) + (has_signs_ ? 1 : 0);
        return header + total_planes * per_plane_max;
    }

   private:
    size_t n_ = 0;
    int num_planes_ = 0;
    bool has_signs_ = false;
};

}  // namespace SZ3

#endif
