#ifndef SZ3_SPERR_ENCODER_HPP
#define SZ3_SPERR_ENCODER_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/thirdparty/sperr/SPERRHeaderOnly.hpp"

namespace SZ3 {

/**
 * @file SPERREncoder.hpp
 * @brief SPECK-based encoder for integer bins.
 */
template <class T, uint N>
class SPERREncoder : public concepts::EncoderInterface<T> {
   public:
    /**
     * @brief Encoder module for integer bins using SPERR SPECK coding.
     *
     * Designed for compressor pipelines that already produce integer bins and
     * need a compact bitstream stage through `EncoderInterface`.
     * `N` selects the coding context (1D/2D/3D).
     */
    static_assert(std::is_integral<T>::value, "SPERREncoder scalar type must be integral.");
    static_assert(N >= 1 && N <= 3, "SPERREncoder only supports N in [1, 3].");

    SPERREncoder() : dims_({0, 0, 0}) {}

    explicit SPERREncoder(const SZ3::SPERR::dims_type &dims) : dims_(dims) {}

    void set_dims(const SZ3::SPERR::dims_type &dims) { dims_ = dims; }

    const SZ3::SPERR::dims_type &get_dims() const { return dims_; }

    void preprocess_encode(const std::vector<T> &bins, int /*stateNum*/) override {
        total_len_ = bins.size();
        (void)resolved_dims(total_len_);
    }

    size_t encode(const std::vector<T> &bins, uchar *&bytes) override {
        if (bins.empty()) {
            throw std::runtime_error("SPERR input bins cannot be empty.");
        }

        const size_t total_len = bins.size();
        if (total_len_ == 0) {
            total_len_ = total_len;
        } else if (total_len_ != total_len) {
            throw std::runtime_error("SPERR bin length mismatch between preprocess and encode.");
        }

        std::vector<uchar> stream = encode_stream(bins);
        std::memcpy(bytes, stream.data(), stream.size());
        bytes += stream.size();
        return stream.size();
    }

    std::vector<T> decode(const uchar *&bytes, size_t targetLength) override {
        const size_t total_len = (total_len_ == 0) ? targetLength : total_len_;
        if (total_len == 0) {
            throw std::runtime_error("SPERR encoder missing target length before decode.");
        }

        size_t consumed = 0;
        std::vector<T> bins = decode_stream(bytes, std::numeric_limits<size_t>::max(), total_len, consumed);
        bytes += consumed;

        if (targetLength != 0 && bins.size() != targetLength) {
            throw std::runtime_error("SPERR decode length mismatch.");
        }
        return bins;
    }

    std::vector<uchar> encode_stream(const std::vector<T> &bins) const {
        if (bins.empty()) {
            throw std::runtime_error("SPERR input bins cannot be empty.");
        }

        const std::vector<int64_t> signed_vals = to_signed_vector(bins);
        const SZ3::SPERR::dims_type dims = resolved_dims(signed_vals.size());

        if (N == 1) {
            return encode_speck_1d(signed_vals);
        }
        if (N == 2) {
            return encode_speck_2d(signed_vals, dims);
        }
        return encode_speck_3d(signed_vals, dims);
    }

    std::vector<T> decode_stream(const uchar *cmpData, size_t cmpSize, size_t total_len, size_t &consumed) const {
        if (total_len == 0) {
            throw std::runtime_error("SPERR decode requires a positive target length.");
        }

        const SZ3::SPERR::dims_type dims = resolved_dims(total_len);
        std::vector<int64_t> signed_vals;

        if (N == 1) {
            signed_vals = decode_speck_1d(total_len, cmpData, cmpSize, consumed);
        } else if (N == 2) {
            signed_vals = decode_speck_2d(dims, cmpData, cmpSize, consumed);
        } else {
            signed_vals = decode_speck_3d(dims, cmpData, cmpSize, consumed);
        }

        return from_signed_vector(signed_vals);
    }

    void save(uchar *&c) override {
        write(total_len_, c);
        write(dims_.data(), dims_.size(), c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        read(total_len_, c, remaining_length);
        read(dims_.data(), dims_.size(), c, remaining_length);
    }

    void postprocess_decode() override {}

    void postprocess_encode() override {}

    void preprocess_decode() override {}

   private:
    static size_t dims_product(const SZ3::SPERR::dims_type &dims) {
        return dims[0] * dims[1] * dims[2];
    }

    SZ3::SPERR::dims_type resolved_dims(size_t total_len) const {
        if (N == 1) {
            return SZ3::SPERR::dims_type{total_len, 1, 1};
        }

        SZ3::SPERR::dims_type dims = dims_;
        if (dims[0] == 0 || dims[1] == 0) {
            throw std::runtime_error("SPERR encoder requires non-zero dims for 2D/3D coding.");
        }

        if (N == 2) {
            dims[2] = 1;
        } else {
            if (dims[2] == 0) {
                throw std::runtime_error("SPERR encoder requires non-zero z dimension for 3D coding.");
            }
        }

        if (dims_product(dims) != total_len) {
            throw std::runtime_error("SPERR encoder dim/product mismatch with bin length.");
        }
        return dims;
    }

    static int64_t to_int64(T value) {
        if (std::numeric_limits<T>::is_signed) {
            return static_cast<int64_t>(value);
        }
        const uint64_t u = static_cast<uint64_t>(value);
        if (u > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
            throw std::runtime_error("SPERR value overflows int64.");
        }
        return static_cast<int64_t>(u);
    }

    static T from_int64(int64_t value) {
        if (std::numeric_limits<T>::is_signed) {
            if (value < static_cast<int64_t>(std::numeric_limits<T>::lowest()) ||
                value > static_cast<int64_t>(std::numeric_limits<T>::max())) {
                throw std::runtime_error("SPERR decoded value cannot fit output scalar type.");
            }
            return static_cast<T>(value);
        }

        if (value < 0) {
            throw std::runtime_error("SPERR decoded negative value cannot fit unsigned output scalar type.");
        }

        const uint64_t u = static_cast<uint64_t>(value);
        if (u > static_cast<uint64_t>(std::numeric_limits<T>::max())) {
            throw std::runtime_error("SPERR decoded value cannot fit unsigned output scalar type.");
        }
        return static_cast<T>(u);
    }

    static uint64_t signed_magnitude(int64_t value) {
        if (value == std::numeric_limits<int64_t>::min()) {
            throw std::runtime_error("SPERR signed value overflow when taking magnitude.");
        }
        return static_cast<uint64_t>(std::llabs(value));
    }

    static uint64_t max_magnitude(const std::vector<int64_t> &signed_vals) {
        uint64_t max_mag = 0;
        for (size_t i = 0; i < signed_vals.size(); i++) {
            max_mag = std::max(max_mag, signed_magnitude(signed_vals[i]));
        }
        return max_mag;
    }

    static std::vector<int64_t> to_signed_vector(const std::vector<T> &bins) {
        std::vector<int64_t> signed_vals(bins.size(), int64_t{0});
        for (size_t i = 0; i < bins.size(); i++) {
            signed_vals[i] = to_int64(bins[i]);
        }
        return signed_vals;
    }

    static std::vector<T> from_signed_vector(const std::vector<int64_t> &signed_vals) {
        std::vector<T> bins(signed_vals.size(), T{0});
        for (size_t i = 0; i < signed_vals.size(); i++) {
            bins[i] = from_int64(signed_vals[i]);
        }
        return bins;
    }

    static std::vector<uchar> encode_speck_1d(const std::vector<int64_t> &signed_vals) {
        const uint64_t max_mag = max_magnitude(signed_vals);
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint8_t>::max())) {
            return encode_speck_1d_impl<uint8_t>(signed_vals);
        }
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint16_t>::max())) {
            return encode_speck_1d_impl<uint16_t>(signed_vals);
        }
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
            return encode_speck_1d_impl<uint32_t>(signed_vals);
        }
        return encode_speck_1d_impl<uint64_t>(signed_vals);
    }

    static std::vector<int64_t> decode_speck_1d(size_t total_len, const uchar *cmpData, size_t cmpSize, size_t &consumed) {
        const uint8_t num_bitplanes = SZ3::SPERR::speck_int_get_num_bitplanes(cmpData);
        if (num_bitplanes <= 8) {
            return decode_speck_1d_impl<uint8_t>(total_len, cmpData, cmpSize, consumed);
        }
        if (num_bitplanes <= 16) {
            return decode_speck_1d_impl<uint16_t>(total_len, cmpData, cmpSize, consumed);
        }
        if (num_bitplanes <= 32) {
            return decode_speck_1d_impl<uint32_t>(total_len, cmpData, cmpSize, consumed);
        }
        return decode_speck_1d_impl<uint64_t>(total_len, cmpData, cmpSize, consumed);
    }

    static std::vector<uchar> encode_speck_2d(const std::vector<int64_t> &signed_vals, const SZ3::SPERR::dims_type &dims) {
        const uint64_t max_mag = max_magnitude(signed_vals);
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint8_t>::max())) {
            return encode_speck_2d_impl<uint8_t>(signed_vals, dims);
        }
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint16_t>::max())) {
            return encode_speck_2d_impl<uint16_t>(signed_vals, dims);
        }
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
            return encode_speck_2d_impl<uint32_t>(signed_vals, dims);
        }
        return encode_speck_2d_impl<uint64_t>(signed_vals, dims);
    }

    static std::vector<int64_t> decode_speck_2d(const SZ3::SPERR::dims_type &dims, const uchar *cmpData, size_t cmpSize,
                                                size_t &consumed) {
        const uint8_t num_bitplanes = SZ3::SPERR::speck_int_get_num_bitplanes(cmpData);
        if (num_bitplanes <= 8) {
            return decode_speck_2d_impl<uint8_t>(dims, cmpData, cmpSize, consumed);
        }
        if (num_bitplanes <= 16) {
            return decode_speck_2d_impl<uint16_t>(dims, cmpData, cmpSize, consumed);
        }
        if (num_bitplanes <= 32) {
            return decode_speck_2d_impl<uint32_t>(dims, cmpData, cmpSize, consumed);
        }
        return decode_speck_2d_impl<uint64_t>(dims, cmpData, cmpSize, consumed);
    }

    static std::vector<uchar> encode_speck_3d(const std::vector<int64_t> &signed_vals, const SZ3::SPERR::dims_type &dims) {
        const uint64_t max_mag = max_magnitude(signed_vals);
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint8_t>::max())) {
            return encode_speck_3d_impl<uint8_t>(signed_vals, dims);
        }
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint16_t>::max())) {
            return encode_speck_3d_impl<uint16_t>(signed_vals, dims);
        }
        if (max_mag <= static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
            return encode_speck_3d_impl<uint32_t>(signed_vals, dims);
        }
        return encode_speck_3d_impl<uint64_t>(signed_vals, dims);
    }

    static std::vector<int64_t> decode_speck_3d(const SZ3::SPERR::dims_type &dims, const uchar *cmpData, size_t cmpSize,
                                                size_t &consumed) {
        const uint8_t num_bitplanes = SZ3::SPERR::speck_int_get_num_bitplanes(cmpData);
        if (num_bitplanes <= 8) {
            return decode_speck_3d_impl<uint8_t>(dims, cmpData, cmpSize, consumed);
        }
        if (num_bitplanes <= 16) {
            return decode_speck_3d_impl<uint16_t>(dims, cmpData, cmpSize, consumed);
        }
        if (num_bitplanes <= 32) {
            return decode_speck_3d_impl<uint32_t>(dims, cmpData, cmpSize, consumed);
        }
        return decode_speck_3d_impl<uint64_t>(dims, cmpData, cmpSize, consumed);
    }

    template <class UI>
    static std::vector<uchar> encode_speck_1d_impl(const std::vector<int64_t> &signed_vals) {
        std::vector<UI> mags(signed_vals.size(), UI{0});
        SZ3::SPERR::Bitmask signs;
        signs.resize(signed_vals.size());
        signs.reset_true();

        for (size_t i = 0; i < signed_vals.size(); i++) {
            const uint64_t mag = signed_magnitude(signed_vals[i]);
            if (mag > static_cast<uint64_t>(std::numeric_limits<UI>::max())) {
                throw std::runtime_error("SPERR SPECK1D magnitude exceeds selected integer type.");
            }
            mags[i] = static_cast<UI>(mag);
            signs.wbit(i, signed_vals[i] >= 0);
        }

        SZ3::SPERR::SPECK1D_INT_ENC<UI> encoder;
        encoder.set_dims({signed_vals.size(), 1, 1});
        const SZ3::SPERR::RTNType rtn = encoder.use_coeffs(std::move(mags), std::move(signs));
        if (rtn != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR SPECK1D setup failed.");
        }
        encoder.encode();

        std::vector<uchar> stream;
        encoder.append_encoded_bitstream(stream);
        return stream;
    }

    template <class UI>
    static std::vector<int64_t> decode_speck_1d_impl(size_t total_len, const uchar *cmpData, size_t cmpSize,
                                                     size_t &consumed) {
        SZ3::SPERR::SPECK1D_INT_DEC<UI> decoder;
        decoder.set_dims({total_len, 1, 1});

        const size_t full_len = static_cast<size_t>(decoder.get_stream_full_len(cmpData));
        if (full_len > cmpSize) {
            throw std::runtime_error("SPERR SPECK1D stream length mismatch.");
        }

        decoder.use_bitstream(cmpData, full_len);
        decoder.decode();

        std::vector<UI> mags = decoder.release_coeffs();
        SZ3::SPERR::Bitmask signs = decoder.release_signs();
        if (mags.size() != total_len || signs.size() != total_len) {
            throw std::runtime_error("SPERR SPECK1D decoded length mismatch.");
        }

        std::vector<int64_t> signed_vals(total_len, int64_t{0});
        for (size_t i = 0; i < total_len; i++) {
            if (static_cast<uint64_t>(mags[i]) > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
                throw std::runtime_error("SPERR SPECK1D decoded magnitude overflows int64.");
            }
            const int64_t mag = static_cast<int64_t>(mags[i]);
            signed_vals[i] = signs.rbit(i) ? mag : -mag;
        }

        consumed = full_len;
        return signed_vals;
    }

    template <class UI>
    static std::vector<uchar> encode_speck_2d_impl(const std::vector<int64_t> &signed_vals,
                                                   const SZ3::SPERR::dims_type &dims) {
        const size_t total_len = dims[0] * dims[1];
        if (signed_vals.size() != total_len) {
            throw std::runtime_error("SPERR SPECK2D encode length mismatch.");
        }

        std::vector<UI> mags(signed_vals.size(), UI{0});
        SZ3::SPERR::Bitmask signs;
        signs.resize(signed_vals.size());
        signs.reset_true();

        for (size_t i = 0; i < signed_vals.size(); i++) {
            const uint64_t mag = signed_magnitude(signed_vals[i]);
            if (mag > static_cast<uint64_t>(std::numeric_limits<UI>::max())) {
                throw std::runtime_error("SPERR SPECK2D magnitude exceeds selected integer type.");
            }
            mags[i] = static_cast<UI>(mag);
            signs.wbit(i, signed_vals[i] >= 0);
        }

        SZ3::SPERR::SPECK2D_INT_ENC<UI> encoder;
        encoder.set_dims(dims);
        const SZ3::SPERR::RTNType rtn = encoder.use_coeffs(std::move(mags), std::move(signs));
        if (rtn != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR SPECK2D setup failed.");
        }
        encoder.encode();

        std::vector<uchar> stream;
        encoder.append_encoded_bitstream(stream);
        return stream;
    }

    template <class UI>
    static std::vector<int64_t> decode_speck_2d_impl(const SZ3::SPERR::dims_type &dims, const uchar *cmpData,
                                                     size_t cmpSize, size_t &consumed) {
        const size_t total_len = dims[0] * dims[1];

        SZ3::SPERR::SPECK2D_INT_DEC<UI> decoder;
        decoder.set_dims(dims);

        const size_t full_len = static_cast<size_t>(decoder.get_stream_full_len(cmpData));
        if (full_len > cmpSize) {
            throw std::runtime_error("SPERR SPECK2D stream length mismatch.");
        }

        decoder.use_bitstream(cmpData, full_len);
        decoder.decode();

        std::vector<UI> mags = decoder.release_coeffs();
        SZ3::SPERR::Bitmask signs = decoder.release_signs();
        if (mags.size() != total_len || signs.size() != total_len) {
            throw std::runtime_error("SPERR SPECK2D decoded length mismatch.");
        }

        std::vector<int64_t> signed_vals(total_len, int64_t{0});
        for (size_t i = 0; i < total_len; i++) {
            if (static_cast<uint64_t>(mags[i]) > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
                throw std::runtime_error("SPERR SPECK2D decoded magnitude overflows int64.");
            }
            const int64_t mag = static_cast<int64_t>(mags[i]);
            signed_vals[i] = signs.rbit(i) ? mag : -mag;
        }

        consumed = full_len;
        return signed_vals;
    }

    template <class UI>
    static std::vector<uchar> encode_speck_3d_impl(const std::vector<int64_t> &signed_vals,
                                                   const SZ3::SPERR::dims_type &dims) {
        const size_t total_len = dims[0] * dims[1] * dims[2];
        if (signed_vals.size() != total_len) {
            throw std::runtime_error("SPERR SPECK3D encode length mismatch.");
        }

        std::vector<UI> mags(signed_vals.size(), UI{0});
        SZ3::SPERR::Bitmask signs;
        signs.resize(signed_vals.size());
        signs.reset_true();

        for (size_t i = 0; i < signed_vals.size(); i++) {
            const uint64_t mag = signed_magnitude(signed_vals[i]);
            if (mag > static_cast<uint64_t>(std::numeric_limits<UI>::max())) {
                throw std::runtime_error("SPERR SPECK3D magnitude exceeds selected integer type.");
            }
            mags[i] = static_cast<UI>(mag);
            signs.wbit(i, signed_vals[i] >= 0);
        }

        SZ3::SPERR::SPECK3D_INT_ENC<UI> encoder;
        encoder.set_dims(dims);
        const SZ3::SPERR::RTNType rtn = encoder.use_coeffs(std::move(mags), std::move(signs));
        if (rtn != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR SPECK3D setup failed.");
        }
        encoder.encode();

        std::vector<uchar> stream;
        encoder.append_encoded_bitstream(stream);
        return stream;
    }

    template <class UI>
    static std::vector<int64_t> decode_speck_3d_impl(const SZ3::SPERR::dims_type &dims, const uchar *cmpData,
                                                     size_t cmpSize, size_t &consumed) {
        const size_t total_len = dims[0] * dims[1] * dims[2];

        SZ3::SPERR::SPECK3D_INT_DEC<UI> decoder;
        decoder.set_dims(dims);

        const size_t full_len = static_cast<size_t>(decoder.get_stream_full_len(cmpData));
        if (full_len > cmpSize) {
            throw std::runtime_error("SPERR SPECK3D stream length mismatch.");
        }

        decoder.use_bitstream(cmpData, full_len);
        decoder.decode();

        std::vector<UI> mags = decoder.release_coeffs();
        SZ3::SPERR::Bitmask signs = decoder.release_signs();
        if (mags.size() != total_len || signs.size() != total_len) {
            throw std::runtime_error("SPERR SPECK3D decoded length mismatch.");
        }

        std::vector<int64_t> signed_vals(total_len, int64_t{0});
        for (size_t i = 0; i < total_len; i++) {
            if (static_cast<uint64_t>(mags[i]) > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
                throw std::runtime_error("SPERR SPECK3D decoded magnitude overflows int64.");
            }
            const int64_t mag = static_cast<int64_t>(mags[i]);
            signed_vals[i] = signs.rbit(i) ? mag : -mag;
        }

        consumed = full_len;
        return signed_vals;
    }

    SZ3::SPERR::dims_type dims_;
    size_t total_len_ = 0;
};

}  // namespace SZ3

#endif
