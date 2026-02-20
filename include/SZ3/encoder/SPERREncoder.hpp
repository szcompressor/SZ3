#ifndef SZ3_SPERR_ENCODER_HPP
#define SZ3_SPERR_ENCODER_HPP

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "SZ3/decomposition/SPERRDecomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"

namespace SZ3 {

/**
 * @file SPERREncoder.hpp
 * @brief SPERR bitstream assembly/parsing mapped to SZ3's encoder interface.
 *
 * Encoded frame layout produced by this module:
 * 1. Conditioner header (`SZ3::SPERR::condi_type`).
 * 2. SPECK3D integer payload (quantized wavelet coefficients + signs).
 * 3. Optional Outlier_Coder payload (PWE mode only when outliers exist).
 */

/**
 * @brief SPERR frame encoder/decoder implementing `EncoderInterface<uchar>`.
 *
 * `EncoderInterface<uchar>` methods are pass-through byte copy methods used by
 * `SPERRCompressor`'s generic module pipeline.
 * `encode(SPERRFrame&)` / `decode(layout,...)` perform actual SPERR stream coding.
 */
template <class T>
class SPERREncoder : public concepts::EncoderInterface<uchar> {
   public:
    void preprocess_encode(const std::vector<uchar> & /*bins*/, int /*stateNum*/) override {}

    size_t encode(const std::vector<uchar> &bins, uchar *&bytes) override {
        std::memcpy(bytes, bins.data(), bins.size());
        bytes += bins.size();
        return bins.size();
    }

    std::vector<uchar> decode(const uchar *&bytes, size_t targetLength) override {
        std::vector<uchar> out(targetLength);
        std::memcpy(out.data(), bytes, targetLength);
        bytes += targetLength;
        return out;
    }

    void save(uchar *& /*c*/) override {}

    void load(const uchar *& /*c*/, size_t & /*remaining_length*/) override {}

    void postprocess_decode() override {}

    void postprocess_encode() override {}

    void preprocess_decode() override {}

    /// Encode SPERR frame into `[conditioner][SPECK][optional outliers]`.
    std::vector<uchar> encode(SPERRFrame &frame) const {
        auto stream = std::vector<uchar>();
        stream.insert(stream.end(), frame.conditioner_header.begin(), frame.conditioner_header.end());

        if (frame.constant_field) {
            return stream;
        }

        switch (frame.uint_flag) {
            case SZ3::SPERR::UINTType::UINT8:
                encode_speck<uint8_t>(frame, stream);
                break;
            case SZ3::SPERR::UINTType::UINT16:
                encode_speck<uint16_t>(frame, stream);
                break;
            case SZ3::SPERR::UINTType::UINT32:
                encode_speck<uint32_t>(frame, stream);
                break;
            case SZ3::SPERR::UINTType::UINT64:
                encode_speck<uint64_t>(frame, stream);
                break;
            default:
                throw std::runtime_error("SPERR invalid integer type for encoding.");
        }

        if (!frame.outliers.empty()) {
            auto outlier_coder = SZ3::SPERR::Outlier_Coder();
            outlier_coder.set_length(frame.dims[0] * frame.dims[1] * frame.dims[2]);
            outlier_coder.set_tolerance(frame.quality);
            outlier_coder.use_outlier_list(std::move(frame.outliers));
            const auto rtn = outlier_coder.encode();
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR outlier encoding failed.");
            }
            outlier_coder.append_encoded_bitstream(stream);
        }

        return stream;
    }

    /// Decode SPERR frame from `[conditioner][SPECK][optional outliers]`.
    SPERRFrame decode(const SPERR3DLayout &layout, const uchar *cmpData, size_t cmpSize) const {
        auto frame = SPERRFrame();
        frame.dims = {layout.dimx, layout.dimy, layout.dimz};

        if (cmpSize < frame.conditioner_header.size()) {
            throw std::runtime_error("SPERR stream too short.");
        }

        std::copy_n(cmpData, frame.conditioner_header.size(), frame.conditioner_header.begin());
        const auto conditioner = SZ3::SPERR::Conditioner();
        frame.constant_field = conditioner.is_constant(frame.conditioner_header[0]);

        if (frame.constant_field) {
            return frame;
        }

        frame.q = conditioner.retrieve_q(frame.conditioner_header);
        if (!(frame.q > 0.0)) {
            throw std::runtime_error("SPERR invalid quantization step in stream.");
        }

        auto pos = cmpData + frame.conditioner_header.size();
        size_t remaining = cmpSize - frame.conditioner_header.size();
        if (remaining < SZ3::SPERR::SPECK_INT<uint8_t>::header_size) {
            throw std::runtime_error("SPERR stream missing SPECK payload.");
        }

        const auto num_bitplanes = SZ3::SPERR::speck_int_get_num_bitplanes(pos);
        if (num_bitplanes <= 8) {
            frame.uint_flag = SZ3::SPERR::UINTType::UINT8;
            remaining -= decode_speck<uint8_t>(frame, pos, remaining);
        } else if (num_bitplanes <= 16) {
            frame.uint_flag = SZ3::SPERR::UINTType::UINT16;
            remaining -= decode_speck<uint16_t>(frame, pos, remaining);
        } else if (num_bitplanes <= 32) {
            frame.uint_flag = SZ3::SPERR::UINTType::UINT32;
            remaining -= decode_speck<uint32_t>(frame, pos, remaining);
        } else {
            frame.uint_flag = SZ3::SPERR::UINTType::UINT64;
            remaining -= decode_speck<uint64_t>(frame, pos, remaining);
        }
        pos = cmpData + cmpSize - remaining;

        if (remaining > 0) {
            if (remaining < SZ3::SPERR::SPECK_INT<uint8_t>::header_size) {
                throw std::runtime_error("SPERR outlier stream is truncated.");
            }
            auto outlier_coder = SZ3::SPERR::Outlier_Coder();
            outlier_coder.set_length(frame.dims[0] * frame.dims[1] * frame.dims[2]);
            outlier_coder.set_tolerance(frame.q / 1.5);
            const auto outlier_len = outlier_coder.get_stream_full_len(pos);
            if (outlier_len != remaining) {
                throw std::runtime_error("SPERR outlier stream length mismatch.");
            }
            auto rtn = outlier_coder.use_bitstream(pos, remaining);
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR outlier stream parse failed.");
            }
            rtn = outlier_coder.decode();
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR outlier decoding failed.");
            }
            frame.outliers = outlier_coder.view_outlier_list();
        }

        return frame;
    }

   private:
    /// SPECK integer encode helper for a specific integer coefficient width.
    template <class UI>
    static void encode_speck(SPERRFrame &frame, std::vector<uchar> &stream) {
        auto encoder = SZ3::SPERR::SPECK3D_INT_ENC<UI>();
        encoder.set_dims(frame.dims);
        auto rtn =
            encoder.use_coeffs(std::move(std::get<std::vector<UI>>(frame.coeffs_ui)), std::move(frame.sign_array));
        if (rtn != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR integer coefficient setup failed.");
        }
        encoder.encode();
        encoder.append_encoded_bitstream(stream);
    }

    /// SPECK integer decode helper for a specific integer coefficient width.
    template <class UI>
    static size_t decode_speck(SPERRFrame &frame, const uchar *ptr, size_t remaining) {
        auto decoder = SZ3::SPERR::SPECK3D_INT_DEC<UI>();
        decoder.set_dims(frame.dims);
        const auto speck_len = decoder.get_stream_full_len(ptr);
        if (speck_len > remaining) {
            throw std::runtime_error("SPERR SPECK stream length mismatch.");
        }
        decoder.use_bitstream(ptr, static_cast<size_t>(speck_len));
        decoder.decode();
        frame.coeffs_ui = decoder.release_coeffs();
        frame.sign_array = decoder.release_signs();
        return static_cast<size_t>(speck_len);
    }
};

}  // namespace SZ3

#endif
