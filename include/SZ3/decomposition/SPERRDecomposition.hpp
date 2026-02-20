#ifndef SZ3_SPERR_DECOMPOSITION_HPP
#define SZ3_SPERR_DECOMPOSITION_HPP

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/thirdparty/sperr/SPERRHeaderOnly.hpp"

namespace SZ3 {

/**
 * @file SPERRDecomposition.hpp
 * @brief SPERR core workflow mapped onto SZ3's decomposition module interface.
 *
 * Compression pipeline implemented in this module:
 * 1. `prepare`: validate config and map SZ3 dims `[z,y,x]` to SPERR dims `[x,y,z]`.
 * 2. `forward`: conditioning + 3D CDF97 forward wavelet transform.
 * 3. `quantize_and_collect_outliers`: compute `q`, midtread quantize, and in PWE mode
 *    detect outliers from inverse-wavelet residuals.
 * 4. `encode_frame`: pack conditioner header + SPECK integers + optional outliers.
 *
 * Decompression pipeline implemented in this module:
 * 1. `decode_frame`: parse conditioner header + SPECK stream + optional outliers.
 * 2. `inverse`: inverse quantization + inverse CDF97 + outlier restore + inverse conditioning.
 */

/// Quality control mode interpreted from SZ3 `errorBoundMode`.
enum class SPERRMode { PSNR_MODE, PWE_MODE };

/// 3D shape in SPERR ordering `[x, y, z]`.
struct SPERR3DLayout {
    size_t dimx = 0;
    size_t dimy = 0;
    size_t dimz = 0;
};

using SPERRUIntVec = std::variant<std::vector<uint8_t>, std::vector<uint16_t>, std::vector<uint32_t>,
                                  std::vector<uint64_t>>;

/// Intermediate state exchanged between decomposition and encoding logic.
struct SPERRFrame {
    SZ3::SPERR::dims_type dims = {0, 0, 0};
    SZ3::SPERR::condi_type conditioner_header{};
    bool constant_field = false;

    SPERRMode mode = SPERRMode::PWE_MODE;
    double quality = 0.0;
    double conditioned_range = 0.0;
    double q = 0.0;

    SZ3::SPERR::UINTType uint_flag = SZ3::SPERR::UINTType::UINT8;
    SPERRUIntVec coeffs_ui;
    SZ3::SPERR::Bitmask sign_array;
    std::vector<double> wavelet_coeffs;
    std::vector<double> conditioned_values;
    std::vector<SZ3::SPERR::Outlier> outliers;

    SPERRFrame() : coeffs_ui(std::vector<uint8_t>{}) {}
};

/**
 * @brief SPERR decomposition implementation through `DecompositionInterface<T, uchar, N>`.
 *
 * Note: for SPERR integration, this decomposition module also owns frame packaging logic so
 * `compress/decompress` can be used directly through the generic SZ3 module interface.
 */
template <class T, uint N>
class SPERRDecomposition : public concepts::DecompositionInterface<T, uchar, N> {
   public:
    /**
     * @brief Full SPERR compression path through decomposition interface.
     *
     * Steps:
     * 1. shape preparation
     * 2. forward conditioning + wavelet
     * 3. quantization and optional outlier collection
     * 4. frame encoding to byte bins
     */
    std::vector<uchar> compress(const Config &conf, T *data) override {
        auto layout = prepare(conf);
        auto frame = forward(conf, layout, data);
        if (!frame.constant_field) {
            quantize_and_collect_outliers(frame);
        }
        return encode_frame(frame);
    }

    /**
     * @brief Full SPERR decompression path through decomposition interface.
     *
     * Steps:
     * 1. frame decoding from bins
     * 2. inverse quantization + inverse wavelet + inverse conditioning
     */
    T *decompress(const Config &conf, std::vector<uchar> &quant_inds, T *dec_data) override {
        auto layout = prepare(conf);
        auto frame = decode_frame(layout, quant_inds.data(), quant_inds.size());
        return inverse(conf, frame, dec_data);
    }

    void save(uchar *& /*c*/) override {}

    void load(const uchar *& /*c*/, size_t & /*remaining_length*/) override {}

    std::pair<uchar, uchar> get_out_range() override { return {0, std::numeric_limits<uchar>::max()}; }

    /// Validate SPERR constraints and convert SZ3 dimensions to SPERR layout.
    SPERR3DLayout prepare(const Config &conf) const {
        if constexpr (!std::is_floating_point<T>::value) {
            throw std::invalid_argument("SPERR supports floating-point data only.");
        }
        if (N != 3 || conf.N != 3 || conf.dims.size() != 3) {
            throw std::invalid_argument("SPERR core integration supports 3D data only.");
        }

        SPERR3DLayout layout;
        // SZ3 stores dims in slowest->fastest order for 3D: [z, y, x].
        // SPERR expects [x, y, z].
        layout.dimx = conf.dims[2];
        layout.dimy = conf.dims[1];
        layout.dimz = conf.dims[0];
        return layout;
    }

    /// Forward stage: conditioning then 3D CDF97 DWT.
    SPERRFrame forward(const Config &conf, const SPERR3DLayout &layout, const T *data) const {
        SPERRFrame frame;
        const auto qc = resolve_quality_control(conf);
        frame.dims = {layout.dimx, layout.dimy, layout.dimz};
        frame.mode = qc.mode;
        frame.quality = qc.quality;

        std::vector<double> values(conf.num);
        std::copy(data, data + conf.num, values.begin());

        auto conditioner = SZ3::SPERR::Conditioner();
        frame.conditioner_header = conditioner.condition(values, frame.dims);
        frame.constant_field = conditioner.is_constant(frame.conditioner_header[0]);
        if (frame.constant_field) {
            return frame;
        }

        if (frame.mode == SPERRMode::PWE_MODE) {
            frame.conditioned_values = values;
        }

        const auto minmax = std::minmax_element(values.begin(), values.end());
        frame.conditioned_range = *minmax.second - *minmax.first;

        auto cdf = SZ3::SPERR::CDF97();
        auto rtn = cdf.take_data(std::move(values), frame.dims);
        if (rtn != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR forward wavelet setup failed.");
        }
        cdf.dwt3d();
        frame.wavelet_coeffs = cdf.release_data();

        return frame;
    }

    /// Quantization stage: q estimation, midtread quantization, optional PWE outlier detection.
    void quantize_and_collect_outliers(SPERRFrame &frame) const {
        if (frame.constant_field) {
            return;
        }
        const auto q = estimate_q(frame);
        if (!(q > 0.0)) {
            throw std::invalid_argument("SPERR quantization step must be positive.");
        }
        frame.q = q;

        auto conditioner = SZ3::SPERR::Conditioner();
        conditioner.save_q(frame.conditioner_header, frame.q);

        midtread_quantize(frame.wavelet_coeffs, frame.q, frame.coeffs_ui, frame.sign_array, frame.uint_flag);

        if (frame.mode == SPERRMode::PWE_MODE) {
            auto reconstructed = midtread_inverse_quantize(frame.coeffs_ui, frame.sign_array, frame.q);
            auto cdf = SZ3::SPERR::CDF97();
            auto rtn = cdf.take_data(std::move(reconstructed), frame.dims);
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR inverse wavelet setup failed during outlier detection.");
            }
            cdf.idwt3d();
            reconstructed = cdf.release_data();

            frame.outliers.clear();
            frame.outliers.reserve(frame.conditioned_values.size() / 20);
            for (size_t i = 0; i < frame.conditioned_values.size(); i++) {
                const auto diff = frame.conditioned_values[i] - reconstructed[i];
                if (std::abs(diff) > frame.quality) {
                    frame.outliers.emplace_back(i, diff);
                }
            }
        }

        frame.wavelet_coeffs.clear();
        frame.wavelet_coeffs.shrink_to_fit();
        frame.conditioned_values.clear();
        frame.conditioned_values.shrink_to_fit();
    }

    /// Inverse stage: inverse quantization, inverse wavelet, outlier restore, inverse conditioning.
    T *inverse(const Config &conf, SPERRFrame &frame, T *dec_data) const {
        std::vector<double> values;
        auto conditioner = SZ3::SPERR::Conditioner();

        if (frame.constant_field) {
            auto rtn = conditioner.inverse_condition(values, frame.dims, frame.conditioner_header);
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR inverse conditioning failed for constant field.");
            }
        } else {
            values = midtread_inverse_quantize(frame.coeffs_ui, frame.sign_array, frame.q);

            auto cdf = SZ3::SPERR::CDF97();
            auto rtn = cdf.take_data(std::move(values), frame.dims);
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR inverse wavelet setup failed.");
            }
            cdf.idwt3d();
            values = cdf.release_data();

            for (const auto &outlier : frame.outliers) {
                if (outlier.pos >= values.size()) {
                    throw std::runtime_error("SPERR outlier index out of bounds.");
                }
                values[outlier.pos] += outlier.err;
            }

            rtn = conditioner.inverse_condition(values, frame.dims, frame.conditioner_header);
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR inverse conditioning failed.");
            }
        }

        if (values.size() != conf.num) {
            throw std::runtime_error("SPERR reconstruction length mismatch.");
        }
        for (size_t i = 0; i < conf.num; i++) {
            dec_data[i] = static_cast<T>(values[i]);
        }
        return dec_data;
    }

   private:
    struct SPERRQualityControl {
        SPERRMode mode = SPERRMode::PWE_MODE;
        double quality = 0.0;
    };

    static SPERRQualityControl resolve_quality_control(const Config &conf) {
        SPERRQualityControl qc;
        if (conf.errorBoundMode == EB_PSNR) {
            qc.mode = SPERRMode::PSNR_MODE;
            qc.quality = conf.psnrErrorBound;
        } else {
            qc.mode = SPERRMode::PWE_MODE;
            qc.quality = conf.absErrorBound;
        }
        if (qc.quality <= 0.0) {
            throw std::invalid_argument("SPERR requires a positive quality target.");
        }
        return qc;
    }

    static double estimate_mse_midtread(const std::vector<double> &vals, double q) {
        const auto len = vals.size();
        constexpr size_t stride_size = 4096;
        const auto num_strides = len / stride_size;
        auto accum = std::vector<double>(num_strides + 1, 0.0);

        for (size_t i = 0; i < num_strides; i++) {
            const auto begin = vals.cbegin() + i * stride_size;
            accum[i] = std::accumulate(begin, begin + stride_size, 0.0, [q](double init, double v) {
                const auto diff = std::remainder(v, q);
                return init + diff * diff;
            });
        }

        accum[num_strides] =
            std::accumulate(vals.cbegin() + num_strides * stride_size, vals.cend(), 0.0, [q](double init, double v) {
                const auto diff = std::remainder(v, q);
                return init + diff * diff;
            });

        const auto total = std::accumulate(accum.cbegin(), accum.cend(), 0.0);
        return total / static_cast<double>(len);
    }

    static double estimate_q(const SPERRFrame &frame) {
        if (frame.mode == SPERRMode::PWE_MODE) {
            return frame.quality * 1.5;
        }
        if (frame.conditioned_range <= 0.0) {
            throw std::invalid_argument("SPERR PSNR mode requires positive conditioned data range.");
        }
        if (frame.wavelet_coeffs.empty()) {
            throw std::invalid_argument("SPERR q-estimation requires wavelet coefficients.");
        }
        const auto target_mse =
            (frame.conditioned_range * frame.conditioned_range) * std::pow(10.0, -frame.quality / 10.0);
        auto q = 2.0 * std::sqrt(target_mse * 3.0);
        while (estimate_mse_midtread(frame.wavelet_coeffs, q) > target_mse) {
            q /= std::exp2(0.25);
        }
        return q;
    }

    template <class UI>
    static void encode_speck(SPERRFrame &frame, std::vector<uchar> &stream) {
        auto encoder = SZ3::SPERR::SPECK3D_INT_ENC<UI>();
        encoder.set_dims(frame.dims);
        const auto rtn =
            encoder.use_coeffs(std::move(std::get<std::vector<UI>>(frame.coeffs_ui)), std::move(frame.sign_array));
        if (rtn != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR integer coefficient setup failed.");
        }
        encoder.encode();
        encoder.append_encoded_bitstream(stream);
    }

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

    static std::vector<uchar> encode_frame(SPERRFrame &frame) {
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

    static SPERRFrame decode_frame(const SPERR3DLayout &layout, const uchar *cmpData, size_t cmpSize) {
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

    static void midtread_quantize(const std::vector<double> &vals_d, double q, SPERRUIntVec &vals_ui,
                                  SZ3::SPERR::Bitmask &signs, SZ3::SPERR::UINTType &uint_flag) {
        std::fesetround(FE_TONEAREST);

        const auto max_itr = std::max_element(vals_d.cbegin(), vals_d.cend(),
                                              [](auto a, auto b) { return std::abs(a) < std::abs(b); });
        const auto maxll = std::llrint(std::abs(*max_itr) / q);

        if (maxll <= std::numeric_limits<uint8_t>::max()) {
            uint_flag = SZ3::SPERR::UINTType::UINT8;
            vals_ui = std::vector<uint8_t>(vals_d.size());
        } else if (maxll <= std::numeric_limits<uint16_t>::max()) {
            uint_flag = SZ3::SPERR::UINTType::UINT16;
            vals_ui = std::vector<uint16_t>(vals_d.size());
        } else if (maxll <= std::numeric_limits<uint32_t>::max()) {
            uint_flag = SZ3::SPERR::UINTType::UINT32;
            vals_ui = std::vector<uint32_t>(vals_d.size());
        } else {
            uint_flag = SZ3::SPERR::UINTType::UINT64;
            vals_ui = std::vector<uint64_t>(vals_d.size());
        }

        signs.resize(vals_d.size());

        std::visit(
            [&vals_d, &signs, q](auto &vec) {
                using VecT = std::decay_t<decltype(vec)>;
                using ElemT = typename VecT::value_type;
                const auto inv = 1.0 / q;
                const auto bits_x64 = vals_d.size() - vals_d.size() % 64;

                for (size_t i = 0; i < bits_x64; i += 64) {
                    auto bits64 = uint64_t{0};
                    for (size_t j = 0; j < 64; j++) {
                        const auto ll = std::llrint(vals_d[i + j] * inv);
                        bits64 |= (static_cast<uint64_t>(ll >= 0) << j);
                        vec[i + j] = static_cast<ElemT>(std::llabs(ll));
                    }
                    signs.wlong(i, bits64);
                }

                for (size_t i = bits_x64; i < vals_d.size(); i++) {
                    const auto ll = std::llrint(vals_d[i] * inv);
                    signs.wbit(i, ll >= 0);
                    vec[i] = static_cast<ElemT>(std::llabs(ll));
                }
            },
            vals_ui);
    }

    static std::vector<double> midtread_inverse_quantize(const SPERRUIntVec &vals_ui, const SZ3::SPERR::Bitmask &signs,
                                                          double q) {
        std::vector<double> vals_d(signs.size());
        const auto sign_scale = std::array<double, 2>{-1.0, 1.0};

        std::visit(
            [&vals_d, &signs, q, &sign_scale](const auto &vec) {
                if (vec.size() != signs.size()) {
                    throw std::runtime_error("SPERR quantized/sign array size mismatch.");
                }
                const auto bits_x64 = vals_d.size() - vals_d.size() % 64;

                for (size_t i = 0; i < bits_x64; i += 64) {
                    const auto bits64 = signs.rlong(i);
                    for (size_t j = 0; j < 64; j++) {
                        const auto bit = (bits64 >> j) & uint64_t{1};
                        vals_d[i + j] = q * static_cast<double>(vec[i + j]) * sign_scale[bit];
                    }
                }

                for (size_t i = bits_x64; i < vals_d.size(); i++) {
                    vals_d[i] = q * static_cast<double>(vec[i]) * sign_scale[signs.rbit(i)];
                }
            },
            vals_ui);

        return vals_d;
    }
};

}  // namespace SZ3

#endif
