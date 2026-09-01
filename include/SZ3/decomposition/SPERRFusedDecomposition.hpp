#ifndef SZ3_SPERR_FUSED_DECOMPOSITION_HPP
#define SZ3_SPERR_FUSED_DECOMPOSITION_HPP

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/encoder/SPERREncoder.hpp"
#include "SZ3/quantizer/ScalarQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/thirdparty/sperr/SPERRTypes.hpp"

namespace SZ3 {

/**
 * @file SPERRFusedDecomposition.hpp
 * @brief SPERR decomposition module.
 */
template <class T, uint N>
class SPERRFusedDecomposition : public concepts::DecompositionInterface<T, int64_t, N> {
   public:
    /**
     * @brief Lossy SPERR decomposition stage for 3D floating-point fields.
     *
     * Designed for compressor pipelines that want SPERR transform+quantization
     * as decomposition while keeping encoder/lossless stages composable.
     */
    std::vector<int64_t> compress(const Config &conf, T *data) override {
        const auto layout = prepare(conf);
        SPERRFrame frame = forward(conf, layout, data);
        if (!frame.constant_field) {
            quantize_and_collect_outliers(frame);
        }

        coeff_stream_ = encode_coefficient_stream(frame);
        return build_outlier_map(conf, frame);
    }

    T *decompress(const Config &conf, std::vector<int64_t> &outlier_bins, T *dec_data) override {
        if (coeff_stream_.empty()) {
            throw std::runtime_error("SPERR missing coefficient stream in decomposition state.");
        }

        const auto layout = prepare(conf);
        SPERRFrame frame = decode_coefficient_stream(layout, coeff_stream_.data(), coeff_stream_.size());
        return inverse(conf, frame, outlier_bins, dec_data);
    }

    void save(uchar *&c) override {
        write(coeff_stream_.size(), c);
        if (!coeff_stream_.empty()) {
            write(coeff_stream_.data(), coeff_stream_.size(), c);
        }
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        size_t coeff_len = 0;
        read(coeff_len, c, remaining_length);
        coeff_stream_.resize(coeff_len);
        if (coeff_len > 0) {
            read(coeff_stream_.data(), coeff_len, c, remaining_length);
        }
    }

    size_t size_est() override { return sizeof(size_t) + coeff_stream_.size(); }

    std::pair<int64_t, int64_t> get_out_range() override {
        // SPERR bins are a SPECK bitstream, not a quantizer domain; 0 means "no bin range".
        return std::make_pair<int64_t, int64_t>(0, 0);
    }

    SPERR3DLayout prepare(const Config &conf) const {
        if (!std::is_floating_point<T>::value) {
            throw std::invalid_argument("SPERR supports floating-point data only.");
        }
        if (N != 3 || conf.N != 3 || conf.dims.size() != 3) {
            throw std::invalid_argument("SPERR decomposition supports 3D data only.");
        }

        SPERR3DLayout layout;
        // SZ3 dims: [z, y, x], SPERR dims: [x, y, z].
        layout.dimx = conf.dims[2];
        layout.dimy = conf.dims[1];
        layout.dimz = conf.dims[0];
        return layout;
    }

    SPERRFrame forward(const Config &conf, const SPERR3DLayout &layout, const T *data) const {
        SPERRFrame frame;
        const SPERRQualityControl qc = resolve_quality_control(conf);
        frame.dims = {layout.dimx, layout.dimy, layout.dimz};
        frame.mode = qc.mode;
        frame.quality = qc.quality;

        std::vector<double> values(conf.num);
        std::copy(data, data + conf.num, values.begin());

        SZ3::SPERR::Conditioner conditioner;
        frame.conditioner_header = conditioner.condition(values, frame.dims);
        frame.constant_field = conditioner.is_constant(frame.conditioner_header[0]);
        if (frame.constant_field) {
            return frame;
        }

        if (frame.mode == SPERRMode::PWE_MODE) {
            frame.conditioned_values = values;
        }

        const std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> minmax =
            std::minmax_element(values.begin(), values.end());
        frame.conditioned_range = *minmax.second - *minmax.first;

        SZ3::SPERR::CDF97 cdf;
        const SZ3::SPERR::RTNType rtn = cdf.take_data(std::move(values), frame.dims);
        if (rtn != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR forward wavelet setup failed.");
        }
        cdf.dwt3d();
        frame.wavelet_coeffs = cdf.release_data();

        return frame;
    }

    void quantize_and_collect_outliers(SPERRFrame &frame) const {
        if (frame.constant_field) {
            return;
        }

        const double q = estimate_q(frame);
        if (!(q > 0.0)) {
            throw std::invalid_argument("SPERR quantization step must be positive.");
        }
        frame.q = q;

        SZ3::SPERR::Conditioner conditioner;
        conditioner.save_q(frame.conditioner_header, frame.q);

        std::fesetround(FE_TONEAREST);

        ScalarQuantizer<double, int64_t> coeff_quantizer = make_coeff_quantizer(frame.q);
        frame.coeff_bins.resize(frame.wavelet_coeffs.size(), int64_t{0});

        std::vector<double> reconstructed_coeffs;
        if (frame.mode == SPERRMode::PWE_MODE) {
            reconstructed_coeffs.resize(frame.wavelet_coeffs.size(), 0.0);
        }

        for (size_t i = 0; i < frame.wavelet_coeffs.size(); i++) {
            double coeff = frame.wavelet_coeffs[i];
            frame.coeff_bins[i] = coeff_quantizer.quantize_and_overwrite(coeff, 0.0);
            if (frame.mode == SPERRMode::PWE_MODE) {
                reconstructed_coeffs[i] = coeff;
            }
        }

        if (frame.mode == SPERRMode::PWE_MODE) {
            SZ3::SPERR::CDF97 cdf;
            const SZ3::SPERR::RTNType rtn = cdf.take_data(std::move(reconstructed_coeffs), frame.dims);
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR inverse wavelet setup failed during outlier detection.");
            }
            cdf.idwt3d();
            std::vector<double> reconstructed = cdf.release_data();

            frame.outliers.clear();
            frame.outliers.reserve(frame.conditioned_values.size() / 20);
            for (size_t i = 0; i < frame.conditioned_values.size(); i++) {
                const double diff = frame.conditioned_values[i] - reconstructed[i];
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

    T *inverse(const Config &conf, SPERRFrame &frame, const std::vector<int64_t> &outlier_bins, T *dec_data) const {
        std::vector<double> values;
        SZ3::SPERR::Conditioner conditioner;

        if (frame.constant_field) {
            validate_outlier_map(conf, outlier_bins, conf.num);
            if (std::any_of(outlier_bins.begin(), outlier_bins.end(), [](int64_t code) { return code != 0; })) {
                throw std::runtime_error("SPERR constant field should not contain outlier payload.");
            }

            const SZ3::SPERR::RTNType rtn = conditioner.inverse_condition(values, frame.dims, frame.conditioner_header);
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR inverse conditioning failed for constant field.");
            }
        } else {
            values = recover_coefficients(frame.coeff_bins, frame.q);

            SZ3::SPERR::CDF97 cdf;
            SZ3::SPERR::RTNType rtn = cdf.take_data(std::move(values), frame.dims);
            if (rtn != SZ3::SPERR::RTNType::Good) {
                throw std::runtime_error("SPERR inverse wavelet setup failed.");
            }
            cdf.idwt3d();
            values = cdf.release_data();

            apply_outlier_map(conf, outlier_bins, values);

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

    static ScalarQuantizer<double, int64_t> make_coeff_quantizer(double step) {
        return ScalarQuantizer<double, int64_t>(step);
    }

    static ScalarQuantizer<double, int64_t> make_outlier_quantizer(double tolerance) {
        return ScalarQuantizer<double, int64_t>(tolerance, 1.1, -0.25);
    }

    static std::vector<double> recover_coefficients(const std::vector<int64_t> &coeff_bins, double q) {
        if (!(q > 0.0)) {
            throw std::runtime_error("SPERR invalid quantization step for coefficient recovery.");
        }

        ScalarQuantizer<double, int64_t> quantizer = make_coeff_quantizer(q);
        std::vector<double> coeffs(coeff_bins.size(), 0.0);
        for (size_t i = 0; i < coeff_bins.size(); i++) {
            coeffs[i] = quantizer.recover(0.0, coeff_bins[i]);
        }
        return coeffs;
    }

    static void validate_outlier_map(const Config &conf, const std::vector<int64_t> &bins, size_t expected_len) {
        if (bins.size() != expected_len) {
            throw std::runtime_error("SPERR outlier map length mismatch.");
        }
        const SPERRQualityControl qc = resolve_quality_control(conf);
        if (qc.mode != SPERRMode::PWE_MODE) {
            for (size_t i = 0; i < bins.size(); i++) {
                if (bins[i] != 0) {
                    throw std::runtime_error("SPERR PSNR mode should not carry non-zero outlier corrections.");
                }
            }
        }
    }

    static std::vector<int64_t> build_outlier_map(const Config &conf, const SPERRFrame &frame) {
        std::vector<int64_t> outlier_map(conf.num, int64_t{0});
        if (frame.constant_field || frame.mode != SPERRMode::PWE_MODE) {
            return outlier_map;
        }

        ScalarQuantizer<double, int64_t> quantizer = make_outlier_quantizer(frame.quality);
        for (size_t i = 0; i < frame.outliers.size(); i++) {
            const SZ3::SPERR::Outlier &outlier = frame.outliers[i];
            if (outlier.pos >= outlier_map.size()) {
                throw std::runtime_error("SPERR outlier index out of bounds while building outlier map.");
            }
            double err = outlier.err;
            const int64_t q = quantizer.quantize_and_overwrite(err, 0.0);
            outlier_map[outlier.pos] = q;
        }

        return outlier_map;
    }

    static void apply_outlier_map(const Config &conf, const std::vector<int64_t> &bins, std::vector<double> &values) {
        validate_outlier_map(conf, bins, values.size());

        const SPERRQualityControl qc = resolve_quality_control(conf);
        if (qc.mode != SPERRMode::PWE_MODE) {
            return;
        }

        ScalarQuantizer<double, int64_t> quantizer = make_outlier_quantizer(qc.quality);
        for (size_t i = 0; i < bins.size(); i++) {
            if (bins[i] == 0) {
                continue;
            }
            values[i] += quantizer.recover(0.0, bins[i]);
        }
    }

    static std::vector<uchar> encode_coefficient_stream(const SPERRFrame &frame) {
        std::vector<uchar> stream;
        stream.insert(stream.end(), frame.conditioner_header.begin(), frame.conditioner_header.end());

        if (frame.constant_field) {
            return stream;
        }

        SPERREncoder<int64_t, 3> coeff_encoder(frame.dims);
        std::vector<uchar> speck_stream = coeff_encoder.encode_stream(frame.coeff_bins);
        stream.insert(stream.end(), speck_stream.begin(), speck_stream.end());
        return stream;
    }

    static SPERRFrame decode_coefficient_stream(const SPERR3DLayout &layout, const uchar *cmpData, size_t cmpSize) {
        SPERRFrame frame;
        frame.dims = {layout.dimx, layout.dimy, layout.dimz};

        if (cmpSize < frame.conditioner_header.size()) {
            throw std::runtime_error("SPERR coefficient stream too short for conditioner header.");
        }

        std::copy_n(cmpData, frame.conditioner_header.size(), frame.conditioner_header.begin());
        const SZ3::SPERR::Conditioner conditioner;
        frame.constant_field = conditioner.is_constant(frame.conditioner_header[0]);

        if (frame.constant_field) {
            if (cmpSize != frame.conditioner_header.size()) {
                throw std::runtime_error("SPERR constant-field coefficient stream has trailing payload.");
            }
            return frame;
        }

        frame.q = conditioner.retrieve_q(frame.conditioner_header);
        if (!(frame.q > 0.0)) {
            throw std::runtime_error("SPERR invalid quantization step in coefficient stream.");
        }

        const uchar *pos = cmpData + frame.conditioner_header.size();
        const size_t remaining = cmpSize - frame.conditioner_header.size();
        if (remaining < SZ3::SPERR::SPECK_INT<uint8_t>::header_size) {
            throw std::runtime_error("SPERR coefficient stream missing SPECK payload.");
        }

        SPERREncoder<int64_t, 3> coeff_encoder(frame.dims);
        size_t consumed = 0;
        frame.coeff_bins = coeff_encoder.decode_stream(pos, remaining, frame.dims[0] * frame.dims[1] * frame.dims[2],
                                                       consumed);
        if (consumed != remaining) {
            throw std::runtime_error("SPERR coefficient stream has trailing payload.");
        }

        return frame;
    }

    static double estimate_mse_midtread(const std::vector<double> &vals, double q) {
        const size_t len = vals.size();
        const size_t stride_size = 4096;
        const size_t num_strides = len / stride_size;
        std::vector<double> accum(num_strides + 1, 0.0);

        for (size_t i = 0; i < num_strides; i++) {
            const std::vector<double>::const_iterator begin = vals.cbegin() + i * stride_size;
            accum[i] = std::accumulate(begin, begin + stride_size, 0.0, [q](double init, double v) {
                const double diff = std::remainder(v, q);
                return init + diff * diff;
            });
        }

        accum[num_strides] =
            std::accumulate(vals.cbegin() + num_strides * stride_size, vals.cend(), 0.0, [q](double init, double v) {
                const double diff = std::remainder(v, q);
                return init + diff * diff;
            });

        const double total = std::accumulate(accum.cbegin(), accum.cend(), 0.0);
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

        const double target_mse =
            (frame.conditioned_range * frame.conditioned_range) * std::pow(10.0, -frame.quality / 10.0);
        double q = 2.0 * std::sqrt(target_mse * 3.0);
        while (estimate_mse_midtread(frame.wavelet_coeffs, q) > target_mse) {
            q /= std::exp2(0.25);
        }
        return q;
    }

    std::vector<uchar> coeff_stream_;
};

}  // namespace SZ3

#endif
