/**
 * @file SPERRDecomposition.hpp
 * @ingroup Decomposition
 */

#ifndef SZ3_SZ3_SPERR_DECOMPOSITION_HPP
#define SZ3_SZ3_SPERR_DECOMPOSITION_HPP

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/preprocessor/SPERRTransform.hpp"
#include "SZ3/quantizer/OutlierQuantizer.hpp"
#include "SZ3/quantizer/ScalarQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/thirdparty/sperr/SPERRTypes.hpp"

namespace SZ3 {

/**
 * @brief SPERR decomposition assembled from separately reusable modules.
 *
 * This is an alternative factoring of `SPERRFusedDecomposition`, not a replacement for it.
 * `SPERRFusedDecomposition` fuses transform, quantization, SPECK entropy coding and outlier
 * handling into one stage and emits only the outlier map, which leaves the pipeline's
 * encoder slot with almost nothing to do. This module keeps the four concerns apart:
 *
 *   - transform  : `SPERRTransform<T, N>`               (conditioner + CDF9/7 wavelet)
 *   - quantizer  : `ScalarQuantizer<double, int64_t>`   (mid-tread, step `q`)
 *   - outliers   : `OutlierQuantizer<double, int64_t>`  (sparse corrections, carried in `save()`)
 *   - encoder    : *not owned here* -- the quantized wavelet coefficients are returned
 *                  from `compress()` and travel to whatever `EncoderInterface` the
 *                  pipeline was wired with (`SPERREncoder` for SPECK, `HuffmanEncoder`,
 *                  `BypassEncoder`, ...).
 *
 * The arithmetic is identical to `SPERRFusedDecomposition`, so reconstructions agree
 * bit-for-bit for the same input and configuration. What changes is *where* the bytes
 * go: the coefficient stream leaves through the encoder rather than through the
 * decomposition's own `save()`, and the (sparse) outlier corrections are held in
 * module state exactly the way `LinearQuantizer` holds its `unpred` list.
 *
 * SPERR is 3D floating-point only.
 *
 * @tparam T Original data type (`float` / `double`)
 * @tparam N Data dimension; only 3 is supported
 */
template <class T, uint N>
class SPERRDecomposition : public concepts::DecompositionInterface<T, int64_t, N> {
   public:
    SPERRDecomposition() = default;

    /**
     * @brief Transform, quantize, and collect outliers.
     * @param conf Compression configuration (3D, ABS or PSNR error bound)
     * @param data Input field
     * @return Quantized wavelet coefficients, one bin per element, for the encoder stage
     */
    std::vector<int64_t> compress(const Config &conf, T *data) override {
        transform_.set_dims(SPERRTransform<T, N>::dims_from_config(conf));

        const QualityControl qc = resolve_quality_control(conf);
        outlier_quantizer_ = OutlierQuantizer<double, int64_t>(qc.quality);

        std::vector<double> conditioned = transform_.condition(data, conf.num);
        if (transform_.constant_field()) {
            // Nothing to transform or quantize: the header alone rebuilds the field.
            // A full-length zero stream keeps the encoder stage's contract intact.
            q_ = 0.0;
            coeff_quantizer_ = ScalarQuantizer<double, int64_t>(1.0);
            return std::vector<int64_t>(conf.num, int64_t{0});
        }

        const std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> minmax =
            std::minmax_element(conditioned.cbegin(), conditioned.cend());
        const double conditioned_range = *minmax.second - *minmax.first;

        const bool pwe = (qc.mode == SPERRMode::PWE_MODE);
        std::vector<double> coefficients =
            pwe ? transform_.forward_wavelet(conditioned) : transform_.forward_wavelet(std::move(conditioned));

        q_ = estimate_q(qc, conditioned_range, coefficients);
        if (!(q_ > 0.0)) {
            throw std::invalid_argument("SPERR quantization step must be positive.");
        }
        coeff_quantizer_ = ScalarQuantizer<double, int64_t>(q_);

        std::fesetround(FE_TONEAREST);

        std::vector<int64_t> coeff_bins(coefficients.size(), int64_t{0});
        std::vector<double> reconstructed_coeffs;
        if (pwe) {
            reconstructed_coeffs.resize(coefficients.size(), 0.0);
        }
        for (size_t i = 0; i < coefficients.size(); i++) {
            double coeff = coefficients[i];
            coeff_bins[i] = coeff_quantizer_.quantize_and_overwrite(coeff, 0.0);
            if (pwe) {
                reconstructed_coeffs[i] = coeff;
            }
        }

        outlier_quantizer_.clear();
        if (pwe) {
            // Replay the decoder's inverse wavelet so outliers are measured against the
            // values the decoder will actually see, in the conditioned domain.
            const std::vector<double> reconstructed = transform_.inverse_wavelet(std::move(reconstructed_coeffs));
            outlier_quantizer_.collect(conditioned, reconstructed);
        }

        return coeff_bins;
    }

    /**
     * @brief Dequantize, inverse-transform, and re-apply outlier corrections.
     * @param conf Compression configuration
     * @param coeff_bins Quantized wavelet coefficients from the encoder stage
     * @param dec_data Output buffer
     */
    T *decompress(const Config &conf, std::vector<int64_t> &coeff_bins, T *dec_data) override {
        transform_.set_dims(SPERRTransform<T, N>::dims_from_config(conf));

        std::vector<double> values;
        if (transform_.constant_field()) {
            if (std::any_of(coeff_bins.begin(), coeff_bins.end(), [](int64_t bin) { return bin != 0; })) {
                throw std::runtime_error("SPERR constant field should not carry coefficient payload.");
            }
        } else {
            if (coeff_bins.size() != conf.num) {
                throw std::runtime_error("SPERR coefficient bin count mismatch.");
            }
            values.resize(coeff_bins.size(), 0.0);
            for (size_t i = 0; i < coeff_bins.size(); i++) {
                values[i] = coeff_quantizer_.recover(0.0, coeff_bins[i]);
            }
            values = transform_.inverse_wavelet(std::move(values));
            outlier_quantizer_.apply(values);
        }

        // Constant fields are expanded here; otherwise this just adds the mean back.
        transform_.inverse_condition(values);

        if (values.size() != conf.num) {
            throw std::runtime_error("SPERR reconstruction length mismatch.");
        }
        for (size_t i = 0; i < conf.num; i++) {
            dec_data[i] = static_cast<T>(values[i]);
        }
        return dec_data;
    }

    /**
     * @brief Serialize the three owned modules, in pipeline order.
     * @param c Buffer pointer; advanced past the written bytes.
     */
    void save(uchar *&c) override {
        transform_.save(c);
        coeff_quantizer_.save(c);
        outlier_quantizer_.save(c);
    }

    /**
     * @brief Deserialize the three owned modules.
     * @param c Buffer pointer; advanced past the consumed bytes.
     * @param remaining_length Remaining buffer length; decremented.
     */
    void load(const uchar *&c, size_t &remaining_length) override {
        transform_.load(c, remaining_length);
        coeff_quantizer_.load(c, remaining_length);
        outlier_quantizer_.load(c, remaining_length);
    }

    size_t size_est() override { return transform_.size_est() + kScalarQuantizerBytes + outlier_quantizer_.size_est(); }

    /**
     * @brief Output range advertised to `SZGenericCompressor`.
     *
     * The coefficient bins are genuinely signed (SPECK codes sign and magnitude
     * separately), so the lower end is nominal. `SPERRFusedDecomposition` reports the same
     * pair for the same reason; shifting the bins to be non-negative would defeat the
     * encoder's sign/magnitude split.
     */
    std::pair<int64_t, int64_t> get_out_range() override {
        // SPERR bins are a SPECK bitstream, not a quantizer domain; 0 means "no bin range".
        return std::make_pair<int64_t, int64_t>(0, 0);
    }

    void print() override {
        printf("[SPERRDecomposition] q=%.8G, constant_field=%d, outliers=%zu\n", q_,
               static_cast<int>(transform_.constant_field()), outlier_quantizer_.size());
    }

    /// The conditioner + wavelet stage.
    SPERRTransform<T, N> &get_transform() { return transform_; }

    /// The mid-tread quantizer applied to wavelet coefficients.
    ScalarQuantizer<double, int64_t> &get_coefficient_quantizer() { return coeff_quantizer_; }

    /// The sparse outlier-correction stage.
    OutlierQuantizer<double, int64_t> &get_outlier_quantizer() { return outlier_quantizer_; }

    /// Quantization step chosen for the wavelet coefficients (0 for a constant field).
    double get_q() const { return q_; }

   private:
    /// Error-bound mode plus target, resolved from `Config`.
    struct QualityControl {
        SPERRMode mode = SPERRMode::PWE_MODE;
        double quality = 0.0;
    };

    static QualityControl resolve_quality_control(const Config &conf) {
        QualityControl qc;
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

    /**
     * @brief Rate-control policy: pick the coefficient quantization step.
     *
     * Kept here rather than in any of the four modules because it is a fifth concern --
     * it needs the error-bound mode, the conditioned data range and the coefficients at
     * once. Mirrors `SPERRFusedDecomposition`'s private policy exactly.
     */
    static double estimate_q(const QualityControl &qc, double conditioned_range,
                             const std::vector<double> &coefficients) {
        if (qc.mode == SPERRMode::PWE_MODE) {
            return qc.quality * 1.5;
        }
        if (conditioned_range <= 0.0) {
            throw std::invalid_argument("SPERR PSNR mode requires positive conditioned data range.");
        }
        if (coefficients.empty()) {
            throw std::invalid_argument("SPERR q-estimation requires wavelet coefficients.");
        }

        const double target_mse = (conditioned_range * conditioned_range) * std::pow(10.0, -qc.quality / 10.0);
        double q = 2.0 * std::sqrt(target_mse * 3.0);
        while (estimate_mse_midtread(coefficients, q) > target_mse) {
            q /= std::exp2(0.25);
        }
        return q;
    }

    /// `ScalarQuantizer::save()` writes a uid byte plus three doubles.
    static constexpr size_t kScalarQuantizerBytes = sizeof(uchar) + 3 * sizeof(double);

    SPERRTransform<T, N> transform_;
    ScalarQuantizer<double, int64_t> coeff_quantizer_;
    OutlierQuantizer<double, int64_t> outlier_quantizer_;
    double q_ = 0.0;
};

}  // namespace SZ3

#endif  // SZ3_SZ3_SPERR_DECOMPOSITION_HPP
