/**
 * @file SPERRTransform.hpp
 * @ingroup Preprocessor
 */

#ifndef SZ3_SPERR_TRANSFORM_HPP
#define SZ3_SPERR_TRANSFORM_HPP

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/preprocessor/PreProcessor.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/thirdparty/sperr/SPERRTypes.hpp"

namespace SZ3 {

/**
 * @brief SPERR conditioner + CDF9/7 wavelet transform, as a standalone preprocessor.
 *
 * The transform half of the SPERR pipeline, separated from quantization, entropy coding and
 * outlier handling so it can be reused with any quantizer/encoder pair. `SPERRFusedDecomposition`
 * keeps the fused version.
 *
 * The forward direction is two stages, both exposed individually because outlier detection needs
 * the intermediate (conditioned, pre-wavelet) values:
 *   1. `condition()` -- subtract the mean / detect a constant field.
 *   2. `forward_wavelet()` -- 3D CDF9/7 discrete wavelet transform.
 *
 * The inverse must stay separable for the same reason: outlier corrections are applied in the
 * conditioned domain, i.e. between `inverse_wavelet()` and `inverse_condition()`.
 *
 * All intermediate values are `double`, as SPERR's kernels are, regardless of `T`.
 *
 * @tparam T Original (floating-point) data type
 * @tparam N Data dimension; only 3 is supported by SPERR
 */
template <class T, uint N>
class SPERRTransform : public concepts::PreprocessorInterface<T, N> {
   public:
    SPERRTransform() = default;

    /**
     * @brief Construct with an explicit SPERR-ordered shape.
     * @param dims Shape in SPERR ordering `[x, y, z]`
     */
    explicit SPERRTransform(const SZ3::SPERR::dims_type &dims) { set_dims(dims); }

    /**
     * @brief Translate an SZ3 `Config` (dims ordered `[z, y, x]`) into SPERR ordering `[x, y, z]`.
     * @throw std::invalid_argument if the configuration is not 3D
     */
    static SZ3::SPERR::dims_type dims_from_config(const Config &conf) {
        if (conf.N != 3 || conf.dims.size() != 3) {
            throw std::invalid_argument("SPERRTransform supports 3D data only.");
        }
        return SZ3::SPERR::dims_type{conf.dims[2], conf.dims[1], conf.dims[0]};
    }

    /// Set the SPERR-ordered `[x, y, z]` shape this transform operates on.
    void set_dims(const SZ3::SPERR::dims_type &dims) {
        if (N != 3) {
            throw std::invalid_argument("SPERRTransform supports 3D data only.");
        }
        if (dims[0] == 0 || dims[1] == 0 || dims[2] == 0) {
            throw std::invalid_argument("SPERRTransform requires non-zero dimensions.");
        }
        dims_ = dims;
    }

    /// SPERR-ordered `[x, y, z]` shape.
    const SZ3::SPERR::dims_type &get_dims() const { return dims_; }

    /// Number of elements implied by the configured shape.
    size_t num_elements() const { return dims_[0] * dims_[1] * dims_[2]; }

    /**
     * @brief Forward stage 1: condition the field (subtract mean, detect constant field).
     *
     * Updates the internal conditioner header and the constant-field flag.
     *
     * @param data Input field, `num` elements in SZ3 row-major order
     * @param num Number of elements; must match the configured shape
     * @return Conditioned values. Meaningless (but returned unchanged) for a constant field.
     */
    std::vector<double> condition(const T *data, size_t num) {
        validate_scalar_type();
        validate_length(num);

        std::vector<double> values(num);
        std::copy(data, data + num, values.begin());

        SZ3::SPERR::Conditioner conditioner;
        header_ = conditioner.condition(values, dims_);
        constant_field_ = conditioner.is_constant(header_[0]);
        return values;
    }

    /**
     * @brief Forward stage 2: 3D CDF9/7 discrete wavelet transform.
     * @param values Conditioned values (taken by value so callers may `std::move` into it)
     * @return Wavelet coefficients, same length as the input
     */
    std::vector<double> forward_wavelet(std::vector<double> values) const {
        validate_length(values.size());

        SZ3::SPERR::CDF97 cdf;
        if (cdf.take_data(std::move(values), dims_) != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR forward wavelet setup failed.");
        }
        cdf.dwt3d();
        return cdf.release_data();
    }

    /**
     * @brief Full forward transform: condition, then wavelet.
     * @return Wavelet coefficients, or an empty vector if the field is constant.
     */
    std::vector<double> preprocess(const T *data, size_t num) {
        std::vector<double> values = condition(data, num);
        if (constant_field_) {
            return std::vector<double>();
        }
        return forward_wavelet(std::move(values));
    }

    /**
     * @brief Inverse stage 1: inverse 3D CDF9/7 discrete wavelet transform.
     * @param coefficients Wavelet coefficients (taken by value so callers may `std::move` into it)
     * @return Values in the conditioned domain
     */
    std::vector<double> inverse_wavelet(std::vector<double> coefficients) const {
        validate_length(coefficients.size());

        SZ3::SPERR::CDF97 cdf;
        if (cdf.take_data(std::move(coefficients), dims_) != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR inverse wavelet setup failed.");
        }
        cdf.idwt3d();
        return cdf.release_data();
    }

    /**
     * @brief Inverse stage 2: undo conditioning (add the mean back, or expand a constant field).
     *
     * For a constant field `values` may be empty on entry; it is resized and filled from
     * the header. Otherwise `values` must already hold the conditioned-domain reconstruction.
     */
    void inverse_condition(std::vector<double> &values) const {
        SZ3::SPERR::Conditioner conditioner;
        if (conditioner.inverse_condition(values, dims_, header_) != SZ3::SPERR::RTNType::Good) {
            throw std::runtime_error("SPERR inverse conditioning failed.");
        }
    }

    /**
     * @brief Full inverse transform: inverse wavelet, inverse conditioning, cast back to `T`.
     * @param coefficients Wavelet coefficients; ignored (and expected empty) for a constant field
     * @param dec_data Output buffer of at least `num` elements
     * @param num Expected output length
     */
    void postprocess(std::vector<double> coefficients, T *dec_data, size_t num) const {
        validate_scalar_type();
        validate_length(num);

        std::vector<double> values;
        if (!constant_field_) {
            values = inverse_wavelet(std::move(coefficients));
        }
        inverse_condition(values);

        if (values.size() != num) {
            throw std::runtime_error("SPERR reconstruction length mismatch.");
        }
        for (size_t i = 0; i < num; i++) {
            dec_data[i] = static_cast<T>(values[i]);
        }
    }

    /// True when `condition()`/`load()` detected a constant field (no wavelet/quantization needed).
    bool constant_field() const { return constant_field_; }

    /// Raw 17-byte SPERR conditioner header.
    const SZ3::SPERR::condi_type &header() const { return header_; }

    /// Install a conditioner header (e.g. one produced by another instance).
    void set_header(const SZ3::SPERR::condi_type &header) {
        header_ = header;
        const SZ3::SPERR::Conditioner conditioner;
        constant_field_ = conditioner.is_constant(header_[0]);
    }

    /**
     * @brief Serialize the transform state (the conditioner header).
     * @param c Buffer pointer; advanced past the written bytes.
     */
    void save(uchar *&c) const { write(header_.data(), header_.size(), c); }

    /**
     * @brief Deserialize the transform state.
     * @param c Buffer pointer; advanced past the consumed bytes.
     * @param remaining_length Remaining buffer length; decremented.
     */
    void load(const uchar *&c, size_t &remaining_length) {
        read(header_.data(), header_.size(), c, remaining_length);
        const SZ3::SPERR::Conditioner conditioner;
        constant_field_ = conditioner.is_constant(header_[0]);
    }

    /// Serialized size of `save()`.
    size_t size_est() const { return header_.size(); }

   private:
    static void validate_scalar_type() {
        if (!std::is_floating_point<T>::value) {
            throw std::invalid_argument("SPERRTransform supports floating-point data only.");
        }
    }

    void validate_length(size_t num) const {
        if (num_elements() == 0) {
            throw std::invalid_argument("SPERRTransform requires dimensions to be set.");
        }
        if (num != num_elements()) {
            throw std::invalid_argument("SPERRTransform length does not match configured dimensions.");
        }
    }

    SZ3::SPERR::dims_type dims_ = {0, 0, 0};
    SZ3::SPERR::condi_type header_{};
    bool constant_field_ = false;
};

}  // namespace SZ3

#endif  // SZ3_SPERR_TRANSFORM_HPP
