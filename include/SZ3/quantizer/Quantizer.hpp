#ifndef SZ3_QUANTIZER_HPP
#define SZ3_QUANTIZER_HPP

#include <SZ3/def.hpp>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace SZ3::concepts {

/**
 * quantize a single data point with precision loss
 * E.g: float->int
 * @tparam Ti original data type
 * @tparam To quantized data type
 */
template <class Ti, class To>
class   QuantizerInterface {
   public:
    virtual ~QuantizerInterface() = default;

    /**
     * quantize the error (error=data-pred) based on error bound, and overwrite the data with reconstructed value
     * @param data single data point
     * @param pred predicted value for this data point
     * @return quantized error
     */
    virtual To quantize_and_overwrite(Ti &data, Ti pred) = 0;

    /**
     * reconstructed the data point
     * @param pred predicted value for the data point
     * @param quant_index quantized error
     * @return reconstructed value of the data point
     */
    virtual Ti recover(Ti pred, To quant_index) = 0;

    virtual To force_save_unpred(Ti ori) = 0;

    /**
     ** serialize the quantizer and store it to a buffer
     * @param c One large buffer is pre-allocated, and the start location of the serialized quantizer in the buffer is
     *indicated by c. After saving the quantizer to the buffer, this function should change c to indicate the next empty
     *location in the buffer
     */
    virtual void save(uchar *&c) const = 0;

    /**
     * deserialize the quantizer from a buffer
     * @param c start location of the quantizer in the buffer
     * @param remaining_length the remaining length of the buffer
     */
    virtual void load(const uchar *&c, size_t &remaining_length) = 0;

    virtual std::pair<To, To> get_out_range() const = 0;

    virtual void precompress_data() {}

    virtual void predecompress_data() {}

    /**
     * this function is always executed before save()
     * DO NOT reset non-temporary variables (such as unpredictable_data) in this function.
     */
    virtual void postcompress_data() {}

    /**
     * DO NOT reset non-temporary variables (such as unpredictable_data) in this function.
     */
    virtual void postdecompress_data() {}

    virtual void print() {}
};
}  // namespace SZ3::concepts

namespace SZ3 {
/// Detects the optional (non-virtual) `size_est()` some quantizers expose.
template <class Q, class = void>
struct quantizer_has_size_est : std::false_type {};
template <class Q>
struct quantizer_has_size_est<Q, std::void_t<decltype(std::declval<Q &>().size_est())>> : std::true_type {};

/// The quantizer's serialized-size estimate, or 0 when it does not expose one.
template <class Q>
size_t quantizer_size_est(Q &q) {
    if constexpr (quantizer_has_size_est<Q>::value) {
        return q.size_est();
    } else {
        return 0;
    }
}
}  // namespace SZ3

#endif
