#ifndef SZ3_QUANTIZER_HPP
#define SZ3_QUANTIZER_HPP

#include <cstddef>
#include <type_traits>
#include <utility>

#include "SZ3/def.hpp"

namespace SZ3::concepts {

/**
 * @brief Interface for quantization
 * 
 * Quantizers reduce the precision of data points, mapping them to a smaller set of values (integers).
 * E.g., float -> int.
 * 
 * @tparam Ti Original data type (Input)
 * @tparam To Quantized data type (Output)
 */
template <class Ti, class To>
class QuantizerInterface {
   public:
    virtual ~QuantizerInterface() = default;

    /**
     * @brief Quantize the prediction error and overwrite data with reconstructed value
     * 
     * @param data Single data point (input/output). Will be overwritten by reconstructed value.
     * @param pred Predicted value for this data point
     * @return To Quantized index/bin
     */
    virtual To quantize_and_overwrite(Ti &data, Ti pred) = 0;

    /**
     * @brief Reconstruct the data point from quantized index
     * 
     * @param pred Predicted value
     * @param quant_index Quantized index
     * @return Ti Reconstructed value
     */
    virtual Ti recover(Ti pred, To quant_index) = 0;

    virtual To force_save_unpred(Ti ori) = 0;

    /**
     * @brief Serialize the quantizer and store it to a buffer
     * 
     * @param c Reference to the buffer pointer. It will be advanced to the next empty location.
     */
    virtual void save(uchar *&c) const = 0;

    /**
     * @brief Deserialize the quantizer from a buffer
     * 
     * @param c Reference to the buffer pointer. It will be advanced after reading.
     * @param remaining_length Remaining length of the buffer
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

/**
 * @brief The bin type a Quantizer emits, deduced from `quantize_and_overwrite`.
 *
 * `concepts::QuantizerInterface<Ti, To>` does not publish `To` as a member typedef, so a
 * Decomposition that wants to forward its Quantizer's bin type cannot name it directly. This
 * deduces it from the only place it is observable - the return type of `quantize_and_overwrite`.
 *
 * Use it as a defaulted trailing template parameter so a Decomposition inherits the right
 * `DecompositionInterface<T, To, N>` instead of hard-coding `int`:
 *
 * @code
 * template <class T, uint N, class Quantizer, class To = quantizer_bin_t<T, Quantizer>>
 * class MyDecomposition : public concepts::DecompositionInterface<T, To, N> { ... };
 * @endcode
 *
 * For the `int`-bin quantizers this resolves to `int`, so existing wirings are unaffected.
 *
 * @tparam T         Data type the quantizer consumes
 * @tparam Quantizer A type implementing `concepts::QuantizerInterface<T, To>`
 */
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

template <class T, class Quantizer>
using quantizer_bin_t =
    decltype(std::declval<Quantizer &>().quantize_and_overwrite(std::declval<T &>(), std::declval<T>()));

}  // namespace SZ3

#endif
