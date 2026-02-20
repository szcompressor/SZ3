#ifndef SZ3_QUANTIZER_HPP
#define SZ3_QUANTIZER_HPP

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
    ALWAYS_INLINE virtual To quantize_and_overwrite(Ti &data, Ti pred) = 0;

    /**
     * @brief Reconstruct the data point from quantized index
     * 
     * @param pred Predicted value
     * @param quant_index Quantized index
     * @return Ti Reconstructed value
     */
    ALWAYS_INLINE virtual Ti recover(Ti pred, To quant_index) = 0;

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

#endif
