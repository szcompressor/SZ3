#ifndef SZ3_PREDICTOR_HPP
#define SZ3_PREDICTOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/utils/BlockwiseIterator.hpp"

namespace SZ3::concepts {

/**
 * @brief Interface for scalar data prediction, used exclusively by `BlockwiseDecomposition`.
 * 
 * A `Predictor` estimates a data value based on its neighboring points. `BlockwiseDecomposition`
 * uses this interface to iterate over data blocks and apply prediction + quantization on each point.
 * 
 * **Note:** This interface applies only when using `BlockwiseDecomposition`. Other decomposition
 * strategies (e.g., `InterpolationDecomposition`, `ZFPDecomposition`) do not use this interface.
 * If your prediction method operates on full data arrays rather than scalar values, implement
 * `DecompositionInterface` directly instead.
 * 
 * @tparam T Original data type
 * @tparam N Original data dimension
 */
template <class T, uint N>
class PredictorInterface {
   public:
    using block_iter = typename block_data<T, N>::block_iterator;

    virtual ~PredictorInterface() = default;

    /**
     * @brief Predict the value for a single data point
     * 
     * @param block Iterator to the current data block
     * @param data Pointer to the single data point
     * @param index Relative index of the data point in the block
     * @return T The predicted value
     */
    ALWAYS_INLINE virtual T predict(const block_iter &block, T *data, const std::array<size_t, N> &index) = 0;

    /**
     * @brief Estimate the prediction error for a single data point
     * 
     * Error = |prediction value - read value|
     * 
     * @param block Iterator to the current data block
     * @param data Pointer to the single data point
     * @param index Relative index of the data point in the block
     * @return T The estimated prediction error
     */
    ALWAYS_INLINE virtual T estimate_error(const block_iter &block, T *data, const std::array<size_t, N> &index) = 0;

    /**
     * @brief Compute auxiliary info (e.g., coefficients) for the given data block
     * 
     * @param block Iterator to the current data block
     * @return true If the predictor is suitable for the block
     * @return false If the predictor is not suitable (e.g., shape mismatch)
     */
    virtual bool precompress(const block_iter &block) = 0;

    /**
     * @brief Store the auxiliary info (e.g., coefficients) to internal storage
     */
    virtual void precompress_block_commit() = 0;

    virtual bool predecompress(const block_iter &) = 0;

    /**
     * @brief Serialize the predictor and store it to a buffer
     * 
     * @param c Reference to the buffer pointer. It will be advanced to the next empty location after writing.
     */
    virtual void save(uchar *&c) = 0;

    /**
     * @brief Deserialize the predictor from a buffer
     * 
     * @param c Reference to the buffer pointer. It will be advanced after reading.
     * @param remaining_length Remaining length of the buffer
     */
    virtual void load(const uchar *&c, size_t &remaining_length) = 0;

    virtual size_t get_padding() { return 0; }

    virtual void print() const = 0;
};

}  // namespace SZ3::concepts

#endif
