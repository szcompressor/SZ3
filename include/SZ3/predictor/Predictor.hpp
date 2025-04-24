#ifndef SZ3_PREDICTOR_HPP
#define SZ3_PREDICTOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/utils/BlockwiseIterator.hpp"

namespace SZ3::concepts {

/**
 * Data prediction interface
 * BlockwiseDecomposition will automatically iterate through multidimensional data to apply the
 * Predictor on each data point.
 * @tparam T original data type
 * @tparam N original data dimension
 *
 */
template <class T, uint N>
class PredictorInterface {
   public:
    using block_iter = typename block_data<T, N>::block_iterator;

    virtual ~PredictorInterface() = default;

    /**
     * predict the value for a single data point
     * @param block_iter the iterator of the block
     * @param T* the pointer to the single data point
     * @param std::array<size_t, N> the relative index of the data point in the block
     * @return the predicted value
     */
    ALWAYS_INLINE virtual T predict(const block_iter &, T *, const std::array<size_t, N> &) = 0;

    /**
     * estimate the prediction error ( |prediction value - read value|)  for a single data point
     * @param iter the iterator of the single data point
     * @return the estimated prediction error
     */
    ALWAYS_INLINE virtual T estimate_error(const block_iter &, T *, const std::array<size_t, N> &) = 0;

    /**
     * compute auxiliary info (e.g., coefficients) for the given data block
     * @param block_iter of the block
     * @return whether the predictor is suitable for the block (e.g., data with 100x1 shape is not suitable for 2D
     * regression)
     */
    virtual bool precompress(const block_iter &) = 0;

    /**
     * store the auxiliary info (e.g., coefficients) to this class's internal storage
     */
    virtual void precompress_block_commit() = 0;

    virtual bool predecompress(const block_iter &) = 0;

    /**
     * serialize the predictor and store it to a buffer
     * @param c One large buffer is pre-allocated, and the start location of the serialized predictor in the buffer is
     * indicated by c. After saving the predictor to the buffer, this function should change c to indicate the next
     * empty location in the buffer
     */
    virtual void save(uchar *&c) = 0;

    /**
     * deserialize the predictor from a buffer
     * @param c start location of the predictor in the buffer
     * @param remaining_length the remaining length of the buffer
     */
    virtual void load(const uchar *&c, size_t &remaining_length) = 0;

    virtual size_t get_padding() { return 0; }

    virtual void print() const = 0;
};

}  // namespace SZ3::concepts

#endif
