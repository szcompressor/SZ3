#ifndef SZ3_DECOMPOSITION_INTERFACE
#define SZ3_DECOMPOSITION_INTERFACE

#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"

namespace SZ3::concepts {

/**
 * @brief Interface for Decomposition modules
 * 
 * Decomposition represents the lossy part of the compression pipeline, converting original data
 * to integers (e.g., prediction + quantization).
 * 
 * @tparam Ti Original data type
 * @tparam To Decomposed data type (usually int)
 * @tparam N Original data dimension
 */
template <class Ti, class To, uint N, class ToAllocator = std::allocator<To>>
class DecompositionInterface {
   public:
    virtual ~DecompositionInterface() = default;

    /**
     * @brief Predict the data and quantize the error
     * 
     * @param conf Compression configuration
     * @param data Original input data
     * @return std::vector<To, ToAllocator> Quantized prediction errors
     */
    virtual std::vector<To, ToAllocator> compress(const Config &conf, Ti *data) = 0;

    /**
     * @brief Reconstruct the data from quantized indices
     * 
     * @param conf Compression configuration
     * @param quant_inds Quantized prediction errors
     * @param dec_data Buffer to write the reconstructed data
     * @return Ti* Pointer to the reconstructed data (same as dec_data)
     */
    virtual Ti *decompress(const Config &conf, std::vector<To, ToAllocator> &quant_inds, Ti *dec_data) = 0;

    /**
     * @brief Serialize the decomposition module and store it to a buffer
     * 
     * @param c Reference to the buffer pointer. It will be advanced to the next empty location.
     */
    virtual void save(uchar *&c) = 0;

    /**
     * @brief Deserialize the decomposition module from a buffer
     * 
     * @param c Reference to the buffer pointer. It will be advanced after reading.
     * @param remaining_length Remaining length of the buffer
     */
    virtual void load(const uchar *&c, size_t &remaining_length) = 0;

    virtual size_t size_est() { return 0; }

    virtual std::pair<To, To> get_out_range() = 0;

    virtual void print() {}

    //        virtual void clear() {};
};

}  // namespace SZ3::concepts

#endif
