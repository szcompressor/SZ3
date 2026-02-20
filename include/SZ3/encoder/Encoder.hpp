/**
 * @file Encoder.hpp
 * @ingroup Encoder
 */

#ifndef SZ3_ENCODER_HPP
#define SZ3_ENCODER_HPP

#include <vector>

#include "SZ3/def.hpp"
namespace SZ3 {

namespace concepts {

/**
 * @brief Interface for encoders
 * 
 * Encoders transform the input data into a more compact representation (usually lossless).
 * Examples: Huffman encoding, Run-length encoding.
 * 
 * @tparam T Input data type (usually int)
 * @tparam TAllocator Allocator type
 */
template <class T, class TAllocator = std::allocator<T>>
class EncoderInterface {
   public:
    virtual ~EncoderInterface() = default;

    /**
     * @brief Initialize the encoder
     * 
     * Perform pre-computation tasks, such as building a Huffman tree.
     * 
     * @param bins Vector of to-be-encoded integers
     * @param stateNum Expected range of values [0, stateNum). If 0, no range guarantee.
     */
    virtual void preprocess_encode(const std::vector<T, TAllocator> &bins, int stateNum) = 0;

    /**
     * @brief Encode the input vector into a byte stream
     * 
     * @param bins Input vector
     * @param bytes Output byte stream (pointer reference will be updated)
     * @return size_t Size of the output in bytes
     */
    virtual size_t encode(const std::vector<T, TAllocator> &bins, uchar *&bytes) = 0;

    /**
     * @brief Decode a byte stream into a vector
     * 
     * @param bytes Input byte stream (pointer reference will be updated)
     * @param targetLength Number of elements to decode
     * @return std::vector<T, TAllocator> Decoded vector
     */
    virtual std::vector<T, TAllocator> decode(const uchar *&bytes, size_t targetLength) = 0;

    /**
     * @brief Serialize the encoder and store it to a buffer
     * 
     * @param c Reference to the buffer pointer. It will be advanced to the next empty location.
     */
    virtual void save(uchar *&c) = 0;

    /**
     * @brief Deserialize the encoder from a buffer
     * 
     * @param c Reference to the buffer pointer. It will be advanced after reading.
     * @param remaining_length Remaining length of the buffer
     */
    virtual void load(const uchar *&c, size_t &remaining_length) = 0;

    virtual void postprocess_decode() = 0;

    virtual void postprocess_encode() = 0;

    virtual void preprocess_decode() = 0;

    // return the size of the encoder itself (such as the tree size of the huffman encoder)
    virtual size_t size_est() { return 0; }
};
}  // namespace concepts
}  // namespace SZ3
#endif
