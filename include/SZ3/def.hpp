#ifndef SZ3_DEF_HPP
#define SZ3_DEF_HPP

namespace SZ3 {

/**
 * @brief Unsigned integer type definition
 */
typedef unsigned int uint;

/**
 * @brief Unsigned char type definition
 */
typedef unsigned char uchar;

/**
 * @brief Error message for insufficient buffer size
 */
#define SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH \
    "The buffer for compressed data is not large enough."

#ifdef _MSC_VER
#define ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

// In nofma(a * b) + c the product is rounded before the add, so the compiler cannot merge them into one
// FMA instruction, which rounds only once; otherwise a build with FMA and one without decompress different values.
template <class T>
ALWAYS_INLINE T nofma(T x) {
    // Not __builtin_assoc_barrier: GCC 13 and 14 drop it when they vectorize the loop, and fuse the product anyway.
    volatile T r = x;
    return r;
}
}  // namespace SZ3

#endif
