#ifndef SZ3_DEF_HPP
#define SZ3_DEF_HPP

namespace SZ3 {

typedef unsigned int uint;
typedef unsigned char uchar;
#define SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH \
    "The buffer for compressed data is not large enough."
}  // namespace SZ3

#ifdef _MSC_VER
#define ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

#ifdef __has_builtin
#if __has_builtin(__builtin_assoc_barrier)
#define SZ3_HAS_ASSOC_BARRIER
#endif
#endif

namespace SZ3 {
// Builds differ in whether they fuse a * b + c into one fma, and a decoder must reproduce its encoder's
// rounding bit for bit.
template <class T>
ALWAYS_INLINE T rounded(T x) {
#ifdef SZ3_HAS_ASSOC_BARRIER
    return __builtin_assoc_barrier(x);
#else
    volatile T r = x;
    return r;
#endif
}
}  // namespace SZ3

#endif
