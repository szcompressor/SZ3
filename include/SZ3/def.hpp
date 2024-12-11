#ifndef _DEF_HPP
#define _DEF_HPP

namespace SZ3 {

typedef unsigned int uint;
typedef unsigned char uchar;
#define SZ_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH \
    "The buffer for compressed data is not large enough."
}  // namespace SZ3

#ifdef _MSC_VER
#define ALWAYS_INLINE __forceinline
#else
#define ALWAYS_INLINE __attribute__((always_inline))
#endif

#endif
