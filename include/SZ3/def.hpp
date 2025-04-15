#ifndef _DEF_HPP
#define _DEF_HPP

#include <cmath>

namespace SZ {
#define force_inline __attribute__((always_inline)) inline
    typedef unsigned int uint;
    typedef unsigned char uchar;
}

#ifdef _MSC_VER
#define ALWAYS_INLINE __forceinline
#else
#define ALWAYS_INLINE __attribute__((always_inline))
#endif


#endif
