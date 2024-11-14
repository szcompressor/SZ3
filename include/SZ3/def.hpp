#ifndef _DEF_HPP
#define _DEF_HPP

namespace SZ3 {

typedef unsigned int uint;
typedef unsigned char uchar;
#define SZ_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH \
    "The buffer for compressed data is not large enough. Ideally, set it at least 2X original data size."
}  // namespace SZ3

#endif
