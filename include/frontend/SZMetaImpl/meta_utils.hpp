#ifndef _utils_h
#define _utils_h

#include <cstdlib>

namespace SZMETA {

    template<typename T>
    inline void
    write_variable_to_dst(unsigned char *&dst, const T &var) {
        memcpy(dst, &var, sizeof(T));
        dst += sizeof(T);
    }

    template<typename T>
    inline void
    write_array_to_dst(unsigned char *&dst, const T *array, size_t length) {
        memcpy(dst, array, length * sizeof(T));
        dst += length * sizeof(T);
    }

    template<typename T>
    inline void
    read_variable_from_src(const unsigned char *&src, T &var) {
        memcpy(&var, src, sizeof(T));
        src += sizeof(T);
    }

    template<typename T>
    inline T *
    read_array_from_src(const unsigned char *&src, size_t length) {
        T *array = (T *) malloc(length * sizeof(T));
        memcpy(array, src, length * sizeof(T));
        src += length * sizeof(T);
        return array;
    }

}

#endif