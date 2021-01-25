#ifndef _SZ_COMPRESSOR_HPP
#define _SZ_COMPRESSOR_HPP

#include "def.hpp"

namespace SZ {
    template<class T>
    class Compressor {
    public:
        virtual T *decompress(uchar *compressed_data, size_t length, bool pre_de_lossless = false) = 0;

        virtual uchar *compress(T *data, size_t &compressed_size) = 0;
    };
}
#endif