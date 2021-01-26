#ifndef SZ_COMPRESSOR_HPP
#define SZ_COMPRESSOR_HPP

#include "def.hpp"

namespace SZ {
    namespace concepts {
        template<class T>
        class CompressorInterface {
        public:
            virtual T *decompress(uchar const *compressed_data, size_t length) = 0;

            virtual uchar *compress(T *data, size_t &compressed_size) = 0;
        };
    }
}
#endif