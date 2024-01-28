#ifndef SZ_COMPRESSOR_HPP
#define SZ_COMPRESSOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"


namespace SZ3::concepts {

    /**
     * Compressor describes the whole compression workflow.
     * Each compressor represent a specific existing compression algorithm (interpolation, lorenzo + regression, etc.)
     * Compressor is implemented using modules (predictor, encoder, lossless, etc.).
     * @tparam T input data type
     */
    template<class T>
    class CompressorInterface {

    public:
        /**
         * decompress data to a new() space created in this function
         *
         * @param cmpData compressed data in bytes
         * @param cmpSize size of compressed data in bytes
         * @param num size of original data
         * @return pointer to the decompress data
         */
        virtual T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) = 0;

        /**
         * decompress data to a pre-allocated space
         * @param cmpData compressed data in bytes
         * @param cmpSize size of compressed data in bytes
         * @param decData pre-allocated space to store the decompress data
         * @return pointer to the decompress data
         */
        virtual T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) = 0;

        /**
         * compress data to a new() space
         * @param conf compression configuration
         * @param data input data in original format
         * @param compressed_size size of the compressed data in bytes
         * @return compressed data in bytes
         */
        virtual uchar *compress(const Config &conf, T *data, size_t &compressed_size) = 0;

    };
}

#endif
