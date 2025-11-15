#ifndef SZ3_IMPL_SZALGOSVD_HPP
#define SZ3_IMPL_SZALGOSVD_HPP

#include "SZ3/utils/Config.hpp"
#include "SZ3/decomposition/SVDDecomposition.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/compressor/SZGenericCompressor.hpp"

namespace SZ3 {

template <class T, uint N>
size_t SZ_compress_SVD(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    if constexpr (std::is_floating_point<T>::value) {
        LinearQuantizer<T> svd_quantizer(conf.absErrorBound * conf.svd_quant_eb_scale);
        LinearQuantizer<T> res_quantizer(conf.absErrorBound);
        auto compressor = make_compressor_sz_generic<T, N>(
            SVDDecomposition<T, N, LinearQuantizer<T>>(conf, svd_quantizer, res_quantizer),
            HuffmanEncoder<int>(),
            Lossless_zstd()
        );
        return compressor->compress(conf, data, cmpData, cmpCap);
    } else {
        throw std::invalid_argument("SVD algorithm only supports floating-point data types.");
    }
}

template <class T, uint N>
void SZ_decompress_SVD(Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    if constexpr (std::is_floating_point<T>::value) {
        LinearQuantizer<T> svd_quantizer(conf.absErrorBound * conf.svd_quant_eb_scale);
        LinearQuantizer<T> res_quantizer(conf.absErrorBound);
        auto compressor = make_compressor_sz_generic<T, N>(
            SVDDecomposition<T, N, LinearQuantizer<T>>(conf, svd_quantizer, res_quantizer),
            HuffmanEncoder<int>(),
            Lossless_zstd()
        );
        compressor->decompress(conf, cmpData, cmpSize, decData);
    } else {
        throw std::invalid_argument("SVD algorithm only supports floating-point data types.");
    }
}

} // namespace SZ3

#endif
