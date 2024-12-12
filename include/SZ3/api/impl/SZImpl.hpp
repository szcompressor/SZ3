#ifndef SZ3_IMPL_SZ_HPP
#define SZ3_IMPL_SZ_HPP

#include "SZ3/api/impl/SZDispatcher.hpp"
#include "SZ3/api/impl/SZImplOMP.hpp"
#include "SZ3/def.hpp"

namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_impl(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
#ifndef _OPENMP
    conf.openmp = false;
#endif
    if (conf.openmp) {
        return SZ_compress_OMP<T, N>(conf, data, cmpData, cmpCap);
    } else {
        return SZ_compress_dispatcher<T, N>(conf, data, cmpData, cmpCap);
    }
}

template <class T, uint N>
void SZ_decompress_impl(Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
#ifndef _OPENMP
    conf.openmp = false;
#endif
    if (conf.openmp) {
        SZ_decompress_OMP<T, N>(conf, cmpData, cmpSize, decData);
    } else {
        SZ_decompress_dispatcher<T, N>(conf, cmpData, cmpSize, decData);
    }
}


template<class T>
size_t SZ_compress_size_bound(const Config &conf) {
    bool omp = conf.openmp;
#ifndef _OPENMP
    omp = false;
#endif
    if (omp) {
        return SZ_compress_size_bound_omp<T>(conf);
    } else {
        return conf.size_est() + ZSTD_compressBound(conf.num * sizeof(T));
    }
}

}  // namespace SZ3
#endif
