#ifndef SZ3_IMPL_SZDISPATCHER_OMP_HPP
#define SZ3_IMPL_SZDISPATCHER_OMP_HPP

#include <cmath>
#include <memory>

#include "SZ3/api/impl/SZDispatcher.hpp"

#ifdef _OPENMP

#include <omp.h>

#endif
namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_OMP(Config& conf, const T* data, uchar* cmpData, size_t cmpCap) {
#ifdef _OPENMP
    unsigned char* buffer_pos = cmpData;

    std::vector<uchar*> compressed_t;
    std::vector<size_t> cmp_size_t, cmp_start_t;
    std::vector<T> min_t, max_t;
    std::vector<Config> conf_t;
    //    Timer timer(true);
    int nThreads = 1;
    // double eb;
#pragma omp parallel
#pragma omp single
    {
        nThreads = omp_get_num_threads();
    }
    if (conf.dims[0] < nThreads) {
        nThreads = conf.dims[0];
        omp_set_num_threads(nThreads);
    }
    printf("OpenMP enabled for compression, threads = %d\n", nThreads);
    compressed_t.resize(nThreads);
    cmp_size_t.resize(nThreads + 1);
    cmp_start_t.resize(nThreads + 1);
    conf_t.resize(nThreads);
    min_t.resize(nThreads);
    max_t.resize(nThreads);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        auto dims_t = conf.dims;
        int lo = tid * conf.dims[0] / nThreads;
        int hi = (tid + 1) * conf.dims[0] / nThreads;
        dims_t[0] = hi - lo;
        auto it = dims_t.begin();
        size_t num_t_base = std::accumulate(++it, dims_t.end(), static_cast<size_t>(1), std::multiplies<size_t>());
        size_t num_t = dims_t[0] * num_t_base;

        const T* data_t = data + lo * num_t_base;
        // std::vector<T> data_t(data + lo * num_t_base, data + lo * num_t_base + num_t);
        if (conf.errorBoundMode != EB_ABS) {
            auto minmax = std::minmax_element(data_t, data_t + num_t);
            min_t[tid] = *minmax.first;
            max_t[tid] = *minmax.second;
#pragma omp barrier
#pragma omp single
            {
                T range = *std::max_element(max_t.begin(), max_t.end()) - *std::min_element(min_t.begin(), min_t.end());
                calAbsErrorBound<T>(conf, data, range);
                //                timer.stop("OMP init");
                //                timer.start();
            }
        }

        conf_t[tid] = conf;
        conf_t[tid].setDims(dims_t.begin(), dims_t.end());
        size_t cmp_size_cap = ZSTD_compressBound(conf_t[tid].num * sizeof(T));
        compressed_t[tid] = static_cast<uchar*>(malloc(cmp_size_cap));
        // we have to use conf_t[tid].N instead of N since each chunk may be a slice of the original data
        if (conf_t[tid].N == 1) {
            cmp_size_t[tid] = SZ_compress_dispatcher<T, 1>(conf_t[tid], data_t, compressed_t[tid], cmp_size_cap);
        } else if (conf_t[tid].N == 2) {
            cmp_size_t[tid] = SZ_compress_dispatcher<T, 2>(conf_t[tid], data_t, compressed_t[tid], cmp_size_cap);
        } else if (conf_t[tid].N == 3) {
            cmp_size_t[tid] = SZ_compress_dispatcher<T, 3>(conf_t[tid], data_t, compressed_t[tid], cmp_size_cap);
        } else if (conf_t[tid].N == 4) {
            cmp_size_t[tid] = SZ_compress_dispatcher<T, 4>(conf_t[tid], data_t, compressed_t[tid], cmp_size_cap);
        } else {
            throw std::invalid_argument("Unsupported N");
        }

#pragma omp barrier
#pragma omp single
        {
            //            timer.stop("OMP compression");
            //            timer.start();
            cmp_start_t[0] = 0;
            for (int i = 1; i <= nThreads; i++) {
                cmp_start_t[i] = cmp_start_t[i - 1] + cmp_size_t[i - 1];
            }
            // size_t bufferSize = sizeof(int) + (nThreads + 1) * Config::size_est() + cmp_start_t[nThreads];
            //                buffer = new uchar[bufferSize];
            //                buffer_pos = buffer;
            write(nThreads, buffer_pos);
            for (int i = 0; i < nThreads; i++) {
                conf_t[i].save(buffer_pos);
            }
            write(cmp_size_t.data(), nThreads, buffer_pos);
        }

        memcpy(buffer_pos + cmp_start_t[tid], compressed_t[tid], cmp_size_t[tid]);
        free(compressed_t[tid]);
    }

    return buffer_pos - cmpData + cmp_start_t[nThreads];
    //    timer.stop("OMP memcpy");

#else
    return SZ_compress_dispatcher<T, N>(conf, data, cmpData, cmpCap);
#endif
}

template <class T, uint N>
void SZ_decompress_OMP(Config& conf, const uchar* cmpData, size_t cmpSize, T* decData) {
#ifdef _OPENMP

    auto cmpr_data_pos = cmpData;
    int nThreads = 1;
    read(nThreads, cmpr_data_pos);
    omp_set_num_threads(nThreads);
    printf("OpenMP enabled for decompression, threads = %d\n", nThreads);

    std::vector<Config> conf_t(nThreads);
    for (int i = 0; i < nThreads; i++) {
        conf_t[i].load(cmpr_data_pos);
    }

    std::vector<size_t> cmp_start_t, cmp_size_t;
    cmp_size_t.resize(nThreads);
    read(cmp_size_t.data(), nThreads, cmpr_data_pos);
    auto cmpr_data_p = cmpr_data_pos;

    cmp_start_t.resize(nThreads + 1);
    cmp_start_t[0] = 0;
    for (int i = 1; i <= nThreads; i++) {
        cmp_start_t[i] = cmp_start_t[i - 1] + cmp_size_t[i - 1];
    }

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto dims_t = conf.dims;
        int lo = tid * conf.dims[0] / nThreads;
        int hi = (tid + 1) * conf.dims[0] / nThreads;
        dims_t[0] = hi - lo;
        auto it = dims_t.begin();
        size_t num_t_base = std::accumulate(++it, dims_t.end(), static_cast<size_t>(1), std::multiplies<size_t>());

        if (conf_t[tid].N == 1) {
            SZ_decompress_dispatcher<T, 1>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid],
                                           decData + lo * num_t_base);
        } else if (conf_t[tid].N == 2) {
            SZ_decompress_dispatcher<T, 2>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid],
                                           decData + lo * num_t_base);
        } else if (conf_t[tid].N == 3) {
            SZ_decompress_dispatcher<T, 3>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid],
                               decData + lo * num_t_base);
        } else if (conf_t[tid].N == 4) {
            SZ_decompress_dispatcher<T, 4>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid],
                               decData + lo * num_t_base);
        } else {
            throw std::invalid_argument("Unsupported N");
        }
    }
#else
    SZ_decompress_dispatcher<T, N>(conf, cmpData, cmpSize, decData);
#endif
}

template <class T>
size_t SZ_compress_size_bound_omp(const Config& conf) {
#ifdef _OPENMP
    int nThreads = 1;
#pragma omp parallel
#pragma omp single
    {
        nThreads = omp_get_num_threads();
    }
    if (conf.dims[0] < nThreads) {
        nThreads = conf.dims[0];
    }
    size_t chunk_size = conf.dims[0] / nThreads * (conf.num / conf.dims[0]);
    size_t last_chunk_size = (conf.dims[0] - conf.dims[0] / nThreads * (nThreads - 1)) * (conf.num / conf.dims[0]);
    //for each thread, we save conf, compressed size, and compressed data
    return sizeof(int) + nThreads * conf.size_est() + nThreads * sizeof(size_t) +
           (nThreads - 1) * ZSTD_compressBound(chunk_size * sizeof(T)) +
           ZSTD_compressBound(last_chunk_size * sizeof(T));
#else
    return conf.size_est() + ZSTD_compressBound(conf.num * sizeof(T));
#endif
}
} // namespace SZ3

#endif