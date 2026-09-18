#ifndef SZ3_IMPL_SZDISPATCHER_OMP_HPP
#define SZ3_IMPL_SZDISPATCHER_OMP_HPP

#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>
#include <new>

#include "SZ3/api/impl/SZDispatcher.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"

#ifdef _OPENMP

#include <omp.h>

#include <stdexcept>

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
    if (conf.dims[0] < static_cast<size_t>(nThreads)) {
        nThreads = static_cast<int>(conf.dims[0]);
    }
    compressed_t.resize(nThreads);
    cmp_size_t.resize(nThreads + 1);
    cmp_start_t.resize(nThreads + 1);
    conf_t.resize(nThreads);
    min_t.resize(nThreads);
    max_t.resize(nThreads);
    // An exception that leaves an OpenMP region does not unwind to the caller -- the runtime calls
    // std::terminate -- so each thread stores its own and the first one is rethrown below, outside
    // the region, where an ordinary catch can see it.
    std::exception_ptr failure;
    // num_threads(nThreads) applies to this region alone. omp_set_num_threads(nThreads) would write
    // the process-wide nthreads-var instead, and that outlives the call: an application running 32
    // threads that compresses one chunk here would go on running 8 afterwards. SZ3 is a library
    // inside someone else's program and has no business changing that.
#pragma omp parallel num_threads(nThreads)
    try {
        int tid = omp_get_thread_num();

        auto dims_t = conf.dims;
        size_t lo = static_cast<size_t>(tid) * conf.dims[0] / nThreads;
        size_t hi = static_cast<size_t>(tid + 1) * conf.dims[0] / nThreads;
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
        // Room for the size header Lossless_zstd::compress writes ahead of the zstd stream.
        size_t cmp_size_cap = sizeof(size_t) + Lossless_zstd::compress_bound(conf_t[tid].num * sizeof(T));
        std::unique_ptr<uchar[]> compressed_owner(new uchar[cmp_size_cap]);
        compressed_t[tid] = compressed_owner.get();
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
    } catch (...) {
#pragma omp critical
        {
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) {
        std::rethrow_exception(failure);
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
    const uchar* const cmp_end = cmpData + cmpSize;
    int nThreads = 1;
    if (static_cast<size_t>(cmp_end - cmpr_data_pos) < sizeof(nThreads))
        throw std::out_of_range("SZ3 OMP: truncated thread count");
    read(nThreads, cmpr_data_pos);
    // Each thread contributes at least a config and a size, so the count cannot exceed the buffer size.
    if (nThreads <= 0 || static_cast<size_t>(nThreads) > cmpSize)
        throw std::out_of_range("SZ3 OMP: invalid thread count");

    std::vector<Config> conf_t(nThreads);
    for (int i = 0; i < nThreads; i++) {
        size_t confRemaining = static_cast<size_t>(cmp_end - cmpr_data_pos);
        conf_t[i].load(cmpr_data_pos, confRemaining);
    }

    if (conf_t[0].sz3MagicNumber != SZ3_MAGIC_NUMBER) {
        throw std::invalid_argument("magic number mismatch, the input data is not compressed by SZ3");
    }
    if (versionStr(conf_t[0].sz3DataVer) != SZ3_DATA_VER) {
        // Same contract as the serial path in sz.hpp: the exception carries the whole diagnostic,
        // so nothing here writes to the caller's stdout.
        std::stringstream ss;
        ss << "SZ3 " << SZ3_VER << " reads data version " << SZ3_DATA_VER << ", but this data is version "
           << versionStr(conf_t[0].sz3DataVer) << ". Use SZ3 v" << versionStr(conf_t[0].sz3DataVer)
           << " to decompress it.";
        throw std::invalid_argument(ss.str());
    }

    std::vector<size_t> cmp_start_t, cmp_size_t;
    cmp_size_t.resize(nThreads);
    if (static_cast<size_t>(cmp_end - cmpr_data_pos) < static_cast<size_t>(nThreads) * sizeof(size_t))
        throw std::out_of_range("SZ3 OMP: truncated per-thread sizes");
    read(cmp_size_t.data(), nThreads, cmpr_data_pos);
    auto cmpr_data_p = cmpr_data_pos;

    cmp_start_t.resize(nThreads + 1);
    cmp_start_t[0] = 0;
    // The payloads follow back-to-back; build the offsets without overflowing so no slice points out.
    const size_t payload_avail = static_cast<size_t>(cmp_end - cmpr_data_p);
    for (int i = 1; i <= nThreads; i++) {
        if (cmp_size_t[i - 1] > payload_avail - cmp_start_t[i - 1])
            throw std::out_of_range("SZ3 OMP: per-thread compressed sizes exceed the buffer");
        cmp_start_t[i] = cmp_start_t[i - 1] + cmp_size_t[i - 1];
    }

    std::exception_ptr failure;
#pragma omp parallel num_threads(nThreads)
    try {
        // nThreads is how many chunks the writer split the data into -- it came out of the stream,
        // not from this machine. num_threads asks for that many threads but nothing guarantees
        // them: OMP_THREAD_LIMIT caps the whole program, dynamic adjustment may lower it, and a
        // region nested inside the caller's own gets exactly one thread unless nested parallelism
        // was turned on, which it is not by default -- the case a filter called from inside an
        // application's parallel region lands in. Indexing the chunks by thread id would then skip
        // every chunk whose id no thread has, leaving that span of the output untouched and
        // reporting success, so each thread walks the chunks in strides of however many arrived.
        const int actual = omp_get_num_threads();
        for (int tid = omp_get_thread_num(); tid < nThreads; tid += actual) {
            auto dims_t = conf.dims;
            size_t lo = static_cast<size_t>(tid) * conf.dims[0] / nThreads;
            size_t hi = static_cast<size_t>(tid + 1) * conf.dims[0] / nThreads;
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
    } catch (...) {
#pragma omp critical
        {
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) {
        std::rethrow_exception(failure);
    }
#else
    throw std::invalid_argument(
        "SZ3: this data was compressed with OpenMP; decompressing it needs an OpenMP-enabled build");
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
    if (conf.dims[0] < static_cast<size_t>(nThreads)) {
        nThreads = static_cast<int>(conf.dims[0]);
    }
    size_t chunk_size = conf.dims[0] / static_cast<size_t>(nThreads) * (conf.num / conf.dims[0]);
    size_t last_chunk_size = (conf.dims[0] - conf.dims[0] / nThreads * (nThreads - 1)) * (conf.num / conf.dims[0]);
    // for each thread, we save conf, compressed size, and compressed data
    // the per-chunk compressed data may carry the size header written by Lossless_zstd::compress
    return sizeof(int) + nThreads * conf.size_est() + 2 * nThreads * sizeof(size_t) +
           (nThreads - 1) * Lossless_zstd::compress_bound(chunk_size * sizeof(T)) +
           Lossless_zstd::compress_bound(last_chunk_size * sizeof(T));
#else
    return conf.size_est() + Lossless_zstd::compress_bound(conf.num * sizeof(T));
#endif
}
}  // namespace SZ3

#endif