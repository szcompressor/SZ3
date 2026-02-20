#ifndef SZ3_IMPL_SZDISPATCHER_HPP
#define SZ3_IMPL_SZDISPATCHER_HPP

/**
 * @file SZDispatcher.hpp
 * @ingroup API
 * @brief Routes compression/decompression calls to the appropriate algorithm implementation.
 *
 * This is the central dispatch layer. Based on `Config::cmprAlgo`, it calls one of:
 * - `ALGO_LORENZO_REG`: `SZAlgoLorenzoReg.hpp` — blockwise prediction (Lorenzo/Regression) with `BlockwiseDecomposition`.
 * - `ALGO_INTERP` / `ALGO_INTERP_LORENZO`: `SZAlgoInterp.hpp` — interpolation-based decomposition.
 * - `ALGO_NOPRED`: `SZAlgoNopred.hpp` — quantization only (no predictor).
 * - `ALGO_BIOMD` / `ALGO_BIOMDXTC`: `SZAlgoBioMD.hpp` — molecular dynamics specific compression.
 * - `ALGO_SVD`: `SZAlgoSVD.hpp` — SVD-based decomposition.
 * - `ALGO_LOSSLESS`: Falls back to Zstd lossless-only compression.
 */

#include "SZ3/api/impl/SZAlgoInterp.hpp"
#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
#include "SZ3/api/impl/SZAlgoNopred.hpp"
#include "SZ3/api/impl/SZAlgoBioMD.hpp"
#include "SZ3/api/impl/SZAlgoSVD.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

/**
 * @brief Dispatch a compress call to the appropriate algorithm.
 *
 * Computes the absolute error bound from the config, then routes to the
 * correct algorithm. If the lossy compressed size provides a poor ratio (<3x),
 * automatically compares with Zstd lossless and picks the smaller result.
 *
 * @tparam T Data type
 * @tparam N Dimension
 * @param conf Configuration (may be mutated, e.g., `cmprAlgo` may be overridden to `ALGO_LOSSLESS`)
 * @param data Input data
 * @param cmpData Output compressed buffer
 * @param cmpCap Buffer capacity
 * @return Compressed size in bytes
 */
template <class T, uint N>
size_t SZ_compress_dispatcher(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    if (conf.cmprAlgo == ALGO_SVD && std::is_integral<T>::value) {
        throw std::invalid_argument("SVD algorithm only supports floating-point data types.");
    }
    assert(N == conf.N);
    calAbsErrorBound(conf, data);
    size_t cmpSize = 0;

    // if absErrorBound is 0, use lossless only mode
    if (conf.absErrorBound == 0) {
        conf.cmprAlgo = ALGO_LOSSLESS;
    }

    // do lossy compression
    bool isCmpCapSufficient = true;
    if (conf.cmprAlgo != ALGO_LOSSLESS) {
        try {
            std::vector<T> dataCopy(data, data + conf.num);
            if (conf.cmprAlgo == ALGO_LORENZO_REG) {
                cmpSize = SZ_compress_LorenzoReg<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
            } else if (conf.cmprAlgo == ALGO_INTERP) {
                cmpSize = SZ_compress_Interp<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
            } else if (conf.cmprAlgo == ALGO_INTERP_LORENZO) {
                cmpSize = SZ_compress_Interp_lorenzo<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
            } else if (conf.cmprAlgo == ALGO_NOPRED) {
                cmpSize = SZ_compress_nopred<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
            } else if (conf.cmprAlgo == ALGO_BIOMD) {
                return SZ_compress_bioMD<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
            } else if (conf.cmprAlgo == ALGO_BIOMDXTC) {
                return SZ_compress_bioMDXtcBased<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
            } else if (conf.cmprAlgo == ALGO_SVD) {
                cmpSize = SZ_compress_SVD<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
            } else {
                throw std::invalid_argument("Unknown compression algorithm");
            }

        } catch (std::length_error &e) {
            if (std::string(e.what()) == SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH) {
                isCmpCapSufficient = false;
                // printf("SZ is downgraded to lossless mode because the buffer for compressed data is not large enough.\n");
            } else {
                throw;
            }
        }
    }

    // do lossless only compression if 1) cmpr algorithm is lossless or 2) compressed buffer not large enough for lossy
    if (conf.cmprAlgo == ALGO_LOSSLESS || !isCmpCapSufficient) {
        conf.cmprAlgo = ALGO_LOSSLESS;
        auto zstd = Lossless_zstd();
        return zstd.compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), cmpData, cmpCap);
    }

    // if lossy compression ratio < 3, test if lossless only mode has a better ratio than lossy
    if (conf.num * sizeof(T) / 1.0 / cmpSize < 3) {
        auto zstd = Lossless_zstd();
        auto zstdCmpCap = ZSTD_compressBound(conf.num * sizeof(T)) + sizeof(size_t);
        auto zstdCmpData = static_cast<uchar *>(malloc(zstdCmpCap));
        size_t zstdCmpSize =
            zstd.compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), zstdCmpData, zstdCmpCap);
        if (zstdCmpSize < cmpSize && zstdCmpSize <= cmpCap) {
            conf.cmprAlgo = ALGO_LOSSLESS;
            memcpy(cmpData, zstdCmpData, zstdCmpSize);
            cmpSize = zstdCmpSize;
        }
        free(zstdCmpData);
    }
    return cmpSize;
}

/**
 * @brief Dispatch a decompress call to the appropriate algorithm.
 *
 * Routes based on `conf.cmprAlgo`, which was saved in the compressed data header.
 *
 * @tparam T Data type
 * @tparam N Dimension
 * @param conf Configuration (read from compressed data)
 * @param cmpData Compressed data buffer
 * @param cmpSize Compressed data size
 * @param decData Output decompressed data buffer
 */
template <class T, uint N>
void SZ_decompress_dispatcher(Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    if (conf.cmprAlgo == ALGO_LOSSLESS) {
        auto zstd = Lossless_zstd();
        size_t decDataSize = 0;
        auto decDataPos = reinterpret_cast<uchar *>(decData);
        zstd.decompress(cmpData, cmpSize, decDataPos, decDataSize);
        if (decDataSize != conf.num * sizeof(T)) {
            throw std::runtime_error("Decompressed data size does not match the original data size");
        }
    } else if (conf.cmprAlgo == ALGO_LORENZO_REG) {
        SZ_decompress_LorenzoReg<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == ALGO_INTERP) {
        SZ_decompress_Interp<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == ALGO_NOPRED) {
        SZ_decompress_nopred<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == ALGO_BIOMD) {
        SZ_decompress_bioMD<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == ALGO_BIOMDXTC) {
        SZ_decompress_bioMDXtcBased<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == ALGO_SVD) {
        SZ_decompress_SVD<T, N>(conf, cmpData, cmpSize, decData);
    } else {
        throw std::invalid_argument("Unknown compression algorithm");
    }
}
}  // namespace SZ3
#endif
