#ifndef SZ3_IMPL_SZDISPATCHER_HPP
#define SZ3_IMPL_SZDISPATCHER_HPP

#include "SZ3/api/impl/SZAlgoInterp.hpp"
#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
#include "SZ3/api/impl/SZAlgoNopred.hpp"
#include "SZ3/api/impl/SZAlgoBioMD.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_dispatcher(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    calAbsErrorBound(conf, data);
    size_t cmpSize = 0;
    bool isCmpCapSufficient = true;

    if (conf.absErrorBound != 0) {
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
            } else {
                printf("SZ_compress_dispatcher, Method not supported\n");
                exit(0);
            }

        } catch (std::length_error &e) {
            if (std::string(e.what()) == SZ_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH) {
                isCmpCapSufficient = false;
                printf(
                    "The buffer for compressed data is not large enough. Ideally, set it as 2X original data size.\n "
                    "SZ is downgraded to lossless mode.\n");
            } else {
                throw;
            }
        }
    }
    if (conf.absErrorBound == 0 || !isCmpCapSufficient) {
        // must use lossless only mode if 1) user sets error bound to 0, or 2) compressed buffer not large enough for
        // lossy
        conf.cmprAlgo = ALGO_LOSSLESS;
        auto zstd = Lossless_zstd();
        cmpSize = zstd.compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), cmpData, cmpCap);
    } else {
        if (conf.num * sizeof(T) / 1.0 / cmpSize < 3) {
            // may use lossless only mode when lossy compression ratio is less than 2.5
            auto zstd = Lossless_zstd();
            auto zstdCmpCap = ZSTD_compressBound(conf.num * sizeof(T));
            auto zstdCmpData = static_cast<uchar *>(malloc(cmpCap));
            size_t zstdCmpSize =
                zstd.compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), zstdCmpData, zstdCmpCap);
            if (zstdCmpSize < cmpSize) {
                conf.cmprAlgo = ALGO_LOSSLESS;
                if (zstdCmpSize > cmpCap) {
                    throw std::length_error(SZ_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
                }
                memcpy(cmpData, zstdCmpData, zstdCmpSize);
                cmpSize = zstdCmpSize;
            }
            free(zstdCmpData);
        }
    }
    return cmpSize;
}

template <class T, uint N>
void SZ_decompress_dispatcher(Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    if (conf.cmprAlgo == ALGO_LOSSLESS) {
        auto zstd = Lossless_zstd();
        size_t decDataSize = 0;
        auto decDataPos = reinterpret_cast<uchar *>(decData);
        zstd.decompress(cmpData, cmpSize, decDataPos, decDataSize);
        if (decDataSize != conf.num * sizeof(T)) {
            throw std::runtime_error("Decompressed data size does not match the original data size\n");
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
    } else {
        printf("SZ_decompress_dispatcher, Method not supported\n");
        exit(0);
    }
}
}  // namespace SZ3
#endif
