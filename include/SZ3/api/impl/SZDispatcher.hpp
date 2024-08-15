#ifndef SZ3_IMPL_SZDISPATCHER_HPP
#define SZ3_IMPL_SZDISPATCHER_HPP

#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/api/impl/SZAlgoInterp.hpp"
#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
#include "SZ3/api/impl/SZAlgo.hpp"
#include <cmath>

namespace SZ3 {
    template<class T, uint N>
    size_t SZ_compress_dispatcher(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
        
        assert(N == conf.N);
        calAbsErrorBound(conf, data);

//        char *cmpData;
        if (conf.absErrorBound == 0) {
            auto zstd = Lossless_zstd();
            return zstd.compress((uchar *) data, conf.num * sizeof(T), cmpData, cmpCap);
        } else if (conf.cmprAlgo == ALGO_LORENZO_REG) {
            return SZ_compress_LorenzoReg<T, N>(conf, data, cmpData, cmpCap);
        } else if (conf.cmprAlgo == ALGO_INTERP) {
            return SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
        } else if (conf.cmprAlgo == ALGO_INTERP_LORENZO) {
            return SZ_compress_Interp_lorenzo<T, N>(conf, data, cmpData, cmpCap);
        } else if (conf.cmprAlgo == ALGO_NOPRED) {
            return SZ_compress_nopred<T, N>(conf, data, cmpData, cmpCap);
        }
        return 0;
//        return cmpData;
    }
    
    template<class T, uint N>
    void SZ_decompress_dispatcher(Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
        if (conf.absErrorBound == 0) {
            auto zstd = Lossless_zstd();
            auto zstdDstCap = conf.num * sizeof(T);
            zstd.decompress(cmpData, cmpSize, (uchar *) decData, zstdDstCap);
        } else if (conf.cmprAlgo == ALGO_LORENZO_REG) {
            SZ_decompress_LorenzoReg<T, N>(conf, cmpData, cmpSize, decData);
        } else if (conf.cmprAlgo == ALGO_INTERP) {
            SZ_decompress_Interp<T, N>(conf, cmpData, cmpSize, decData);
        } else if (conf.cmprAlgo == ALGO_NOPRED) {
            SZ_decompress_nopred<T, N>(conf, cmpData, cmpSize, decData);
        } else {
            printf("SZ_decompress_dispatcher, Method not supported\n");
            exit(0);
        }
        
    }
}
#endif