#ifndef SZ3_SZN_HPP
#define SZ3_SZN_HPP

#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Config.hpp"
#include "szInterp.hpp"
#include "szLorenzoReg.hpp"
#include <cmath>


template<class T, uint N>
char *SZ_compress_N(SZ::Config &conf, T *data, size_t &outSize) {

    assert(N == conf.N);
    SZ::calAbsErrorBound(conf, data);

    char *cmpData;
    if (conf.cmprMethod == METHOD_LORENZO_REG) {
        cmpData = (char *) SZ_compress_LorenzoReg_N<T, N>(conf, data, outSize);
    } else if (conf.cmprMethod == METHOD_INTERP) {
        cmpData = (char *) SZ_compress_Interp_N<T, N>(conf, data, outSize);
    } else if (conf.cmprMethod == METHOD_INTERP_LORENZO) {
        cmpData = (char *) SZ_compress_Interp_lorenzo_N<T, N>(conf, data, outSize);
    }
    return cmpData;
}


template<class T, uint N>
T *SZ_decompress_N(SZ::Config &conf, char *cmpData, size_t cmpSize) {
    if (conf.cmprMethod == METHOD_LORENZO_REG) {
        return SZ_decompress_LorenzoReg_N<T, N>(conf, cmpData, cmpSize);
    } else if (conf.cmprMethod == METHOD_INTERP) {
        return SZ_decompress_Interp_N<T, N>(conf, cmpData, cmpSize);
    }
    return nullptr;
}
#endif