#ifndef SZ3_SZ_HPP
#define SZ3_SZ_HPP

#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZMetaFrontend.hpp"
#include "SZ3/frontend/SZGeneralFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/predictor/PolyRegressionPredictor.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Verification.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "szInterp.hpp"
#include "szLorenzoReg.hpp"
#include <cmath>
#include <memory>


template<class T, uint N>
char *SZ_compress_N(SZ::Config &conf, T *data, size_t &outSize) {

    assert(N == conf.N);
    SZ::calAbsErrorBound(conf, data);

    char *cmpData;
    if (conf.cmprMethod == METHOD_LORENZO_REG || conf.cmprMethod == METHOD_LORENZO_REG_FAST) {
        cmpData = (char *) SZ_compress_LorenzoReg_N<T, N>(conf, data, outSize);
    } else if (conf.cmprMethod == METHOD_INTERP) {
        cmpData = (char *) SZ_compress_Interp_N<T, N>(conf, data, outSize);
    } else if (conf.cmprMethod == METHOD_INTERP_LORENZO) {
        cmpData = (char *) SZ_compress_Interp_lorenzo_N<T, N>(conf, data, outSize);
    }
    memcpy(cmpData + outSize, &conf, sizeof(SZ::Config)); //Config.dims is a vector and only the pointer is saved.
    outSize += sizeof(SZ::Config); //TODO add save() and load() to Config
    return cmpData;
}


template<class T, uint N>
T *SZ_decompress_N(SZ::Config &conf, char *cmpData, size_t cmpSize) {
    auto dims = conf.dims;
    memcpy(&conf, cmpData + (cmpSize - sizeof(SZ::Config)), sizeof(SZ::Config));
    conf.dims = dims;
    conf.update_dims(dims.begin(), dims.end());

    if (conf.cmprMethod == METHOD_LORENZO_REG || conf.cmprMethod == METHOD_LORENZO_REG_FAST) {
        return SZ_decompress_LorenzoReg_N<T, N>(conf, cmpData, cmpSize);
    } else if (conf.cmprMethod == METHOD_INTERP) {
        return SZ_decompress_Interp_N<T, N>(conf, cmpData, cmpSize);
    }
    return nullptr;
}

template<class T>
char *SZ_compress(SZ::Config &conf, T *data, size_t &outSize) {
    if (conf.N == 1) {
        return SZ_compress_N<T, 1>(conf, data, outSize);
    } else if (conf.N == 2) {
        return SZ_compress_N<T, 2>(conf, data, outSize);
    } else if (conf.N == 3) {
        return SZ_compress_N<T, 3>(conf, data, outSize);
    } else if (conf.N == 4) {
        return SZ_compress_N<T, 4>(conf, data, outSize);
    } else {
        for (int i = 4; i < conf.N; i++) {
            conf.dims[3] *= conf.dims[i];
        }
        conf.dims.resize(4);
        conf.N = 4;
        return SZ_compress_N<T, 4>(conf, data, outSize);
    }
}

template<class T>
T *SZ_decompress(SZ::Config &conf, char *cmpData, size_t cmpSize) {
    if (conf.N == 1) {
        return SZ_decompress_N<T, 1>(conf, cmpData, cmpSize);
    } else if (conf.N == 2) {
        return SZ_decompress_N<T, 2>(conf, cmpData, cmpSize);
    } else if (conf.N == 3) {
        return SZ_decompress_N<T, 3>(conf, cmpData, cmpSize);
    } else if (conf.N == 4) {
        return SZ_decompress_N<T, 4>(conf, cmpData, cmpSize);
    } else {
        return SZ_decompress_N<T, 4>(conf, cmpData, cmpSize);
    }
}


#endif