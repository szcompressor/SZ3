#ifndef SZ3_SZ_BIOMD_HPP
#define SZ3_SZ_BIOMD_HPP

#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/decomposition/SZBioMDXtcDecomposition.hpp"
#include "SZ3/decomposition/SZBioMDDecomposition.hpp"
#include "SZ3/encoder/XtcBasedEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/def.hpp"

namespace SZ3 {
    
    template<class T, uint N>
    size_t SZ_compress_bioMD(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
        assert(N == conf.N);
        assert(conf.cmprAlgo == ALGO_BIOMD);
        calAbsErrorBound(conf, data);
        
        auto quantizer = LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
        auto sz = make_compressor_sz_generic<T, N>(make_decomposition_biomd<T, N>(conf, quantizer), HuffmanEncoder<int>(),
                                                   Lossless_zstd());
        return sz->compress(conf, data, cmpData, cmpCap);
    }
    
    template<class T, uint N>
    void SZ_decompress_bioMD(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
        assert(conf.cmprAlgo == ALGO_BIOMD);
        
        LinearQuantizer<T> quantizer;
        auto sz = make_compressor_sz_generic<T, N>(make_decomposition_biomd<T, N>(conf, quantizer),
                                                   HuffmanEncoder<int>(), Lossless_zstd());
        sz->decompress(conf, cmpData, cmpSize, decData);
    }
    
    template<class T, uint N>
    size_t SZ_compress_bioMDXtcBased(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
        assert(N == conf.N);
        assert(conf.cmprAlgo == ALGO_BIOMDXTC);
        calAbsErrorBound(conf, data);
        
        auto sz = make_compressor_sz_generic<T, N>(SZBioMDXtcDecomposition<T, N>(conf), XtcBasedEncoder<int>(),
                                                   Lossless_bypass());
        return sz->compress(conf, data, cmpData, cmpCap);
    }
    
    template<class T, uint N>
    void SZ_decompress_bioMDXtcBased(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
        assert(conf.cmprAlgo == ALGO_BIOMDXTC);
        
        auto sz = make_compressor_sz_generic<T, N>(SZBioMDXtcDecomposition<T, N>(conf),
                                                   XtcBasedEncoder<int>(), Lossless_bypass());
        sz->decompress(conf, cmpData, cmpSize, decData);
    }
    
}
#endif
