#ifndef SZ3_COMPRESSOR_ZFP_HPP
#define SZ3_COMPRESSOR_ZFP_HPP

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/zfp/zfpcodec1.h"
#include "SZ3/utils/zfp/zfpcodec2.h"
#include "SZ3/utils/zfp/zfpcodec3.h"

namespace SZ3 {

template <class T, uint N>
class ZFPCompressor : public concepts::CompressorInterface<T> {
   public:
    size_t compress(const Config &conf, T *data, uchar *cmpData, size_t cmpCap) override {
        ZFP::MemoryBitStream stream;
        stream.open(cmpData, cmpCap);
        int expmin = INT_MIN;
        if (conf.absErrorBound > 0) {
            frexp(conf.absErrorBound, &expmin);
            expmin--;
        }
        const T *p = data;

        if constexpr (N == 1) {
            ZFP::Codec1<ZFP::MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
            for (auto x = 0; x < conf.dims[0]; x += 4, p += 4) {
                codec.encode(p, 1, codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u)));
            }
        } else if constexpr (N == 2) {
            ZFP::Codec2<ZFP::MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
            for (uint y = 0; y < conf.dims[1]; y += 4, p += 4 * (conf.dims[0] - (conf.dims[0] + 3) / 4)) {
                for (uint x = 0; x < conf.dims[0]; x += 4, p += 4) {
                    codec.encode(p, 1, conf.dims[0],
                                 codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u),
                                            std::min(static_cast<uint>(conf.dims[1] - y), 4u)));
                }
            }
        } else if constexpr (N == 3) {
            ZFP::Codec3<ZFP::MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
            for (uint z = 0; z < conf.dims[2];
                 z += 4, p += 4 * conf.dims[0] * (conf.dims[1] - (conf.dims[1] + 3) / 4)) {
                for (uint y = 0; y < conf.dims[1]; y += 4, p += 4 * (conf.dims[0] - (conf.dims[0] + 3) / 4)) {
                    for (uint x = 0; x < conf.dims[0]; x += 4, p += 4) {
                        codec.encode(p, 1, conf.dims[0], conf.dims[0] * conf.dims[1],
                                     codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u),
                                                std::min(static_cast<uint>(conf.dims[1] - y), 4u),
                                                std::min(static_cast<uint>(conf.dims[2] - z), 4u)));
                    }
                }
            }
        }
        stream.flush();
        return stream.size();
    }

    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) override {
        int expmin = INT_MIN;
        if (conf.absErrorBound > 0) {
            frexp(conf.absErrorBound, &expmin);
            expmin--;
        }
        ZFP::MemoryBitStream stream;
        stream.open(const_cast<uchar *>(cmpData), cmpSize);
        float *p = decData;

        if constexpr (N == 1) {
            ZFP::Codec1<ZFP::MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
            for (auto x = 0; x < conf.dims[0]; x += 4, p += 4) {
                codec.decode(p, 1, codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u)));
            }
        } else if constexpr (N == 2) {
            ZFP::Codec2<ZFP::MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
            for (uint y = 0; y < conf.dims[1]; y += 4, p += 4 * (conf.dims[0] - (conf.dims[0] + 3) / 4)) {
                for (uint x = 0; x < conf.dims[0]; x += 4, p += 4) {
                    codec.decode(p, 1, conf.dims[0],
                                 codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u),
                                            std::min(static_cast<uint>(conf.dims[1] - y), 4u)));
                }
            }
        } else if constexpr (N == 3) {
            ZFP::Codec3<ZFP::MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
            for (uint z = 0; z < conf.dims[2];
                 z += 4, p += 4 * conf.dims[0] * (conf.dims[1] - (conf.dims[1] + 3) / 4)) {
                for (uint y = 0; y < conf.dims[1]; y += 4, p += 4 * (conf.dims[0] - (conf.dims[0] + 3) / 4)) {
                    for (uint x = 0; x < conf.dims[0]; x += 4, p += 4) {
                        codec.decode(p, 1, conf.dims[0], conf.dims[0] * conf.dims[1],
                                     codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u),
                                                std::min(static_cast<uint>(conf.dims[1] - y), 4u),
                                                std::min(static_cast<uint>(conf.dims[2] - z), 4u)));
                    }
                }
            }
        }
        return decData;
    }
};

template <class T, uint N>
std::shared_ptr<ZFPCompressor<T, N>> make_compressor_zfp() {
    return std::make_shared<ZFPCompressor<T, N>>();
}

// template <class T, uint N>
// void zfp_compress_test(Config &conf, T *data, char *dst, size_t &outSize) {
//     typedef typename ZFP::ScalarTraits<T>::Fixed Fixed;
//     typedef typename ZFP::ScalarTraits<T>::Int Int;
//     typedef typename ZFP::ScalarTraits<T>::UInt UInt;
//
//     static const uint ebits = ZFP::ScalarTraits<T>::ebits;  // number of exponent bits
//     static const int ebias = (1 << (ebits - 1)) - 1;        // floating-point exponent bias
//
//     calAbsErrorBound(conf, data);
//
//     ZFP::MemoryBitStream stream;
//     stream.open(dst, outSize);
//     const T *p = data;
//
//     int expmin = INT_MIN;
//     if (conf.absErrorBound > 0) {
//         frexp(conf.absErrorBound, &expmin);
//         expmin--;
//     }
//     uint maxprec = CHAR_BIT * sizeof(Int);
//     int minexp = std::max(expmin, std::numeric_limits<T>::min_exponent - std::numeric_limits<T>::digits);
//
//     std::vector<uint> dims(conf.dims.begin(), conf.dims.end());
//     std::vector<uint> s(N, 1);  // stride
//     for (uint i = 1; i < N; i++) {
//         s[i] = s[i - 1] * dims[i - 1];
//     }
//     if constexpr (N == 1) {
//         for (auto x = 0; x < dims[0]; x += 4, p += 4) {
//             int block_size = 4;
//             Fixed q[block_size];
//             int emax = 0;
//             uint nx = std::min(dims[0] - x, 4u);
//             if (nx == 4) {
//                 // convert to fixed-point
//                 emax = ZFP::Codec1<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0]);
//                 // perform block transform
//                 ZFP::Codec1<ZFP::MemoryBitStream, T>::fwd_xform(q);
//             } else {
//                 emax = ZFP::Codec1<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], nx);
//                 ZFP::Codec1<ZFP::MemoryBitStream, T>::fwd_xform(q, nx);
//             }
//             // reorder and convert to integer
//             Int buffer[block_size];
//             for (uint i = 0; i < block_size; i++) {
//                 buffer[i] = q[i].reinterpret();
//             }
//             stream.write(emax + ebias, ebits);
//             uint precision = std::min(maxprec, static_cast<uint>(std::max(0, emax - minexp + 4)));
//             ZFP::IntCodec04<ZFP::MemoryBitStream, Int, UInt>::encode(stream, buffer, 0, UINT_MAX, precision);
//         }
//     } else if constexpr (N == 2) {
//         for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
//             for (uint x = 0; x < dims[0]; x += 4, p += 4) {
//                 int block_size = 16;
//                 Fixed q[block_size];
//                 int emax = 0;
//                 uint nx = std::min(dims[0] - x, 4u);
//                 uint ny = std::min(dims[1] - y, 4u);
//                 if (nx == 4 && ny == 4) {
//                     // convert to fixed-point
//                     emax = ZFP::Codec2<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1]);
//                     // perform block transform
//                     ZFP::Codec2<ZFP::MemoryBitStream, T>::fwd_xform(q);
//                 } else {
//                     emax = ZFP::Codec2<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], nx, ny);
//                     ZFP::Codec2<ZFP::MemoryBitStream, T>::fwd_xform(q, nx, ny);
//                 }
//                 // reorder and convert to integer
//                 Int buffer[block_size];
//                 for (uint i = 0; i < block_size; i++)
//                     buffer[i] = q[ZFP::Codec2<ZFP::MemoryBitStream, T>::perm[i]].reinterpret();
//                 stream.write(emax + ebias, ebits);
//                 auto precision = std::min(maxprec, static_cast<uint>(std::max(0, emax - minexp + 6)));
//                 ZFP::IntCodec16<ZFP::MemoryBitStream, Int, UInt>::encode(stream, buffer, 0, UINT_MAX, precision);
//             }
//         }
//     } else if constexpr (N == 3) {
//         for (uint z = 0; z < dims[2]; z += 4, p += 4 * dims[0] * (dims[1] - (dims[1] + 3) / 4)) {
//             for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
//                 for (uint x = 0; x < dims[0]; x += 4, p += 4) {
//                     int block_size = 64;
//                     Fixed q[block_size];
//                     int emax = 0;
//                     uint nx = std::min(dims[0] - x, 4u);
//                     uint ny = std::min(dims[1] - y, 4u);
//                     uint nz = std::min(dims[2] - z, 4u);
//                     if (nx == 4 && ny == 4 && nz == 4) {
//                         // convert to fixed-point
//                         emax = ZFP::Codec3<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], s[2]);
//                         // perform block transform
//                         ZFP::Codec3<ZFP::MemoryBitStream, T>::fwd_xform(q);
//                     } else {
//                         emax = ZFP::Codec3<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], s[2], nx, ny, nz);
//                         ZFP::Codec3<ZFP::MemoryBitStream, T>::fwd_xform(q, nx, ny, nz);
//                     }
//                     // reorder and convert to integer
//                     Int buffer[block_size];
//                     for (uint i = 0; i < block_size; i++)
//                         buffer[i] = q[ZFP::Codec3<ZFP::MemoryBitStream, T>::perm[i]].reinterpret();
//                     stream.write(emax + ebias, ebits);
//                     auto precision = std::min(maxprec, static_cast<uint>(std::max(0, emax - minexp + 8)));
//                     ZFP::IntCodec64<ZFP::MemoryBitStream, Int, UInt>::encode(stream, buffer, 0, UINT_MAX, precision);
//                 }
//             }
//         }
//     }
//     stream.flush();
//     outSize = stream.size();
//     printf("outSize=%zu\n", outSize);
// }

}  // namespace SZ3
#endif
