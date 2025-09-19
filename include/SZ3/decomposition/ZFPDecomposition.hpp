#ifndef SZ3_ZFP_DECOMPOSITION_HPP
#define SZ3_ZFP_DECOMPOSITION_HPP

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/zfp/zfpcodec.h"
#include "SZ3/utils/zfp/zfpcodec1.h"
#include "SZ3/utils/zfp/zfpcodec2.h"
#include "SZ3/utils/zfp/zfpcodec3.h"

namespace SZ3 {
template <class T, class To, uint N>
class ZFPDecomposition : public concepts::DecompositionInterface<T, To, N> {
   public:
    typedef typename ZFP::ScalarTraits<T>::Fixed Fixed;

    std::vector<To> compress(const Config &conf, T *data) override {
        const T *p = data;
        size_t n_blocks = 1;
        int block_size = 1;
        for (uint i = 0; i < N; i++) {
            n_blocks *= (conf.dims[i] + 3) / 4;
            block_size *= 4;
        }
        std::vector<int> output(1 + n_blocks + n_blocks * block_size);
        output[0] = n_blocks;
        auto output_emax = &output[1];
        auto output_data = &output[n_blocks + 1];

        std::vector<uint> dims(conf.dims.begin(), conf.dims.end());
        std::vector<uint> s(N, 1);  // stride
        for (uint i = 1; i < N; i++) {
            s[i] = s[i - 1] * dims[i - 1];
        }

        if constexpr (N == 1) {
            for (auto x = 0; x < dims[0]; x += 4, p += 4) {
                Fixed q[block_size];
                int emax = 0;
                uint nx = std::min(dims[0] - x, 4u);
                if (nx == 4) {
                    // convert to fixed-point
                    emax = ZFP::Codec1<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0]);
                    // perform block transform
                    ZFP::Codec1<ZFP::MemoryBitStream, T>::fwd_xform(q);
                } else {
                    emax = ZFP::Codec1<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], nx);
                    ZFP::Codec1<ZFP::MemoryBitStream, T>::fwd_xform(q, nx);
                }
                // reorder and convert to integer
                *output_emax++ = emax;
                for (uint i = 0; i < block_size; i++) {
                    *output_data++ = q[i].reinterpret();
                }
            }
        } else if constexpr (N == 2) {
            for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
                for (uint x = 0; x < dims[0]; x += 4, p += 4) {
                    Fixed q[block_size];
                    int emax = 0;
                    uint nx = std::min(dims[0] - x, 4u);
                    uint ny = std::min(dims[1] - y, 4u);
                    if (nx == 4 && ny == 4) {
                        // convert to fixed-point
                        emax = ZFP::Codec2<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1]);
                        // perform block transform
                        ZFP::Codec2<ZFP::MemoryBitStream, T>::fwd_xform(q);
                    } else {
                        emax = ZFP::Codec2<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], nx, ny);
                        ZFP::Codec2<ZFP::MemoryBitStream, T>::fwd_xform(q, nx, ny);
                    }
                    // reorder and convert to integer
                    *output_emax++ = emax;
                    for (uint i = 0; i < block_size; i++) {
                        *output_data++ = q[ZFP::Codec2<ZFP::MemoryBitStream, T>::perm[i]].reinterpret();
                    }
                }
            }
        } else if constexpr (N == 3) {
            for (uint z = 0; z < dims[2]; z += 4, p += 4 * dims[0] * (dims[1] - (dims[1] + 3) / 4)) {
                for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
                    for (uint x = 0; x < dims[0]; x += 4, p += 4) {
                        Fixed q[block_size];
                        int emax = 0;
                        uint nx = std::min(dims[0] - x, 4u);
                        uint ny = std::min(dims[1] - y, 4u);
                        uint nz = std::min(dims[2] - z, 4u);
                        if (nx == 4 && ny == 4 && nz == 4) {
                            // convert to fixed-point
                            emax = ZFP::Codec3<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], s[2]);
                            // perform block transform
                            ZFP::Codec3<ZFP::MemoryBitStream, T>::fwd_xform(q);
                        } else {
                            emax = ZFP::Codec3<ZFP::MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], s[2], nx, ny, nz);
                            ZFP::Codec3<ZFP::MemoryBitStream, T>::fwd_xform(q, nx, ny, nz);
                        }
                        *output_emax++ = emax;
                        // reorder and convert to integer
                        for (uint i = 0; i < block_size; i++) {
                            *output_data++ = q[ZFP::Codec3<ZFP::MemoryBitStream, T>::perm[i]].reinterpret();
                        }
                    }
                }
            }
        }
        return output;
    }
    T *decompress(const Config &conf, std::vector<To> &transformed, T *dec_data) override {
        T *p = dec_data;
        size_t n_blocks = transformed[0];
        auto emax_pos = &transformed[1];
        auto transformed_pos = &transformed[1 + n_blocks];

        std::vector<uint> dims(conf.dims.begin(), conf.dims.end());
        std::vector<uint> s(N, 1);  // stride
        for (uint i = 1; i < N; i++) {
            s[i] = s[i - 1] * dims[i - 1];
        }
        if constexpr (N == 1) {
            for (auto x = 0; x < dims[0]; x += 4, p += 4) {
                int block_size = 4;
                int emax = *emax_pos++;
                Fixed q[block_size];
                for (uint i = 0; i < block_size; i++) {
                    q[i] = Fixed::reinterpret(*transformed_pos++);
                }
                ZFP::Codec1<ZFP::MemoryBitStream, T>::inv_xform(q);
                uint nx = std::min(dims[0] - x, 4u);
                if (nx == 4) {
                    ZFP::Codec1<ZFP::MemoryBitStream, T>::inv_cast(q, p, s[0], emax);
                } else {
                    ZFP::Codec1<ZFP::MemoryBitStream, T>::inv_cast(q, p, s[0], nx, emax);
                }
            }
        } else if constexpr (N == 2) {
            for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
                for (uint x = 0; x < dims[0]; x += 4, p += 4) {
                    int block_size = 16;
                    int emax = *emax_pos++;
                    Fixed q[block_size];
                    for (uint i = 0; i < block_size; i++) {
                        q[ZFP::Codec2<ZFP::MemoryBitStream, T>::perm[i]] = Fixed::reinterpret(*transformed_pos++);
                    }
                    ZFP::Codec2<ZFP::MemoryBitStream, T>::inv_xform(q);
                    uint nx = std::min(dims[0] - x, 4u);
                    uint ny = std::min(dims[1] - y, 4u);
                    if (nx == 4 && ny == 4) {
                        ZFP::Codec2<ZFP::MemoryBitStream, T>::inv_cast(q, p, s[0], s[1], emax);
                    } else {
                        ZFP::Codec2<ZFP::MemoryBitStream, T>::inv_cast(q, p, s[0], s[1], nx, ny, emax);
                    }
                }
            }
        } else if constexpr (N == 3) {
            for (uint z = 0; z < dims[2]; z += 4, p += 4 * dims[0] * (dims[1] - (dims[1] + 3) / 4)) {
                for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
                    for (uint x = 0; x < dims[0]; x += 4, p += 4) {
                        int block_size = 64;
                        int emax = *emax_pos++;
                        Fixed q[block_size];
                        for (uint i = 0; i < block_size; i++) {
                            q[ZFP::Codec3<ZFP::MemoryBitStream, T>::perm[i]] = Fixed::reinterpret(*transformed_pos++);
                        }
                        ZFP::Codec3<ZFP::MemoryBitStream, T>::inv_xform(q);
                        uint nx = std::min(dims[0] - x, 4u);
                        uint ny = std::min(dims[1] - y, 4u);
                        uint nz = std::min(dims[2] - z, 4u);
                        if (nx == 4 && ny == 4 && nz == 4) {
                            ZFP::Codec3<ZFP::MemoryBitStream, T>::inv_cast(q, p, s[0], s[1], s[2], emax);
                        } else {
                            ZFP::Codec3<ZFP::MemoryBitStream, T>::inv_cast(q, p, s[0], s[1], s[2], nx, ny, nz, emax);
                        }
                    }
                }
            }
        }
        return dec_data;
    }
    void save(uchar *&c) override {}

    void load(const uchar *&c, size_t &remaining_length) override {}

    std::pair<To, To> get_out_range() override { return std::make_pair(0, 0); }

};


template <class T, class To, uint N>
ZFPDecomposition<T, To, N> make_decomposition_zfp() {
    return ZFPDecomposition<T, To, N>();
}

}  // namespace SZ3

#endif
