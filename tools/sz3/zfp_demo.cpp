/**
 * This is a demo for creating your own compressor in SZ3.
 * Essentially, there are two main steps: first, you need to implement your own compressor.
 * Please check out the 4 examples provided in the main function for implementing your own compressor.
 *
 * Second, you need to run your compressor with SZ3 executable or a new executable.
 * This demon code is a new executable, so you have control over the parameter parsing, IO, etc.
 * If you want to use the SZ3 executable to run your own compressor, please follow the steps below:
 * 1. Add a new ALGO in SZ3 Config
 * 2. Add a new hpp file in "include/SZ3/api/impl/" to assemble your compressor, one example easy to follow is
 * "include/SZ3/api/impl/SZAlgoNopred.hpp"
 * 3. Add the corresponding code in "include/SZ3/api/impl/SZDispatcher.hpp" to dispatch the new compressor.
 * 4. When executing the SZ3 executable, use -c to specify the config file. The config file should use the new ALGO.
 *  Example config file is "tools/sz3/sz3.config".
 */

#include "SZ3/api/sz.hpp"
#include "SZ3/utils/zfp/zfpcodec.h"
#include "SZ3/utils/zfp/zfpcodec1.h"
#include "SZ3/utils/zfp/zfpcodec2.h"
#include "SZ3/utils/zfp/zfpcodec3.h"

using namespace SZ3;

template <class T, uint N>
void SZ3_interpolation_compress(Config &conf, T *data, char *dst, size_t &outSize) {
    calAbsErrorBound(conf, data);

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    sz->compress(conf, data, reinterpret_cast<uchar *>(dst), outSize);
}

template <class T, uint N>
void SZ3_interpolation_decompress(const Config &conf, const char *cmpData, size_t cmpSize, T *decData) {
    auto cmpDataPos = reinterpret_cast<const uchar *>(cmpData);
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

template <class T, uint N>
class MyDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    MyDecomposition(const Config &conf) : num(conf.num), quantizer(conf.absErrorBound) {}

    std::vector<int> compress(const Config &conf, T *data) override {
        // write your own logic here
        std::vector<int> output(num);
        for (size_t i = 0; i < num; i++) {
            output[i] = quantizer.quantize_and_overwrite(data[i], 0);  // this will replace the value of data[i]
        }
        return output;
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        // write your own logic here
        for (size_t i = 0; i < num; i++) {
            dec_data[i] = quantizer.recover(0, quant_inds[i]);
        }
        return dec_data;
    }

    void save(uchar *&c) override {
        write(num, c);
        quantizer.save(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        read(num, c, remaining_length);
        quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return std::make_pair(0, 0); }

   private:
    size_t num;
    LinearQuantizer<T> quantizer;
};

template <class Int, uint N>
class ZFPEncoder : public concepts::EncoderInterface<Int> {
   public:
    using UInt = std::make_unsigned_t<Int>;
    static const uint ebits = sizeof(Int) == 4 ? 8 : 11;  // number of exponent bits
    static const int ebias = (1 << (ebits - 1)) - 1;      // floating-point exponent bias

    ZFPEncoder(const Config &conf) {
        int expmin = INT_MIN;
        if (conf.absErrorBound > 0) {
            frexp(conf.absErrorBound, &expmin);
            expmin--;
        }
        minexp = std::max(expmin, std::numeric_limits<Int>::min_exponent - std::numeric_limits<Int>::digits);
    }

    size_t encode(const std::vector<Int> &data, uchar *&bytes) override {
        MemoryBitStream stream;
        stream.open(bytes, data.size() * sizeof(Int));

        size_t n_blocks = data[0];
        uint maxprec = CHAR_BIT * sizeof(Int);
        auto emax_pos = &data[1];
        auto int_pos = &data[1 + n_blocks];
        int block_size = N == 3 ? 64 : (N == 2 ? 16 : 4);

        stream.write(n_blocks, sizeof(size_t) * 8);
        for (auto i = 0; i < n_blocks; i++) {
            stream.write(*emax_pos + ebias, ebits);
            if constexpr (N == 1) {
                uint precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 4)));
                IntCodec04<MemoryBitStream, Int, UInt>::encode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                emax_pos++;
            } else if constexpr (N == 2) {
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 6)));
                IntCodec16<MemoryBitStream, Int, UInt>::encode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                emax_pos++;
            } else if constexpr (N == 3) {
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 8)));
                IntCodec64<MemoryBitStream, Int, UInt>::encode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                emax_pos++;
            }
        }
        stream.flush();
        bytes += stream.size();
        printf("outSize=%zu\n", stream.size());
        return stream.size();
    }

    std::vector<Int> decode(const uchar *&bytes, size_t targetLength) override {
        MemoryBitStream stream;
        stream.open(const_cast<uchar *>(bytes), targetLength);

        size_t n_blocks = stream.read(sizeof(size_t) * 8);
        uint maxprec = CHAR_BIT * sizeof(Int);
        int block_size = N == 3 ? 64 : (N == 2 ? 16 : 4);

        std::vector<int> output(1 + n_blocks + n_blocks * block_size);
        output[0] = n_blocks;
        auto emax_pos = &output[1];
        auto int_pos = &output[1 + n_blocks];

        for (auto i = 0; i < n_blocks; i++) {
            *emax_pos = stream.read(ebits) - ebias;
            if constexpr (N == 1) {
                uint precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 4)));
                IntCodec04<MemoryBitStream, Int, UInt>::decode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                emax_pos++;
            } else if constexpr (N == 2) {
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 6)));
                IntCodec16<MemoryBitStream, Int, UInt>::decode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                emax_pos++;
            } else if constexpr (N == 3) {
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 8)));
                IntCodec64<MemoryBitStream, Int, UInt>::decode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                emax_pos++;
            }
        }
        return output;
    }

    int minexp = INT_MIN;
    void preprocess_encode(const std::vector<Int> &bins, int stateNum) override {}
    void postprocess_encode() override {}

    void preprocess_decode() override {}

    void postprocess_decode() override {}

    void save(uchar *&c) override {}

    void load(const uchar *&c, size_t &remaining_length) override {}
};

template <class T, class To, uint N>
class ZFPDecomposition {
   public:
    typedef typename ScalarTraits<T>::Fixed Fixed;

    std::vector<To> compress(Config &conf, const T *data) {
        calAbsErrorBound<T>(conf, data);

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
                    emax = ZFP::Codec1<MemoryBitStream, T>::fwd_cast(q, p, s[0]);
                    // perform block transform
                    ZFP::Codec1<MemoryBitStream, T>::fwd_xform(q);
                } else {
                    emax = ZFP::Codec1<MemoryBitStream, T>::fwd_cast(q, p, s[0], nx);
                    ZFP::Codec1<MemoryBitStream, T>::fwd_xform(q, nx);
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
                        emax = ZFP::Codec2<MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1]);
                        // perform block transform
                        ZFP::Codec2<MemoryBitStream, T>::fwd_xform(q);
                    } else {
                        emax = ZFP::Codec2<MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], nx, ny);
                        ZFP::Codec2<MemoryBitStream, T>::fwd_xform(q, nx, ny);
                    }
                    // reorder and convert to integer
                    *output_emax++ = emax;
                    for (uint i = 0; i < block_size; i++) {
                        *output_data++ = q[ZFP::Codec2<MemoryBitStream, T>::perm[i]].reinterpret();
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
                            emax = ZFP::Codec3<MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], s[2]);
                            // perform block transform
                            ZFP::Codec3<MemoryBitStream, T>::fwd_xform(q);
                        } else {
                            emax = ZFP::Codec3<MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], s[2], nx, ny, nz);
                            ZFP::Codec3<MemoryBitStream, T>::fwd_xform(q, nx, ny, nz);
                        }
                        *output_emax++ = emax;
                        // reorder and convert to integer
                        for (uint i = 0; i < block_size; i++) {
                            *output_data++ = q[ZFP::Codec3<MemoryBitStream, T>::perm[i]].reinterpret();
                        }
                    }
                }
            }
        }
        return output;
    }
    T *decompress(const Config &conf, std::vector<To> &transformed, T *dec_data) {
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
                ZFP::Codec1<MemoryBitStream, T>::inv_xform(q);
                uint nx = std::min(dims[0] - x, 4u);
                if (nx == 4) {
                    ZFP::Codec1<MemoryBitStream, T>::inv_cast(q, p, s[0], emax);
                } else {
                    ZFP::Codec1<MemoryBitStream, T>::inv_cast(q, p, s[0], nx, emax);
                }
            }
        } else if constexpr (N == 2) {
            for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
                for (uint x = 0; x < dims[0]; x += 4, p += 4) {
                    int block_size = 16;
                    int emax = *emax_pos++;
                    Fixed q[block_size];
                    for (uint i = 0; i < block_size; i++) {
                        q[ZFP::Codec2<MemoryBitStream, T>::perm[i]] = Fixed::reinterpret(*transformed_pos++);
                    }
                    ZFP::Codec2<MemoryBitStream, T>::inv_xform(q);
                    uint nx = std::min(dims[0] - x, 4u);
                    uint ny = std::min(dims[1] - y, 4u);
                    if (nx == 4 && ny == 4) {
                        ZFP::Codec2<MemoryBitStream, T>::inv_cast(q, p, s[0], s[1], emax);
                    } else {
                        ZFP::Codec2<MemoryBitStream, T>::inv_cast(q, p, s[0], s[1], nx, ny, emax);
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
                            q[ZFP::Codec3<MemoryBitStream, T>::perm[i]] = Fixed::reinterpret(*transformed_pos++);
                        }
                        ZFP::Codec3<MemoryBitStream, T>::inv_xform(q);
                        uint nx = std::min(dims[0] - x, 4u);
                        uint ny = std::min(dims[1] - y, 4u);
                        uint nz = std::min(dims[2] - z, 4u);
                        if (nx == 4 && ny == 4 && nz == 4) {
                            ZFP::Codec3<MemoryBitStream, T>::inv_cast(q, p, s[0], s[1], s[2], emax);
                        } else {
                            ZFP::Codec3<MemoryBitStream, T>::inv_cast(q, p, s[0], s[1], s[2], nx, ny, nz, emax);
                        }
                    }
                }
            }
        }
        return dec_data;
    }
};

template <class T, uint N>
void SZ3_customized_compress1(Config &conf, T *data, char *dst, size_t &outSize) {
    typedef typename ScalarTraits<T>::Fixed Fixed;
    typedef typename ScalarTraits<T>::Int Int;
    typedef typename ScalarTraits<T>::UInt UInt;

    static const uint ebits = ScalarTraits<T>::ebits;  // number of exponent bits
    static const int ebias = (1 << (ebits - 1)) - 1;   // floating-point exponent bias

    calAbsErrorBound(conf, data);

    MemoryBitStream stream;
    stream.open(dst, outSize);
    const T *p = data;

    int expmin = INT_MIN;
    if (conf.absErrorBound > 0) {
        frexp(conf.absErrorBound, &expmin);
        expmin--;
    }
    uint maxprec = CHAR_BIT * sizeof(Int);
    int minexp = std::max(expmin, std::numeric_limits<T>::min_exponent - std::numeric_limits<T>::digits);

    std::vector<uint> dims(conf.dims.begin(), conf.dims.end());
    std::vector<uint> s(N, 1);  // stride
    for (uint i = 1; i < N; i++) {
        s[i] = s[i - 1] * dims[i - 1];
    }
    if constexpr (N == 1) {
        for (auto x = 0; x < dims[0]; x += 4, p += 4) {
            int block_size = 4;
            Fixed q[block_size];
            int emax = 0;
            uint nx = std::min(dims[0] - x, 4u);
            if (nx == 4) {
                // convert to fixed-point
                emax = ZFP::Codec1<MemoryBitStream, T>::fwd_cast(q, p, s[0]);
                // perform block transform
                ZFP::Codec1<MemoryBitStream, T>::fwd_xform(q);
            } else {
                emax = ZFP::Codec1<MemoryBitStream, T>::fwd_cast(q, p, s[0], nx);
                ZFP::Codec1<MemoryBitStream, T>::fwd_xform(q, nx);
            }
            // reorder and convert to integer
            Int buffer[block_size];
            for (uint i = 0; i < block_size; i++) {
                buffer[i] = q[i].reinterpret();
            }
            stream.write(emax + ebias, ebits);
            uint precision = std::min(maxprec, static_cast<uint>(std::max(0, emax - minexp + 4)));
            IntCodec04<MemoryBitStream, Int, UInt>::encode(stream, buffer, 0, UINT_MAX, precision);
        }
    } else if constexpr (N == 2) {
        for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
            for (uint x = 0; x < dims[0]; x += 4, p += 4) {
                int block_size = 16;
                Fixed q[block_size];
                int emax = 0;
                uint nx = std::min(dims[0] - x, 4u);
                uint ny = std::min(dims[1] - y, 4u);
                if (nx == 4 && ny == 4) {
                    // convert to fixed-point
                    emax = ZFP::Codec2<MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1]);
                    // perform block transform
                    ZFP::Codec2<MemoryBitStream, T>::fwd_xform(q);
                } else {
                    emax = ZFP::Codec2<MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], nx, ny);
                    ZFP::Codec2<MemoryBitStream, T>::fwd_xform(q, nx, ny);
                }
                // reorder and convert to integer
                Int buffer[block_size];
                for (uint i = 0; i < block_size; i++)
                    buffer[i] = q[ZFP::Codec2<MemoryBitStream, T>::perm[i]].reinterpret();
                stream.write(emax + ebias, ebits);
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, emax - minexp + 6)));
                IntCodec16<MemoryBitStream, Int, UInt>::encode(stream, buffer, 0, UINT_MAX, precision);
            }
        }
    } else if constexpr (N == 3) {
        for (uint z = 0; z < dims[2]; z += 4, p += 4 * dims[0] * (dims[1] - (dims[1] + 3) / 4)) {
            for (uint y = 0; y < dims[1]; y += 4, p += 4 * (dims[0] - (dims[0] + 3) / 4)) {
                for (uint x = 0; x < dims[0]; x += 4, p += 4) {
                    int block_size = 64;
                    Fixed q[block_size];
                    int emax = 0;
                    uint nx = std::min(dims[0] - x, 4u);
                    uint ny = std::min(dims[1] - y, 4u);
                    uint nz = std::min(dims[2] - z, 4u);
                    if (nx == 4 && ny == 4 && nz == 4) {
                        // convert to fixed-point
                        emax = ZFP::Codec3<MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], s[2]);
                        // perform block transform
                        ZFP::Codec3<MemoryBitStream, T>::fwd_xform(q);
                    } else {
                        emax = ZFP::Codec3<MemoryBitStream, T>::fwd_cast(q, p, s[0], s[1], s[2], nx, ny, nz);
                        ZFP::Codec3<MemoryBitStream, T>::fwd_xform(q, nx, ny, nz);
                    }
                    // reorder and convert to integer
                    Int buffer[block_size];
                    for (uint i = 0; i < block_size; i++)
                        buffer[i] = q[ZFP::Codec3<MemoryBitStream, T>::perm[i]].reinterpret();
                    stream.write(emax + ebias, ebits);
                    auto precision = std::min(maxprec, static_cast<uint>(std::max(0, emax - minexp + 8)));
                    IntCodec64<MemoryBitStream, Int, UInt>::encode(stream, buffer, 0, UINT_MAX, precision);
                }
            }
        }
    }
    stream.flush();
    outSize = stream.size();
    printf("outSize=%zu\n", outSize);
}

template <class T, uint N>
void SZ3_customized_compress(Config &conf, T *data, char *dst, size_t &outSize) {
    calAbsErrorBound(conf, data);

    MemoryBitStream stream;
    stream.open(dst, outSize);
    int expmin = INT_MIN;
    if (conf.absErrorBound > 0) {
        frexp(conf.absErrorBound, &expmin);
        expmin--;
    }
    const T *p = data;

    if constexpr (N == 1) {
        ZFP::Codec1<MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
        for (auto x = 0; x < conf.dims[0]; x += 4, p += 4) {
            codec.encode(p, 1, codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u)));
        }
    } else if constexpr (N == 2) {
        ZFP::Codec2<MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
        for (uint y = 0; y < conf.dims[1]; y += 4, p += 4 * (conf.dims[0] - (conf.dims[0] + 3) / 4)) {
            for (uint x = 0; x < conf.dims[0]; x += 4, p += 4) {
                codec.encode(p, 1, conf.dims[0],
                             codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u),
                                        std::min(static_cast<uint>(conf.dims[1] - y), 4u)));
            }
        }
    } else if constexpr (N == 3) {
        ZFP::Codec3<MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
        for (uint z = 0; z < conf.dims[2]; z += 4, p += 4 * conf.dims[0] * (conf.dims[1] - (conf.dims[1] + 3) / 4)) {
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
    outSize = stream.size();
    printf("outSize=%zu\n", outSize);
}

template <class T, uint N>
void SZ3_customized_decompress(const Config &conf, const char *cmpData, size_t cmpSize, T *decData) {
    int expmin = INT_MIN;
    if (conf.absErrorBound > 0) {
        frexp(conf.absErrorBound, &expmin);
        expmin--;
    }
    MemoryBitStream stream;
    stream.open((void *)cmpData, cmpSize);
    float *p = decData;

    if constexpr (N == 1) {
        ZFP::Codec1<MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
        for (auto x = 0; x < conf.dims[0]; x += 4, p += 4) {
            codec.decode(p, 1, codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u)));
        }
    } else if constexpr (N == 2) {
        ZFP::Codec2<MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
        for (uint y = 0; y < conf.dims[1]; y += 4, p += 4 * (conf.dims[0] - (conf.dims[0] + 3) / 4)) {
            for (uint x = 0; x < conf.dims[0]; x += 4, p += 4) {
                codec.decode(p, 1, conf.dims[0],
                             codec.dims(std::min(static_cast<uint>(conf.dims[0] - x), 4u),
                                        std::min(static_cast<uint>(conf.dims[1] - y), 4u)));
            }
        }
    } else if constexpr (N == 3) {
        ZFP::Codec3<MemoryBitStream, T> codec(stream, 0, UINT_MAX, 0, expmin);
        for (uint z = 0; z < conf.dims[2]; z += 4, p += 4 * conf.dims[0] * (conf.dims[1] - (conf.dims[1] + 3) / 4)) {
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
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "SZ v" << SZ3_VER << std::endl;
        std::cout << "usage: " << argv[0] << " data_file -num_dim dim0 .. dimn ABS" << std::endl;
        std::cout << "example: " << argv[0] << " qmcpack.dat -3 33120 69 69 1e-3" << std::endl;
        return 0;
    }

    std::string src_file_name = argv[1];

    // read dimensions
    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    // if (dim != 3) {
    //     printf("This demo only supports 3D data.\n");
    //     return 0;
    // }
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    // create config with data dimensions
    Config conf;
    conf.setDims(dims.begin(), dims.end());

    // prepare input and output buffers
    std::vector<float> input_data(conf.num);
    SZ3::readfile<float>(argv[1], conf.num, input_data.data());
    std::vector<float> input_data_copy(input_data);
    std::vector<float> dec_data(conf.num);
    auto dec_data_pos = dec_data.data();
    std::vector<char> cmpData(conf.num * sizeof(float) * 2);
    size_t cmpSize = cmpData.size();

    // set error bound (you can also use relative error bound EB_REL etc.)
    conf.errorBoundMode = EB_ABS;
    conf.absErrorBound = atof(argv[argp++]);

    Timer timer;
    // ================================================================================
    // There are 4 customization examples, you just need to pick up one.

    // Example 1: using SZ3 API defined in SZ3/API/sz.hpp
    // SZ3 has some built-in compressors with fixed modules, for example, the interpolation compressor uses
    // interpolation decomposition, linear quantizer, and huffman encoder.
    // The API takes the SZ3::Config as parameter which allows some basic customization, such as interpolation
    // method (linear or cubic) and the quantization levels.
    // SZ_compress<float>(conf, input_data.data(), cmpData.data(), cmpSize);
    // SZ_decompress<float>(conf, cmpData.data(), cmpSize, dec_data_pos);

    // Example 2: using provided modules to assemble a generic compressor
    // The generic compressor follow the decomposition -> encoding -> lossless pipeline.
    // You can use any decomposition, encoder, and lossless modules to form such a compressor, including bypass
    // modules. This example assembles the interpolation compressor with interpolation decomposition, linear
    // quantizer, and huffman encoder. You can change any of the modules to get a new compressor that is not
    // built-in in SZ3.
    // SZ3_interpolation_compress<float, 3>(conf, input_data.data(), cmpData.data(), cmpSize);
    // SZ3_interpolation_decompress<float, 3>(conf, cmpData.data(), cmpSize, dec_data_pos);

    // Example 3: using new modules to assemble a generic compressor
    // This is an extension of Example 2, which uses a custom decomposition module instead of built-in ones.
    // MyDecomposition is the new decomposition module, which implements the DecompositionInterface.
    // We put the source code of MyDecomposition here for demonstration purpose.
    // You should put your new modules in the corresponding module folders (include/SZ3/decomposition for Decomposition
    // modules).

    if (dim == 1) {
        timer.start();
        ZFPDecomposition<float, int32_t, 1> zfp_transform;
        ZFPEncoder<int32_t, 1> zfp_encoder(conf);
        std::vector<uchar> cmp_data(conf.num * sizeof(float));
        auto cmp_data_pos = cmp_data.data();
        auto transformed = zfp_transform.compress(conf, input_data.data());
        zfp_encoder.encode(transformed, cmp_data_pos);
        cmpSize = cmp_data_pos - cmp_data.data();
        timer.stop("decompression");

        timer.start();
        const uchar *cmp_data_pos1 = cmp_data.data();
        transformed = zfp_encoder.decode(cmp_data_pos1, cmpSize);
        zfp_transform.decompress(conf, transformed, dec_data_pos);
        timer.stop("decompression");

        // timer.start();
        // SZ3_customized_compress<float, 1>(conf, input_data.data(), cmpData.data(), cmpSize);
        // timer.stop("decompression");
        //
        // timer.start();
        // SZ3_customized_decompress<float, 1>(conf, cmpData.data(), cmpSize, dec_data_pos);
        // timer.stop("decompression");
    } else if (dim == 2) {
        SZ3_customized_compress<float, 2>(conf, input_data.data(), cmpData.data(), cmpSize);
        SZ3_customized_decompress<float, 2>(conf, cmpData.data(), cmpSize, dec_data_pos);
    } else if (dim == 3) {
        timer.start();
        ZFPDecomposition<float, int32_t, 3> zfp_transform;
        ZFPEncoder<int32_t, 3> zfp_encoder(conf);
        std::vector<uchar> cmp_data(conf.num * sizeof(float));
        auto cmp_data_pos = cmp_data.data();
        auto transformed = zfp_transform.compress(conf, input_data.data());
        zfp_encoder.encode(transformed, cmp_data_pos);
        cmpSize = cmp_data_pos - cmp_data.data();
        timer.stop("decompression");

        timer.start();
        const uchar *cmp_data_pos1 = cmp_data.data();
        transformed = zfp_encoder.decode(cmp_data_pos1, cmpSize);
        zfp_transform.decompress(conf, transformed, dec_data_pos);
        timer.stop("decompression");

        // timer.start();
        // SZ3_customized_compress<float, 3>(conf, input_data.data(), cmpData.data(), cmpSize);
        // timer.stop("decompression");
        //
        // timer.start();
        // SZ3_customized_decompress<float, 3>(conf, cmpData.data(), cmpSize, dec_data_pos);
        // timer.stop("decompression");

    }

    // Example 4: assemble a specialized compressor
    // If your compressor doesn't follow the decomposition -> encoding -> lossless pipeline,
    // you can implement your own compressor by implementing the CompressorInterface
    // Please check out "SZ3/compressor/specialized/SZExaaltCompressor.hpp". It contains two separate encoding process.

    // ================================================================================

    // In the end, you can use the SZ3 API to verify the decompressed data.
    SZ3::verify<float>(input_data_copy.data(),  dec_data.data(), conf.num);

    printf("compression ratio = %.3f\n", (float)(conf.num * sizeof(float)) / cmpSize);
    // write compressed data and decompressed data to the same folder as original data
    writefile((src_file_name + ".demo").data(), cmpData.data(), cmpSize);
    writefile((src_file_name + ".demo.out").data(), dec_data.data(), conf.num);

    return 0;
}
