//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ3_LOSSLESS_ZSTD_HPP
#define SZ3_LOSSLESS_ZSTD_HPP

#include <memory>
#include <stdexcept>

#include "SZ3/def.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

// Do not replace these declarations with #include <zstd.h>. This is an installed public header,
// so that would oblige every consumer to have zstd.h on their include path, which under Homebrew,
// conda and Spack they do not -- they have Zstd's library, not its header.
// test_zstd_decl compiles these against the real header in both include orders, so a drift between
// them is a build failure. Define SZ3_USE_ZSTD_HEADER to include the header instead.
#ifdef SZ3_USE_ZSTD_HEADER
#include <zstd.h>
#else
#include <stddef.h>
extern "C" {
size_t ZSTD_compress(void *dst, size_t dstCapacity, const void *src, size_t srcSize, int compressionLevel);
size_t ZSTD_decompress(void *dst, size_t dstCapacity, const void *src, size_t compressedSize);
size_t ZSTD_compressBound(size_t srcSize);
unsigned ZSTD_isError(size_t code);
}
#endif

namespace SZ3 {
class Lossless_zstd : public concepts::LosslessInterface {
   public:
    Lossless_zstd() = default;

    Lossless_zstd(int comp_level) : compression_level(comp_level) {}

    /**
     * worst-case compressed size for srcLen bytes, for sizing the buffer handed to compress().
     * Keeps ZSTD_compressBound in this one file, which is what test_zstd_decl covers.
     */
    static size_t compress_bound(size_t srcLen) { return ZSTD_compressBound(srcLen); }

    /**
     * Attention
     * When dstCap is smaller than the space needed, ZSTD will not throw any errors.
     * Instead, it will write a portion of the compressed data to dst and stops.
     * This behavior is not desirable in SZ, as we need the whole compressed data for decompression.
     * Therefore, we need to check if the dst buffer (dstCap) is large enough for zstd
     */
    /**
     * compress data with lossless compressors
     * @param src  data to be compressed
     * @param srcLen length (in bytes) of the data to be compressed
     * @param dst compressed data
     * @param dstCap capacity (in bytes) for storing the compressed data
     * @return length (in bytes) of the data compressed
     */
    size_t compress(const uchar *src, size_t srcLen, uchar *dst, size_t dstCap) override {
        write(srcLen, dst);
        dstCap -= sizeof(size_t);  // reserve space for srcLen
        if (dstCap < ZSTD_compressBound(srcLen)) {
            throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        }
        size_t dstLen = ZSTD_compress(dst, dstCap, src, srcLen, compression_level);
        return dstLen + sizeof(size_t);
    }

    /**
     * reverse of compress(), decompress the data with lossless compressors
     * @param src data to be decompressed
     * @param srcLen length (in bytes) of that data
     * @param dst buffer to decompress into; when null on entry the callee allocates it with malloc()
     *            and the caller frees it
     * @param dstCap the capacity of dst, ignored when dst is null
     * @return length (in bytes) of the data decompressed
     */
    size_t decompress(const uchar *src, size_t srcLen, uchar *&dst, size_t dstCap) override {
        // The stream is a decompressed-size field followed by the zstd frame, all untrusted.
        if (srcLen < sizeof(size_t)) {
            throw std::out_of_range("SZ3 lossless: compressed data is smaller than the size header");
        }
        size_t dstLen = 0;
        read(dstLen, src);

        // malloc, because the caller frees what it gets back with free().
        std::unique_ptr<uchar, void (*)(void *)> owner(nullptr, &free);
        if (dst == nullptr) {
            owner.reset(static_cast<uchar *>(malloc(dstLen)));
            if (owner == nullptr) {
                throw std::runtime_error("SZ3 lossless: can not allocate the decompression buffer");
            }
        } else if (dstLen > dstCap) {
            throw std::out_of_range("SZ3 lossless: declared decompressed size exceeds the allowed capacity");
        }
        uchar *out = (dst != nullptr) ? dst : owner.get();

        // A short frame would leave the tail of the output uninitialized for the caller to read.
        size_t res = ZSTD_decompress(out, dstLen, src, srcLen - sizeof(size_t));
        if (ZSTD_isError(res) || res != dstLen) {
            throw std::runtime_error("SZ3 lossless: stream does not decompress to the size it declares");
        }

        dst = out;
        owner.release();
        return res;
    }

   private:
    int compression_level = 3;  // default setting of level is 3
};
}  // namespace SZ3
#endif  // SZ_LOSSLESS_ZSTD_HPP
