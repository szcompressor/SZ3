// ALGO_BIOMD built with GCC's -funsafe-math-optimizations and -ffinite-math-only on an FMA target (each allowed on
// its own, unlike -ffast-math): a NaN in a frame that is not trailing fill must come back as NaN, the chunk stored
// losslessly. Exits 0 when it does.

#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>

#include "SZ3/api/sz.hpp"

int main() {
    const size_t frames = 2, atoms = 300, at = 3 * 17 + 1;
    std::vector<float> x(frames * atoms * 3);
    for (size_t i = 0; i < x.size(); i++) x[i] = 1.0f + 0.01f * float(i % 900) + 0.001f * float(i / 900);
    const unsigned nan_bits = 0x7fc00000u;  // set as bits, so the compiler cannot fold it away
    std::memcpy(&x[at], &nan_bits, sizeof(float));
    SZ3::Config conf(frames, atoms, 3);
    conf.cmprAlgo = SZ3::ALGO_BIOMD;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 5e-4;
    size_t size = 0;
    std::unique_ptr<char[]> cmp(SZ_compress(conf, x.data(), size));
    SZ3::Config dec;
    std::unique_ptr<float[]> out(SZ_decompress<float>(dec, cmp.get(), size));
    unsigned bits;
    std::memcpy(&bits, &out[at], sizeof(float));
    const bool kept = std::memcmp(out.get(), x.data(), x.size() * sizeof(float)) == 0;
    std::printf("algorithm %d, the NaN came back as bits %08x: %s\n", int(dec.cmprAlgo), bits, kept ? "kept" : "LOST");
    return kept ? 0 : 1;
}
