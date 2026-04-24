---
name: fz-add-algorithm
description: Use when the user wants to expose a new compression pipeline as an `ALGO_*` enum entry, dispatched by `Config::cmprAlgo`. Triggers include "add an algorithm to sz3", "new ALGO entry", "wire up a compressor", "register a new algorithm".
---

# Add a new ALGO_* to fz

A 4-edit checklist. Use an existing wiring file as your template — `SZAlgoMGARD.hpp` and `SZAlgoNopred.hpp` are the cleanest references.

## Step 1 — Create the wiring header

`include/SZ3/api/impl/SZAlgoFoo.hpp`:

```cpp
#ifndef SZ3_SZALGO_FOO_HPP
#define SZ3_SZALGO_FOO_HPP

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/FooDecomposition.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

template <class T, uint N>
size_t SZ_compress_Foo(Config& conf, T* data, uchar* cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_FOO);
    calAbsErrorBound(conf, data);  // populates conf.absErrorBound from REL/PSNR/L2NORM if set
    auto sz = make_compressor_sz_generic<T, N>(
        FooDecomposition<T, N>(/* ... */),
        HuffmanEncoder<int>(),
        Lossless_zstd());
    return sz->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_Foo(const Config& conf, const uchar* cmpData, size_t cmpSize, T* decData) {
    assert(conf.cmprAlgo == ALGO_FOO);
    auto cmpDataPos = cmpData;
    auto sz = make_compressor_sz_generic<T, N>(
        FooDecomposition<T, N>(/* ... */),
        HuffmanEncoder<int>(),
        Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

}  // namespace SZ3

#endif
```

The compress and decompress functions MUST construct an identically-configured pipeline — modules are stateful and rely on save/load round-tripping their internal state.

## Step 2 — Add the enum entry

In `include/SZ3/utils/Config.hpp`:

```cpp
enum ALGO {
    ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP, ALGO_NOPRED, ALGO_LOSSLESS,
    ALGO_BIOMD, ALGO_BIOMDXTC, ALGO_SVD, ALGO_ZFP, ALGO_SPERR, ALGO_MGARD,
    ALGO_FOO   // <-- append at end (preserves serialized values for older compressed files)
};

const std::map<std::string, ALGO> ALGO_MAP = {
    /* ... */
    {"ALGO_FOO", ALGO_FOO},   // <-- so INI files can name your algo
};
```

**Always append** to the enum — never insert in the middle. Old compressed files store the integer value, so reordering breaks decompression.

## Step 3 — Wire into the dispatcher

In `include/SZ3/api/impl/SZDispatcher.hpp`:

Add the include near the top (alphabetical with the other `SZAlgo*` includes):
```cpp
#include "SZ3/api/impl/SZAlgoFoo.hpp"
```

Add a compress branch (mirroring the existing pattern, around line 108):
```cpp
} else if (conf.cmprAlgo == ALGO_FOO) {
    if constexpr (std::is_floating_point<T>::value && N >= 1 && N <= 3) {
        cmpSize = SZ_compress_Foo<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
    } else {
        throw std::invalid_argument("FOO algorithm supports 1D/2D/3D floating-point data only.");
    }
}
```

Add a matching decompress branch in `SZ_decompress_dispatcher` (around line 165 in the same file).

Adjust the `if constexpr` guards to match what your algo actually supports (drop the `N >= 1 && N <= 3` if it works for any N; drop the floating-point guard if it works on integer types too).

## Step 4 — Build + smoke-test

```bash
cmake --build build --target sz3 -j$(nproc)

# CLI smoke (use the in-tree dataset; -3 takes dims FASTEST-FIRST)
printf '[GlobalSettings]\nCmprAlgo = ALGO_FOO\nErrorBoundMode = ABS\nAbsErrorBound = 1e-3\n' > /tmp/foo.cfg
build/tools/sz3/sz3 -f -i tools/sz3/testfloat_8_8_128.dat -z /tmp/foo.sz \
  -3 128 8 8 -M ABS 1e-3 -c /tmp/foo.cfg
build/tools/sz3/sz3 -f -s /tmp/foo.sz -x /tmp/foo.out.f32 \
  -i tools/sz3/testfloat_8_8_128.dat -3 128 8 8 -a
```

`-a` prints CR + max-abs error + PSNR. Confirm max abs error ≤ your bound.

For a multi-bound benchmark + comparison against another algo, see the `fz-bench-multibound` skill.

## Optional — common guards

| Constraint | Snippet |
|---|---|
| Floating-point only | `if constexpr (std::is_floating_point<T>::value) { ... } else throw ...;` |
| 3D only | `if constexpr (N == 3) { ... } else throw ...;` |
| Specific platform exclusion | `#if !defined(__MINGW32__)` around the dispatch branch (see how SPERR is gated in `SZDispatcher.hpp:30-31`) |

## Don't forget

- **Header guards** — distinct `SZ3_SZALGO_<NAME>_HPP`
- **`namespace SZ3`** — match the rest of the codebase
- **Both compress AND decompress** branches in the dispatcher — easy to forget one
