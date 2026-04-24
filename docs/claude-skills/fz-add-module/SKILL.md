---
name: fz-add-module
description: Use when the user wants to add a brand-new module to fz — a new Decomposition, Quantizer, Encoder, or Lossless backend. Triggers include "add a quantizer", "new encoder", "implement a custom decomposition".
---

# Add a new module to fz

fz modules are header-only `.hpp` files. Add one file, implement the interface, that's the module.

## File location & naming

| Module kind | Folder | Naming |
|---|---|---|
| Decomposition | `include/SZ3/decomposition/` | `<Name>Decomposition.hpp` (e.g., `MGARDDecomposition.hpp`) |
| Quantizer | `include/SZ3/quantizer/` | `<Name>Quantizer.hpp` (e.g., `LinearQuantizer.hpp`) |
| Encoder | `include/SZ3/encoder/` | `<Name>Encoder.hpp` (e.g., `HuffmanEncoder.hpp`) |
| Lossless | `include/SZ3/lossless/` | `Lossless_<name>.hpp` (e.g., `Lossless_zstd.hpp`) |

Class lives in `namespace SZ3`. Header guard: `SZ3_<UPPER_NAME>_HPP`.

## Save/load contract (all modules)

Pointer-based streaming. The buffer pointer is advanced by the call.

```cpp
void save(uchar*& c) const;                              // advances c, no length check
void load(const uchar*& c, size_t& remaining_length);    // advances c, decrements remaining_length
```

Use `write(value, c)` and `read(value, c, remaining_length)` from `include/SZ3/utils/MemoryUtil.hpp`. Never read past `remaining_length`; throw `std::runtime_error` on size mismatch.

## Decomposition skeleton

Interface at `include/SZ3/decomposition/Decomposition.hpp:21`. `Ti` = input float type, `To` = bin int type, `N` = data dim.

```cpp
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

template <class T, uint N>
class FooDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    std::vector<int> compress(const Config& conf, T* data) override {
        std::vector<int> bins(conf.num);
        // ... fill bins from data ...
        return bins;
    }
    T* decompress(const Config& conf, std::vector<int>& bins, T* dec_data) override {
        // ... reconstruct dec_data from bins ...
        return dec_data;
    }
    void save(uchar*& c) override { /* write any state */ }
    void load(const uchar*& c, size_t& remaining_length) override { /* read state */ }
    std::pair<int, int> get_out_range() override { return {0, max_bin}; }  // first MUST be 0
    size_t size_est() override { return /* metadata bytes */ ; }
};

}  // namespace SZ3
```

For pre-existing patterns: `NoPredictionDecomposition.hpp` (simplest), `MGARDDecomposition.hpp` (transform + quantize), `SPERRDecomposition.hpp` (transform + custom-state serialization).

## Quantizer skeleton

Interface at `include/SZ3/quantizer/Quantizer.hpp:17`. Inherit `QuantizerInterface<Ti, To>`.

```cpp
template <class T>
class FooQuantizer : public concepts::QuantizerInterface<T, int> {
   public:
    ALWAYS_INLINE int quantize_and_overwrite(T& data, T pred) override {
        // map data → integer bin, write reconstructed value back into data
    }
    ALWAYS_INLINE T recover(T pred, int q) override { /* inverse */ }
    int force_save_unpred(T ori) override { /* sentinel + push to unpred buffer */ }
    void save(uchar*& c) const override { /* write uid + state */ }
    void load(const uchar*& c, size_t& remaining_length) override { /* verify uid + read state */ }
    std::pair<int, int> get_out_range() const override { return {0, range_max}; }
};
```

Unique `uid` byte per quantizer (current values: `LinearQuantizer=0b10`, `QuadraticLevelQuantizer=0b11`, `FixedPointQuantizer=0b11` — pick a free one). Saving/loading uid lets the bitstream catch mismatched-quantizer errors.

Reference implementations: `LinearQuantizer.hpp` (default), `FixedPointQuantizer.hpp` (calibrated), `ScalarQuantizer.hpp` (asymmetric reconstruction tweaks).

## Encoder skeleton

Interface at `include/SZ3/encoder/Encoder.hpp:25`.

```cpp
template <class T>
class FooEncoder : public concepts::EncoderInterface<T> {
   public:
    void preprocess_encode(const std::vector<T>& bins, int stateNum) override { /* set up */ }
    size_t encode(const std::vector<T>& bins, uchar*& bytes) override { /* write payload, advance bytes */ }
    std::vector<T> decode(const uchar*& bytes, size_t targetLength) override { /* inverse */ }
    void save(uchar*& c) override { /* metadata */ }
    void load(const uchar*& c, size_t& remaining_length) override { /* metadata */ }
    void preprocess_decode() override {}
    void postprocess_encode() override {}
    void postprocess_decode() override {}
    size_t size_est() override { return /* estimated total payload + metadata bytes */ ; }
};
```

`size_est()` should be a realistic upper bound on the encoded payload size. `SZGenericCompressor` uses it (plus `sizeof(T) * bins.size()`) to pre-size its intermediate buffer; under-reporting can mask overruns. See `BitplaneEncoder.hpp` and `BitplaneRLEEncoder.hpp` for non-trivial examples.

References: `BypassEncoder.hpp` (simplest), `HuffmanEncoder.hpp` (default), `BitplaneEncoder.hpp` (transform-friendly).

## Lossless skeleton

Interface at `include/SZ3/lossless/Lossless.hpp:18`. The simplest existing impl is `Lossless_bypass.hpp`. Most users will not write a new one — `Lossless_zstd` covers the general case.

## Testing

Add a TEST_F to the matching file:

| Module | Test file |
|---|---|
| Encoder | `tools/test/modules/test_encoder.cpp` |
| Quantizer | `tools/test/modules/test_quantizer.cpp` |
| Lossless | `tools/test/modules/test_lossless.cpp` |
| Decomposition | (no per-decomp test file yet — write `tools/test/modules/test_decomposition.cpp` mirroring `test_quantizer.cpp`) |

Build + run tests:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j$(nproc)
ctest --test-dir build
```

## Quick syntax check before wiring

Before plumbing your module into an algo, sanity-check it compiles:

```bash
g++ -std=c++17 -I include -I build/include -fsyntax-only -x c++ - <<'EOF'
#include "SZ3/encoder/FooEncoder.hpp"
int main() { SZ3::FooEncoder<int> e; (void)e; return 0; }
EOF
```
