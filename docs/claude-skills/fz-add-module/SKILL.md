---
name: fz-add-module
description: Use when the user wants to add a brand-new module to fz — a new Decomposition, Quantizer, Encoder, or Lossless backend. Triggers include "add a quantizer", "new encoder", "implement a custom decomposition".
---

# Add a new module to fz

fz modules are header-only `.hpp` files. Add one file, implement the interface, that's the module.

## File location & naming

| Module kind | Folder | Naming |
|---|---|---|
| Decomposition | `include/SZ3/decomposition/` | `<Name>Decomposition.hpp` (e.g., `MGARDFusedDecomposition.hpp`) |
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

For pre-existing patterns: `NoPredictionDecomposition.hpp` (simplest), `MGARDFusedDecomposition.hpp` (transform + quantize), `SPERRFusedDecomposition.hpp` (transform + custom-state serialization).

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

Unique `uid` byte per quantizer (current values: `LinearQuantizer=0b10`, `FixedPointQuantizer=0b11`, `ScalarQuantizer=0b100`, `BitTruncationQuantizer=0b101`, `LogDomainQuantizer=0b111`, `GranularBitRoundQuantizer=0b1000`, `ClusterQuantizer=0b1001`, `LevelQuantizer=0b1010` — pick a free one). Saving/loading uid lets the bitstream catch mismatched-quantizer errors.

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

A module enters the library only after it passes the checks below. Each exists because the
corresponding defect has actually shipped here; passing means a module is known not to repeat one,
not that it is correct.

### Its header compiles on its own

```bash
python3 tools/test/check_headers.py --build-dir build
```

A header that only works because another header came first breaks the moment a translation unit
changes its include order, or a standard library stops providing something transitively. Run it
under both toolchains -- libstdc++ and libc++ provide different things transitively, so a clean run
under one proves nothing about the other:

```bash
python3 tools/test/check_headers.py --build-dir build --compiler g++
```

CI runs it on the Linux and macOS jobs.

### It passes its group contract

```bash
ctest --test-dir build -R SZ3_ModuleContract
```

`include/SZ3/testing/ModuleContract.hpp` holds the checks that apply to every module in a group. It
is installed with the library, so a module developed outside this tree runs the same gate:

```cpp
#include "SZ3/testing/ModuleContract.hpp"

TEST(MyModule, Contract) {
    SZ3_test::expectQuantizerContract<float>("MyQuantizer", [] { return MyQuantizer<float>(1e-3); });
}
```

In-tree modules are added to `tools/test/modules/test_module_contract.cpp`.

**Quantizer** — `expectQuantizerContract<T>(name, factory)`

- `get_out_range().first == 0`, and no bin below it. A signed bin type whose values fill the sign
  bit fails here.
- `recover()` reproduces what `quantize_and_overwrite()` wrote back, in encode order, after a
  `save()`/`load()` round-trip. A desynchronized unpredictable-value list fails here.
- `load()` on a truncated stream throws rather than reading past `remaining_length`.
- Determinism.
- Adversarial values throughout: `0`, `-0`, denormals, `min`, `max`, `±inf`, NaN.

Add `expectAbsErrorBound<T>(name, factory, eb)` when the module claims an absolute bound. The
bound is checked against the value actually written back, not against the idealized level.

**Decomposition** — `expectDecompositionContract<T>(name, conf, factory, eb, fields)`

- `get_out_range().first == 0`, and `.second` fits in the `int` that `preprocess_encode` takes
  (or is 0 for "no bin range"). Every bin sits inside the advertised range.
- `size_est()` covers what `save()` writes, verified with guard bytes.
- A **fresh instance loaded from that state** reconstructs within `eb`. Going through
  `save()`/`load()` between compress and decompress is what the compressor does, and it is where a
  decomposition that silently depends on its own compress-side state fails.
- `load()` on a truncated stream throws.
- Determinism.
- Adversarial fields: constant, ramp, single spike, alternating extremes, smooth.

Pass `eb = 0` only if the module claims no absolute bound. Pass explicit `fields` when the module
honours its bound on a restricted input regime, and say in a comment what that regime is —
`MGARDFusedDecomposition` is the one case in the tree.

**Encoder** — `expectEncoderContract<Bin>(name, factory, streams, stateNum)`

- `size_est()` covers what `encode()` writes, verified with guard bytes past the buffer.
- `encode()` does not assume a zeroed destination.
- Round-trip through `save`/`encode`/`load`/`decode` reproduces the bins exactly.
- Adversarial streams: one element, one run per sample, constant, alternating.

Pass a real `stateNum` for an encoder that cannot derive its own bin range — `ArithmeticEncoder`
sizes its frequency model from it and produces garbage when given 0.

**Lossless** — `expectLosslessContract(name, factory)`

- Round-trip through a generously sized destination reproduces the payload exactly.
- A destination smaller than the payload is refused, and nothing is written past it, verified with
  guard bytes.

**Preprocessor** — no contract yet: `concepts::PreprocessorInterface` declares no members (its only
method is commented out), so there is nothing to check against. Defining it is a prerequisite.

### It has its own tests

One `tools/test/modules/test_*.cpp` per module, covering whatever the group contract cannot know:
the module's error guarantee, its serialized format, its parameters, its documented edge cases.
CMake globs the directory, so a new file needs no build change.

### It works in a composition

A module that round-trips in isolation can still fail inside a pipeline. Add a case to
`tools/test/modules/test_composition.cpp` wiring it through `make_compressor_sz_generic` with a real
encoder and lossless stage. Most of the defects below were found this way, not by per-module tests.

### It is registered

- `include/SZ3/api/sz_dev.hpp` -- the developer header carries every module.
- `docs/site/modules.json` -- the metadata the module site reads. Add a `composition` block for any
  constraint a composition engine would need: an encoder that cannot derive its own bin range
  (`requires_state_num`), a decomposition that only works with one encoder (`requires_encoder`), a
  quantizer whose bins are not a countable domain (`out_range`), the regime an error bound holds in
  (`error_bound_regime`), or a pairing measured to be bad (`avoid_pairing`, with the number). Only
  record what a test or the code already enforces. `ctest -R SZ3_ModulesJson` checks the file
  against the tree and keeps the vocabulary closed.
- A `uid` distinct from every other module in its group, if the group serializes one.
  `ctest -R UidsAreDistinct` checks this by reading the first byte each `save()` writes.

### Clean build, both toolchains

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j$(nproc)
ctest --test-dir build
clang-format -i <changed files>
```

Zero warnings. CI builds Linux, macOS, Windows and Windows MinGW, and compares Linux/macOS digests;
PRs to `master` or `fz` additionally run the SDRBench datasets.

### Defect classes these checks came from

| Check | Defect it would have caught |
|---|---|
| header self-containment | `std::make_shared` without `<memory>`: broke the Linux build once a new test changed the include order |
| bin `>=` range start | `BitTruncationQuantizer<double>` emitted negative bins |
| `recover()` round-trip | `recover_unpred()` read past the unpredictable list with no bounds check |
| truncated `load()` throws | `MemoryUtil::read` guarded `remaining_length` with `assert`, which Release compiles out |
| absolute bound on written value | `LevelQuantizer`'s quadratic curve exceeded its bound by 150000x |
| `size_est()` covers `encode()` | `RunlengthEncoder` had no `size_est()`, overrunning the caller's buffer by 69 bytes |
| dirty destination buffer | `BitshuffleEncoder` only ORed bits into a `malloc`'d buffer |
| adversarial streams | `ArithmeticEncoder` sign-extended a shift, corrupting about half of all streams |
| `uid` distinctness | two quantizers shared `0b101`, so each accepted the other's stream |
| destination capacity | `Lossless_bypass::compress` memcpy'd the payload regardless of the capacity it was given |
| composition test | `TimeSeriesDecomposition` violated its bound 1.97x on the null-reference path |
| `size_est()` covers `save()` | six decompositions did not override `size_est()`, so the compressor sized its buffer from 0 |
| bound after `save()`/`load()` | `RegressionPredictor::load` debited `remaining_length` by the decoded rather than the encoded size |

## Quick syntax check before wiring

Before plumbing your module into an algo, sanity-check it compiles:

```bash
g++ -std=c++17 -I include -I build/include -fsyntax-only -x c++ - <<'EOF'
#include "SZ3/encoder/FooEncoder.hpp"
int main() { SZ3::FooEncoder<int> e; (void)e; return 0; }
EOF
```
