---
name: fz-compose-pipeline
description: Use when the user wants to combine existing fz modules into a custom compression pipeline, or asks "how do I plug in a different encoder/quantizer", "compose modules", "build a custom compressor from fz modules".
---

# Compose a pipeline from fz modules

A pipeline composition lives in a single `SZAlgo*.hpp` header under `include/SZ3/api/impl/` and uses the factory `make_compressor_sz_generic<T, N>(decomp, encoder, lossless)`.

## The factory

```cpp
// include/SZ3/compressor/SZGenericCompressor.hpp:151
template <class T, uint N, class Decomposition, class Encoder, class Lossless>
auto make_compressor_sz_generic(Decomposition d, Encoder e, Lossless l) {
    return std::make_shared<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>>(d, e, l);
}
```

## Type-matching contract

The integer "bin" type must agree across the Decomposition and Encoder:

```
Decomposition : DecompositionInterface<Ti, To, N>
Encoder       : EncoderInterface<To>
```

Common combos:

| Decomposition output | Pair with Encoder |
|---|---|
| `int` (default — Lorenzo, Interp, NoPred, MGARD) | `HuffmanEncoder<int>`, `BitplaneEncoder<int>`, `BitplaneRLEEncoder<int>`, `RunlengthEncoder<int>` |
| `int64_t` (SPERR, FixedPoint-style) | `HuffmanEncoder<int64_t>`, `BitplaneEncoder<int64_t>`, `BitplaneRLEEncoder<int64_t>`, `SPERREncoder<int64_t, N>` |

## Output-range contract

`SZGenericCompressor` enforces `decomposition.get_out_range().first == 0` at `SZGenericCompressor.hpp:78`. Quantizers that produce signed bins (e.g., `LinearQuantizer`) handle this by shifting to `[0, 2*radius]`. If you write a new quantizer/decomposition with a signed natural range, you must shift before returning.

## Example 1 — quantize-only baseline

```cpp
// include/SZ3/api/impl/SZAlgoNopred.hpp
auto sz = make_compressor_sz_generic<T, N>(
    make_decomposition_noprediction<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
    HuffmanEncoder<int>(),
    Lossless_zstd());
```

`NoPredictionDecomposition` (`include/SZ3/decomposition/NoPredictionDecomposition.hpp`) wraps a quantizer and applies it to every value with `pred=0`.

## Example 2 — MGARD multigrid + bit-plane

```cpp
// include/SZ3/api/impl/SZAlgoMGARD.hpp
const double coef_eb = conf.absErrorBound / mgard_coef_eb_scale(N);  // L∞ amplification factor
auto sz = make_compressor_sz_generic<T, N>(
    make_decomposition_mgard<T, N>(conf, LinearQuantizer<T>(coef_eb, conf.quantbinCnt / 2)),
    BitplaneEncoder<int>(),
    Lossless_zstd());
```

Here `MGARDFusedDecomposition` runs the multigrid forward transform first, then drives `LinearQuantizer` over the resulting coefficients. The bit-plane encoder is a swap-in alternative to Huffman that often wins on transformed coefficient distributions.

## Example 3 — SPERR (bypass lossless)

```cpp
// include/SZ3/api/impl/SZAlgoSPERR.hpp
auto sz = make_compressor_sz_generic<T, N>(
    SPERRFusedDecomposition<T, N>(),
    SPERREncoder<int64_t, N>(make_sperr_dims_from_conf(conf)),
    Lossless_bypass());  // SPERR's SPECK encoder already produces a compressed bitstream
```

When the encoder already produces a compressed bitstream, use `Lossless_bypass` to avoid double-compressing.

## Where to put a new pipeline

If your composition is permanent (you want it as `ALGO_FOO`), wire it via the `fz-add-algorithm` skill — that's a 4-edit checklist. If you just want to experiment locally, you can call `make_compressor_sz_generic` directly from `tools/sz3/sz3_customized_demo.cpp` (see `tools/sz3/demo/sz3_customized_demo.cpp` for the pattern).

## Things that bite

- **Mismatched integer types**: passing `HuffmanEncoder<int>` to a decomposition with `To = int64_t` fails to compile with a confusing template error from inside `SZGenericCompressor::compress`. Match types explicitly.
- **`Lossless_bypass` semantics**: it does NOT no-op — it copies bytes and writes the size header. Use it only when the encoder's output is already compressed.
- **Quantizer state lifecycle**: a quantizer constructed once is used for one full compress (or one full decompress). Don't share an instance between calls.
