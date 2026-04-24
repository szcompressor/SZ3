---
name: fz-overview
description: Use when the user asks what fz is, lists/discovers algorithms or modules, or needs a quick architecture map. Triggers include "what algorithms does fz support", "what modules are in fz", "fz architecture", "fz pipeline".
---

# fz overview

fz is a header-only C++17 library that vendors several SOTA error-bounded scientific compressors (SZ3, ZFP, SPERR, MGARD) into one repo and exposes them through a **modular pipeline** so developers can compose new compressors from a shared pool of building blocks.

## The 4-stage pipeline

Every composition built on `SZGenericCompressor` follows this shape:

```
Input → Decomposition → (Quantizer) → Encoder → Lossless → Output
```

Module interfaces (file references):

| Stage | Interface | Header |
|---|---|---|
| Decomposition | `concepts::DecompositionInterface<Ti, To, N>` | `include/SZ3/decomposition/Decomposition.hpp:21` |
| Quantizer | `concepts::QuantizerInterface<Ti, To>` | `include/SZ3/quantizer/Quantizer.hpp:17` |
| Encoder | `concepts::EncoderInterface<T>` | `include/SZ3/encoder/Encoder.hpp:25` |
| Lossless | `concepts::LosslessInterface` | `include/SZ3/lossless/Lossless.hpp:18` |

Composition factory: `make_compressor_sz_generic<T, N>(decomp, encoder, lossless)` at `include/SZ3/compressor/SZGenericCompressor.hpp:151`.

**Output-range contract**: `decomposition.get_out_range().first` MUST be `0` (enforced at `SZGenericCompressor.hpp:78`). Quantizers that produce signed bins must shift to non-negative.

## Available algorithms (`enum ALGO`)

Defined in `include/SZ3/utils/Config.hpp`. Each has a `SZAlgo*.hpp` wiring file in `include/SZ3/api/impl/`:

| Algo | Wiring | Notes |
|---|---|---|
| `ALGO_LORENZO_REG` | `SZAlgoLorenzoReg.hpp` | Blockwise Lorenzo + regression |
| `ALGO_INTERP_LORENZO` *(default)* | `SZAlgoInterp.hpp` | Auto-tuned interpolation + Lorenzo |
| `ALGO_INTERP` | `SZAlgoInterp.hpp` | Interpolation only |
| `ALGO_NOPRED` | `SZAlgoNopred.hpp` | Quantize-only baseline |
| `ALGO_LOSSLESS` | `SZDispatcher.hpp` | Zstd passthrough |
| `ALGO_BIOMD`, `ALGO_BIOMDXTC` | `SZAlgoBioMD.hpp` | Molecular dynamics |
| `ALGO_SVD` | `SZAlgoSVD.hpp` | Tucker/SVD core (3D+) |
| `ALGO_ZFP` | `SZAlgoZFP.hpp` | Bundled ZFP block transform |
| `ALGO_SPERR` | `SZAlgoSPERR.hpp` | Bundled SPERR wavelet+SPECK (3D FP only) |
| `ALGO_MGARD` | `SZAlgoMGARD.hpp` | Bundled MGARD multigrid (1D/2D/3D FP) |

Dispatch happens in `include/SZ3/api/impl/SZDispatcher.hpp` (`SZ_compress_dispatcher` and `SZ_decompress_dispatcher`).

## Modules

- **Decomposition** (`include/SZ3/decomposition/`): `BlockwiseDecomposition`, `InterpolationDecomposition`, `NoPredictionDecomposition`, `SVDDecomposition`, `ZFPDecomposition`, `SPERRDecomposition`, `MGARDDecomposition`, `SZBioMDDecomposition`, `SZBioMDXtcDecomposition`, `TimeSeriesDecomposition`.
- **Quantizer** (`include/SZ3/quantizer/`): `LinearQuantizer<T>` (default), `ScalarQuantizer<T,To>` (SPERR — asymmetric reconstruction), `FixedPointQuantizer<T>` (`ldexp` fixed-point, int64_t out, requires `calibrate(max_abs)`), `QuadraticLevelQuantizer<T>` (quadratic-LUT, opt-in via `sz_dev.hpp`).
- **Encoder** (`include/SZ3/encoder/`): `HuffmanEncoder<T>` (default), `ArithmeticEncoder<T>`, `BypassEncoder<T>`, `RunlengthEncoder<T>` (value-RLE), `BitplaneEncoder<T>` (MSB→LSB packed bit-planes), `BitplaneRLEEncoder<T>` (per-plane RLE w/ raw fallback), `BitshuffleEncoder<T>`, `SPERREncoder<T,N>`, `XtcBasedEncoder`, `ZFPEncoder<T>`.
- **Lossless** (`include/SZ3/lossless/`): `Lossless_zstd` (default), `Lossless_bypass`.

## Build & CLI

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target sz3 -j$(nproc)
```

CLI binary: `build/tools/sz3/sz3`. **Important — sz3 CLI takes dimensions in fastest-first order**: `-3 nx ny nz` means data is laid out as `data[nz][ny][nx]` in C row-major (so `nz` is outermost/slowest, `nx` is innermost/fastest). This is the **reverse** of how a C declaration reads. See the `fz-bench-multibound` skill for full CLI reference.

In-tree smoke dataset: `tools/sz3/testfloat_8_8_128.dat` (8×8×128 float32; pass as `-3 128 8 8`).

## Sibling skills

- `fz-compose-pipeline` — wire existing modules into a new composition
- `fz-add-module` — add a brand-new Decomposition / Encoder / Quantizer / Lossless module
- `fz-add-algorithm` — add a new `ALGO_*` enum + dispatch + wiring
- `fz-bench-multibound` — multi-bound benchmark recipe + sz3 CLI reference
