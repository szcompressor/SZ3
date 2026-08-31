# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What is fz

fz is a **modular toolkit for error-bounded scientific data compression**. It vendors SOTA compressors — SZ3 (interpolation/Lorenzo), ZFP (block transform), SPERR (wavelet/SPECK), MGARD (multigrid) — into one C++17 header-only library so developers can **compose new pipelines from a shared pool of modules** rather than picking one monolithic compressor.

The default API surface (the `sz3` CLI, `SZ_compress`/`SZ_decompress` entry points) selects an algorithm via `Config::cmprAlgo`. Each algorithm is a *composition* of modular stages; the modules themselves are reusable across compositions.

## The 4-stage pipeline

Every composition built on `SZGenericCompressor` follows this shape:

```
Input  →  Decomposition  →  (Quantizer)  →  Encoder  →  Lossless  →  Output
```

| Stage | Interface | Header | What it does |
|---|---|---|---|
| Decomposition | `concepts::DecompositionInterface<Ti, To, N>` | `include/SZ3/decomposition/Decomposition.hpp` | Lossy transform from `Ti` (float) to `To` (typically int bins). Often wraps + drives an internal Quantizer. |
| Quantizer | `concepts::QuantizerInterface<Ti, To>` | `include/SZ3/quantizer/Quantizer.hpp` | Maps a real value to an integer bin (and back). Used by most decompositions internally. |
| Encoder | `concepts::EncoderInterface<T>` | `include/SZ3/encoder/Encoder.hpp` | Bin sequence → byte stream (entropy coding). |
| Lossless | `concepts::LosslessInterface` | `include/SZ3/lossless/Lossless.hpp` | Final byte-stream compression (e.g., Zstd). |

Compose with `make_compressor_sz_generic<T, N>(decomp, encoder, lossless)` from `include/SZ3/compressor/SZGenericCompressor.hpp:151`. Output-range contract: `decomposition.get_out_range().first` MUST be `0` (the compressor enforces this).

## Available algorithms (`enum ALGO`)

Defined in `include/SZ3/utils/Config.hpp`. Each has a `SZAlgo*.hpp` wiring file in `include/SZ3/api/impl/`:

| Algo | Wiring | Pipeline |
|---|---|---|
| `ALGO_LORENZO_REG` | `SZAlgoLorenzoReg.hpp` | Blockwise Lorenzo + regression predictor |
| `ALGO_INTERP_LORENZO` *(default)* | `SZAlgoInterp.hpp` | Auto-tuned interpolation + Lorenzo |
| `ALGO_INTERP` | `SZAlgoInterp.hpp` | Interpolation only |
| `ALGO_NOPRED` | `SZAlgoNopred.hpp` | Quantize-only baseline |
| `ALGO_LOSSLESS` | `SZDispatcher.hpp` | Zstd passthrough |
| `ALGO_BIOMD`, `ALGO_BIOMDXTC` | `SZAlgoBioMD.hpp` | Molecular-dynamics specializations |
| `ALGO_SVD` | `SZAlgoSVD.hpp` | Tucker / SVD core decomposition (3D+) |
| `ALGO_ZFP` | `SZAlgoZFP.hpp` | Bundled ZFP block transform (FP-only) |
| `ALGO_SPERR` | `SZAlgoSPERR.hpp` | Bundled SPERR wavelet + SPECK (3D FP) |
| `ALGO_MGARD` | `SZAlgoMGARD.hpp` | Bundled MGARD multigrid + LinearQuantizer + HuffmanEncoder + Zstd (1D/2D/3D FP) |

## Modules (compose your own pipeline)

**Decomposition** (`include/SZ3/decomposition/`): `BlockwiseDecomposition`, `InterpolationDecomposition`, `NoPredictionDecomposition`, `SVDDecomposition`, `ZFPDecomposition`, `SPERRDecomposition` / `SPERRFusedDecomposition`, `MGARDDecomposition` (alias of `MultiLevelDecomposition`) / `MGARDFusedDecomposition`, `MultiLevelDecomposition`, `PaSTRIDecomposition`, `SZBioMDDecomposition`, `SZBioMDXtcDecomposition`, `TimeSeriesDecomposition`.

**Quantizer** (`include/SZ3/quantizer/`): `LinearQuantizer<T>` (linear delta, default), `ScalarQuantizer<T,To>` (delta with reconstruction tweaks, used by SPERR), `FixedPointQuantizer<T>` (`ldexp` fixed-point, int64_t out), `LevelQuantizer<T>` (non-uniform LUT, quadratic or log level curve chosen at construction, opt-in via `sz_dev.hpp`), `TimeIntQuantizer<T>` (time-series specialization).

**Encoder** (`include/SZ3/encoder/`): `HuffmanEncoder<T>` (default), `ArithmeticEncoder<T>`, `BypassEncoder<T>`, `RunlengthEncoder<T>` (value-RLE), `BitplaneEncoder<T>` (MSB→LSB packed bit-planes), `BitplaneRLEEncoder<T>` (per-plane RLE w/ raw fallback), `BitshuffleEncoder<T>`, `SPERREncoder<T,N>` (SPECK bitstream), `XtcBasedEncoder` (BioMD), `ZFPEncoder<T>`.

**Lossless** (`include/SZ3/lossless/`): `Lossless_zstd` (default), `Lossless_bypass`.

## Build

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --target sz3 -j$(nproc)
```

| Option | Default | Notes |
|---|---|---|
| `BUILD_SZ3_BINARY` | ON | `tools/sz3/sz3` CLI |
| `BUILD_TESTING` | OFF | GTest unit tests under `tools/test/modules/` |
| `BUILD_H5Z_FILTER` | OFF | HDF5 filter plugin |
| `BUILD_MDZ` | OFF | Molecular-dynamics CLI |
| `SZ3_DEBUG_TIMINGS` | OFF | Print per-stage timings |
| `SZ3_USE_BUNDLED_ZSTD` | OFF (ON for MSVC) | Bundled Zstd vs system |

**Use Release for any timing measurements** — Debug builds are 5-10× slower and unrepresentative.

System deps on Linux: `cmake`, `g++`, `libzstd-dev`, `libeigen3-dev`, `pkg-config`. For optional components: `libhdf5-dev` (H5Z filter), `libgsl-dev` (extra stats).

## Testing

- **Unit tests**: `cmake -DBUILD_TESTING=ON ..; make; ctest` — one executable per `.cpp` under `tools/test/modules/`.
- **Integration tests**: `tools/test/integration/` — Python scripts that exercise the `sz3` CLI on reference data.
- **Smoke dataset** committed in-tree: `tools/sz3/testfloat_8_8_128.dat` (8×8×128 float32).
- **CLI**: `tools/sz3/sz3` — sample invocation:
  ```
  sz3 -f -i in.f32 -z out.sz -3 nx ny nz -M ABS 1e-3
  sz3 -f -s out.sz -x out.dec.f32 -i in.f32 -3 nx ny nz -a   # -a prints CR + error metrics
  ```
- For end-to-end timing, build with `-DCMAKE_BUILD_TYPE=Release` (Debug is 5–10× slower).

## Directory layout

```
include/SZ3/
  api/                    # public entry points (sz.hpp), dispatch (impl/SZDispatcher.hpp), per-algo wirings (impl/SZAlgo*.hpp)
  compressor/             # SZGenericCompressor and other top-level compositions
  decomposition/          # all Decomposition modules
  quantizer/              # all Quantizer modules
  encoder/                # all Encoder modules
  lossless/               # all Lossless modules
  predictor/              # building blocks for blockwise predictor decompositions
  preprocessor/           # Transpose, Wavelet, ...
  utils/                  # Config, iterators, statistics, sperr utility types
  utils/thirdparty/       # bundled SOTA: sperr/, mgard/, zfp/, ska_hash/
tools/
  sz3/                    # main CLI (sz3.cpp, sz3.config, testfloat_8_8_128.dat)
  sz3/demo/               # standalone demos (zfp_demo, sz3_customized_demo)
  sz3c/                   # C API wrapper (SZ2-compatible)
  pysz/                   # Python (Cython) bindings
  H5Z-SZ3/                # HDF5 filter
  test/modules/           # GTest unit tests (per-module)
  test/integration/       # Python integration tests (sz3 CLI vs reference)
```

## Conventions

- **Header-only library** (`include/SZ3/`). New modules are `.hpp` files only, no `.cpp`. The `tools/` tree builds executables.
- **Namespace**: `SZ3` for everything in `include/SZ3/`. Bundled third-party namespaces are nested under `SZ3` (e.g., `SZ3::SPERR`, `SZ3::MGARD`, `SZ3::MGARDX`).
- **Save/load contract**: pointer-based streaming. `void save(uchar*& c)` advances `c`. `void load(const uchar*& c, size_t& remaining_length)` advances `c` and decrements `remaining_length`. Never read/write past `remaining_length`.
- **Output-range contract**: `decomposition.get_out_range().first == 0` is required by `SZGenericCompressor`. Quantizers that produce signed bins must shift to non-negative (see `LinearQuantizer::quantize_and_overwrite` for the pattern).
- **Code style**: clang-format (Google base, 120 col, 4-space indent). Run `clang-format -i <file>` before committing.
- **C++17** minimum (`cxx_std_17`). MSVC requires `/bigobj` (set automatically).

## Adding a new module or algorithm

The `docs/claude-skills/` directory ships installable Claude Code skills that walk through these workflows:
- `fz-overview` — pipeline + module map
- `fz-compose-pipeline` — wire existing modules into a new composition
- `fz-add-module` — add a brand-new Decomposition / Quantizer / Encoder / Lossless
- `fz-add-algorithm` — add a brand-new `ALGO_*` enum + dispatch + wiring
- `fz-bench-multibound` — multi-bound benchmark recipe

See `docs/claude-skills/README.md` for install instructions (symlink or copy into `.claude/skills/`).

Concise checklist for adding a new ALGO:
1. Add to `enum ALGO` and `ALGO_MAP` in `include/SZ3/utils/Config.hpp` (append at end to preserve serialized values)
2. Create `include/SZ3/api/impl/SZAlgoXxx.hpp` with templated `SZ_compress_Xxx` / `SZ_decompress_Xxx` following an existing wiring as template (e.g., `SZAlgoMGARD.hpp` or `SZAlgoNopred.hpp`)
3. Add `#include` and dispatch branches in `include/SZ3/api/impl/SZDispatcher.hpp` (compress + decompress)
4. Build, smoke-test on `tools/sz3/testfloat_8_8_128.dat`, then run a multi-bound on a real dataset
5. Optionally add a verify script under `tmp/`

## Bundled third-party (under `include/SZ3/utils/thirdparty/`)

| Source | Vendored from | Wraps to |
|---|---|---|
| `sperr/` | SPERR project | `SPERRFusedDecomposition`, `SPERREncoder` |
| `mgard/` | MGARDx (lightweight portable) | `MGARDFusedDecomposition` |
| `zfp/` | LLNL ZFP | `ZFPDecomposition`, `ZFPEncoder` |
| `ska_hash/` | skarupke flat_hash_map | shared utility |

## Compressed format

Little-endian: `[Header 16B: magic(4) + version(4) + compressed_size(8)] [Payload] [Config trailer]`. Magic: `0xF342F310`.

## Platform notes

- Targets Linux, macOS, Windows (MSVC + MinGW)
- `ALGO_SPERR` and `ALGO_MGARD` are disabled on MinGW (`#if !defined(__MINGW32__)` guards in `SZDispatcher.hpp`). The reason is link-time, not compile-time: libstdc++ autolink fails on the `libhdf5sz3.dll` link. The modules compile fine on MinGW and the unit tests build them there.
- macOS Mach-O binaries don't run in Linux containers — verify scripts auto-fall-back from `build-release/` to `build/` if the release binary fails `--help`
