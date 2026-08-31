# Module acceptance

What a module must satisfy to enter the FZ module library. Every item is enforced by something
that runs, not by review: the command to run it is given with each one.

Each check exists because the corresponding defect has actually shipped in this tree. Passing
means a module is known not to repeat one of them — it is not a proof of correctness.

## 1. It inherits its group's interface

| Group | Interface | Header |
|---|---|---|
| Decomposition | `concepts::DecompositionInterface<Ti, To, N>` | `include/SZ3/decomposition/Decomposition.hpp` |
| Quantizer | `concepts::QuantizerInterface<Ti, To>` | `include/SZ3/quantizer/Quantizer.hpp` |
| Encoder | `concepts::EncoderInterface<T>` | `include/SZ3/encoder/Encoder.hpp` |
| Lossless | `concepts::LosslessInterface` | `include/SZ3/lossless/Lossless.hpp` |
| Preprocessor | `concepts::PreprocessorInterface<T, N>` | `include/SZ3/preprocessor/PreProcessor.hpp` |

A new module never changes an interface definition. Widening one breaks every out-of-tree
implementation of it: an override with the old signature stops overriding, and the derived class
silently becomes abstract.

## 2. Its header compiles on its own

```bash
python3 tools/test/check_headers.py --build-dir build
```

A header that only works because another header came first breaks the moment a translation unit
changes its include order, or a standard library stops providing something transitively. Run it
against both toolchains -- libstdc++ and libc++ provide different things transitively, so a clean
run under one proves nothing about the other. CI runs it on the Linux and macOS jobs.

```bash
python3 tools/test/check_headers.py --build-dir build --compiler g++
```

## 3. It passes its group contract

```bash
ctest -R SZ3_ModuleContract
```

`tools/test/modules/ModuleContract.hpp` holds the checks that apply to every module in a group;
a new module is added to `test_module_contract.cpp` in a few lines. What the contract covers:

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

## 4. It has its own tests

One `tools/test/modules/test_*.cpp` per module, covering whatever the group contract cannot know:
the module's error guarantee, its serialized format, its parameters, its documented edge cases.
CMake globs the directory, so a new file needs no build change.

## 5. It works in a composition

A module that round-trips in isolation can still fail inside a pipeline. Add a case to
`tools/test/modules/test_composition.cpp` wiring it through `make_compressor_sz_generic` with a
real encoder and lossless stage. Most of the defects listed below were found this way, not by the
per-module tests.

## 6. It is registered

- `include/SZ3/api/sz_dev.hpp` — the developer header carries every module.
- `docs/site/modules.json` — the metadata the module site reads.
- A `uid` distinct from every other module in its group, if the group serializes one.
  `ctest -R UidsAreDistinct` checks this by reading the first byte each `save()` writes.

## 7. Clean build, both toolchains

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build
clang-format -i <changed files>
```

Zero warnings. CI builds Linux, macOS, Windows and Windows MinGW, and compares Linux/macOS
digests; PRs to `master` or `fz` additionally run the SDRBench datasets.

## Defect classes these checks came from

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
| composition test | `TimeSeriesDecomposition` violated its bound 1.97x on the null-reference path |
| `size_est()` covers `save()` | six decompositions did not override `size_est()`, so the compressor sized its buffer from 0 |
| bound after `save()`/`load()` | `RegressionPredictor::load` debited `remaining_length` by the decoded rather than the encoded size |
