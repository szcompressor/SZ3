---
name: fz-bench-multibound
description: Use when the user wants to benchmark a compressor across multiple error bounds, capture compression ratio / time / max-abs error, or compare two algorithms head-to-head. Triggers include "benchmark sz3", "multi-bound comparison", "compare algorithms across error bounds", "measure compression ratio at multiple eb".
---

# Multi-bound benchmark

Loop the `sz3` CLI over a list of error bounds (eb), capture metrics, render a table.

## Build for accurate timing

**Always use Release builds for timing.** Debug is 5–10× slower and meaningless for benchmarking.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target sz3 -j$(nproc)
```

## sz3 CLI cheat-sheet

```
sz3 -f -i <in.f32> -z <out.sz> -3 <nx> <ny> <nz> -M ABS <eb> [-c <cfg.ini>]   # compress
sz3 -f -s <out.sz> -x <out.dec.f32> -i <in.f32> -3 <nx> <ny> <nz> -a          # decompress + print -a metrics
```

Flags:

| Flag | Meaning |
|---|---|
| `-f` / `-d` / `-I 32|64` | data type (float / double / int32 / int64) |
| `-i <path>` | input data file (binary) |
| `-z <path>` | compressed output (used during compression) |
| `-s <path>` | compressed input (used during decompression) |
| `-x <path>` | decompressed output |
| `-1 nx` / `-2 nx ny` / `-3 nx ny nz` / `-4 nx ny nz np` | data dimensions — **see ordering note below** |
| `-M <mode> [bound]` | error bound mode (`ABS`, `REL`, `PSNR`, `NORM`, `ABS_AND_REL`, `ABS_OR_REL`) |
| `-A <eb>` / `-R <eb>` / `-S <psnr>` / `-N <l2err>` | bound value alternatives (mode-specific) |
| `-c <cfg.ini>` | INI config file (lets you set `CmprAlgo`, etc. — see below) |
| `-a` | print metrics: max-abs / max-rel / max-pw-rel error, PSNR, NRMSE, L2 error, CR |

### ⚠️ Dimension order is FASTEST-FIRST (reverse of C declaration)

`sz3 -3 NX NY NZ` describes a 3D buffer laid out in C row-major as `data[NZ][NY][NX]`.

- `NX` = innermost = **fastest-varying** dim (smallest stride)
- `NZ` = outermost = **slowest-varying** dim (largest stride)

So if you're translating from a C/Python declaration that reads slowest-first (e.g., NumPy `arr.shape == (256, 384, 384)` for `arr[256][384][384]`), **reverse the order** when passing to sz3:

| Source declaration | sz3 invocation |
|---|---|
| `float data[256][384][384]` (256 slowest, 384 fastest) | `-3 384 384 256` |
| NumPy `arr.shape == (Nz, Ny, Nx)` | `-3 Nx Ny Nz` |
| `nyx-512x512x512-velocity_x.dat` (cubic) | `-3 512 512 512` |

The same convention applies to `-2` (sz3 `-2 NX NY` = C `data[NY][NX]`) and `-4`. Mixing this up doesn't crash — it silently transposes the field, giving meaningless CR and inflated reconstruction error.

## Selecting the algorithm

Either pass `-c <cfg.ini>` with `CmprAlgo = ALGO_FOO`, or omit `-c` to use the default `ALGO_INTERP_LORENZO`.

A minimal config:
```ini
[GlobalSettings]
CmprAlgo = ALGO_MGARD
ErrorBoundMode = ABS
AbsErrorBound = 1e-3
```

Valid `CmprAlgo` values: see `enum ALGO` in `include/SZ3/utils/Config.hpp` and the `ALGO_MAP` strings.

## Multi-bound bash pattern

```bash
#!/usr/bin/env bash
set -uo pipefail

DATA=/path/to/your/dataset.f32
SZ3=$(pwd)/build/tools/sz3/sz3
WORK=$(mktemp -d /tmp/bench.XXXXXX)
ORIG=$(stat -c%s "$DATA")

# Dim order is FASTEST-FIRST. Adjust to your dataset.
NX=384 ; NY=384 ; NZ=256

printf "%-25s %-8s %-12s %-12s %-12s\n" "case" "CR" "compress(s)" "decompress(s)" "max_abs_err"
echo "------------------------------------------------------------------------------------"

for ALGO in ALGO_INTERP_LORENZO ALGO_MGARD; do
  for EB in 1e-2 1e-3 1e-4 1e-5; do
    CFG=$WORK/$ALGO.$EB.cfg
    printf '[GlobalSettings]\nCmprAlgo=%s\nErrorBoundMode=ABS\nAbsErrorBound=%s\n' \
      "$ALGO" "$EB" > "$CFG"

    CMP_LOG=$WORK/cmp.$ALGO.$EB.log
    DEC_LOG=$WORK/dec.$ALGO.$EB.log

    /usr/bin/time -p "$SZ3" -f -c "$CFG" -i "$DATA" \
      -z "$WORK/$ALGO.$EB.sz" -3 $NX $NY $NZ -M ABS $EB > "$CMP_LOG" 2>&1
    /usr/bin/time -p "$SZ3" -f -s "$WORK/$ALGO.$EB.sz" \
      -x "$WORK/$ALGO.$EB.out" -i "$DATA" -3 $NX $NY $NZ -a > "$DEC_LOG" 2>&1

    SZ=$(stat -c%s "$WORK/$ALGO.$EB.sz")
    CR=$(awk -v o=$ORIG -v c=$SZ 'BEGIN{printf "%.2f", o/c}')
    CMP=$(awk '/^real/{print $2}' "$CMP_LOG")
    DEC=$(awk '/^real/{print $2}' "$DEC_LOG")
    ERR=$(awk -F'= *' '/Max absolute error/{print $2; exit}' "$DEC_LOG")
    printf "%-25s %-8s %-12s %-12s %-12s\n" "$ALGO eb=$EB" "${CR}x" "$CMP" "$DEC" "$ERR"
  done
  echo
done
echo "Logs in $WORK"
```

## Common gotchas

- **Stale binaries on a host with mixed builds**. If you have both `build-release/` (e.g., from a different OS) and `build/`, the OS-mismatched one will pass `[ -x ]` but fail to exec. Set `SZ3=$(realpath build/tools/sz3/sz3)` explicitly, or probe with `"$SZ3" --help >/dev/null 2>&1` first.
- **`/usr/bin/time` not installed**. On minimal containers, `apt-get install -y time` (or fall back to bash's `time` builtin with `TIMEFORMAT="real %R"`).
- **`-a` field names are case-sensitive in awk**. The exact lines printed by sz3 are `Max absolute error = ...`, `Max relative error = ...`, `PSNR = ...`. Anchor your awk on the full prefix or you'll match the wrong line (e.g., `Max relative` matches `Max absolute` if you use `/Max/`).
- **Forgetting `-i <orig>` on decompress**. Without it, sz3 still decompresses but `-a` cannot compare against ground truth and prints zeros / NaNs.
