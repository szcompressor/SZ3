# pysz

Python bindings for SZ3 - Error-bounded lossy compression for scientific data.

## Overview

pysz provides a clean Python interface to SZ3, a fast error-bounded lossy compressor for scientific data. Built with Cython 3.0+ for high performance.

## Installation

### From PyPI (recommended)

**Pre-built binary wheels** are available for most platforms (Linux, macOS, Windows) and Python versions (3.8-3.13):

```bash
pip install pysz
```

This is the easiest method - no build tools required! The binary wheels include everything you need.

### Building from source

If pre-built wheels aren't available for your platform, or if you want to build from source:

**You need to first install the following tools:**
- **CMake ≥ 3.13**
- **C++ compiler** (C++17-compatible: g++, clang++, MSVC)
- **Git**
- **Python development headers** (python3-dev or python3-devel)

**Then build pysz from pip:**
```bash
pip install pysz
```
**Or build pysz from source:**
```bash
git clone https://github.com/szcompressor/SZ3.git
cd SZ3/tools/pysz
pip install -e .
```

**What happens during source installation:**
1. SZ3 is automatically downloaded from GitHub
2. SZ3 is built with CMake (zstd is bundled)
3. Python bindings are compiled against the built SZ3
4. Everything is packaged together


## Quick Start

```python
import numpy as np
from pysz import sz, pyConfig

# Create test data
data = np.random.rand(8, 8, 128).astype(np.float32)

# Create config
config = pyConfig(data.shape)
config.errorBoundMode = pyConfig.EB.ABS
config.absErrorBound = 1e-3
config.cmprAlgo = pyConfig.ALGO.INTERP_LORENZO

# Compress
compressed, ratio = sz.compress(data, config)
print(f"Compression ratio: {ratio:.2f}x")

# Decompress
decompressed = sz.decompress(compressed, np.float32, data.shape, config)

# Verify
max_err, psnr, nrmse = sz.verify(data, decompressed)
print(f"Max error: {max_err:.2e}, PSNR: {psnr:.2f} dB, NRMSE: {nrmse:.2e}")
```

## API Reference

### pyConfig Class

Configuration object for SZ3 compression. It mirrors the C++ `Config` class.

```python
# Flexible initialization - all work the same:
config = pyConfig(data.shape)       # Pass shape tuple directly (recommended)
config = pyConfig((100, 200, 300))  # Pass tuple
config = pyConfig([100, 200, 300])  # Pass list
config = pyConfig(100, 200, 300)    # Pass individual dimensions
```

### sz.compress()

```python
sz.compress(data, config) -> (compressed, ratio)
```

Compress a NumPy array. Dimensions are automatically inferred from data shape.

**Parameters:**
- `data` (ndarray): NumPy array (float32, float64, int32, or int64)
- `config` (pyConfig or str): Config object or path to config file

**Returns:**
- `compressed` (ndarray): Compressed data as uint8 array
- `ratio` (float): Compression ratio (original_size / compressed_size)

### sz.decompress()

```python
sz.decompress(compressed, dtype, shape, config) -> data
```

Decompress data back to NumPy array.

**Parameters:**
- `compressed` (ndarray): Compressed uint8 array from `compress()`
- `dtype` (type): NumPy dtype (np.float32, np.float64, np.int32, or np.int64)
- `shape` (tuple): Shape of the original data
- `config` (pyConfig or str): Config object or path to config file

**Returns:**
- `data` (ndarray): Decompressed data with the specified shape

### sz.verify()

```python
sz.verify(src_data, dec_data) -> (max_diff, psnr, nrmse)
```

Compare decompressed data with original data and calculate quality metrics.

**Parameters:**
- `src_data` (ndarray): Original data before compression
- `dec_data` (ndarray): Decompressed data to verify

**Returns:**
- `max_diff` (float): Maximum absolute difference
- `psnr` (float): Peak Signal-to-Noise Ratio in dB
- `nrmse` (float): Normalized Root Mean Square Error

**Example:**
```python
>>> max_err, psnr, nrmse = sz.verify(data, decompressed)
>>> print(f"Max error: {max_err:.2e}, PSNR: {psnr:.2f} dB")
```

## Usage Examples

### Basic Usage

```python
import numpy as np
from pysz import sz, pyConfig

# Create data
data = np.random.randn(100, 200, 300).astype(np.float32)

# Create config - use enums
config = pyConfig(data.shape)
config.errorBoundMode = pyConfig.EB.ABS
config.absErrorBound = 1e-3

# Compress
compressed, ratio = sz.compress(data, config)
print(f"Compressed {data.nbytes} → {compressed.size} bytes ({ratio:.2f}x)")

# Decompress
decompressed = sz.decompress(compressed, np.float32, data.shape, config)

# Verify quality
max_err, psnr, nrmse = sz.verify(data, decompressed)
print(f"Max error: {max_err:.2e}, PSNR: {psnr:.2f} dB, NRMSE: {nrmse:.2e}")
```

### Loading from Config File (Optional)

```python
import numpy as np
from pysz import sz, pyConfig

data = np.fromfile('testdata.dat', dtype=np.float32).reshape(8, 8, 128)

# Load config from file
config = pyConfig(data.shape)
config.loadcfg('sz3.config')
compressed, ratio = sz.compress(data, config)
decompressed = sz.decompress(compressed, np.float32, data.shape, config)

# Or pass file path directly
compressed, ratio = sz.compress(data, 'sz3.config')
decompressed = sz.decompress(compressed, np.float32, data.shape, 'sz3.config')
```

### Relative Error Bound

```python
config = pyConfig(data.shape)
config.errorBoundMode = pyConfig.EB.REL
config.relErrorBound = 1e-4  # 0.01% relative error
compressed, ratio = sz.compress(data, config)
```


### Different Algorithms

```python
# Default: INTERP_LORENZO (best quality)
config = pyConfig(data.shape)
config.errorBoundMode = pyConfig.EB.ABS
config.absErrorBound = 1e-3

# Or try other algorithms
config.cmprAlgo = pyConfig.ALGO.INTERP        # Interpolation only
config.cmprAlgo = pyConfig.ALGO.LORENZO_REG   # Lorenzo/regression
config.cmprAlgo = pyConfig.ALGO.LOSSLESS      # Lossless only
```

### Double Precision

```python
data_double = np.random.randn(50, 50, 50).astype(np.float64)
config = pyConfig(data_double.shape)
config.errorBoundMode = pyConfig.EB.ABS
config.absErrorBound = 1e-6
compressed, ratio = sz.compress(data_double, config)
decompressed = sz.decompress(compressed, np.float64, data_double.shape, config)
```

### Save/Load Compressed Data

```python
# Compress and save
compressed, ratio = sz.compress(data, config)
compressed.tofile('data.sz')

# Later: load and decompress
compressed = np.fromfile('data.sz', dtype=np.uint8)
decompressed = sz.decompress(compressed, np.float32, (8, 8, 128), config)
```

## Troubleshooting

### Import Error

```
ImportError: cannot import name 'sz' from 'pysz'
```

**Solution:**
```bash
cd tools/pysz
pip install -e .
```

### Library Not Found (Linux)

```
OSError: libzstd.so: cannot open shared object file
```

**Solution:**
```bash
export LD_LIBRARY_PATH="../../build/tools/zstd:$LD_LIBRARY_PATH"
```


## Links

- **Repository:** https://github.com/szcompressor/SZ3
- **Issues:** https://github.com/szcompressor/SZ3/issues

## License

See `../../copyright-and-BSD-license.txt`
