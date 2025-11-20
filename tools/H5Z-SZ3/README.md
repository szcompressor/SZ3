# H5Z-SZ3 Filter Integration

The H5Z-SZ3 filter integrates the SZ3 compression library with HDF5, providing an efficient way to compress and decompress data within HDF5 files.

## Table of Contents
- [Installation](#installation)
- [H5Z-SZ3 cd_values](#h5z-sz3-cd_values)
- [Usage](#usage)
  - [HDF5 Executables](#hdf5-executables)
  - [Python (h5py)](#python-h5py)
  - [C/C++](#cc)

## Installation

### Step 1: Build and Install SZ3
Compile and install SZ3 with the H5Z-SZ3 filter enabled:
```bash
Add -DBUILD_H5Z_FILTER=true to the CMake command to enable H5Z-SZ3 filter in SZ3
```

### Step 2: Configure Environment
Add the directory containing the H5Z-SZ3 library (`libhdf5sz3.so` or `libhdf5sz3.dylib`) to your `HDF5_PLUGIN_PATH` and `LD_LIBRARY_PATH` environment variables:
```bash
export HDF5_PLUGIN_PATH=${HDF5_PLUGIN_PATH}:<PATH_TO_YOUR_H5Z-SZ3_LIB>/
```

## H5Z-SZ3 cd_values
* HDF5 restricts the parameters that can be passed to filters through an integers array called `cd_values`.
* H5Z-SZ3 uses `cd_values` to pass the desired compression settings (e.g., algorithm, error bounds) to the compression process. It serializes the `Config` object to `cd_values` using the `save()` function.
* During decompression, H5Z-SZ3 does not rely on `cd_values`. Instead, it reads all the configuration directly from the compressed data.

## Usage

### HDF5 Executables

Note: if the HDF5 in your system contains both `h5repack-shared` and `h5repack`, you need to use `h5repack-shared` for filters, as outlined in [H5Z-ZFP Issue #137](https://github.com/LLNL/H5Z-ZFP/issues/137).

#### Compression

**Method 1: Compression with default settings (see SZ3/utils/Config.hpp for defaults)**
```bash
h5repack-shared -f UD=32024,0 data.h5 data.sz3.h5
```

**Method 2: Compression with customized settings in configuration file**

Step 1: Generate H5Z-SZ3 cd_values from SZ3 configuration file:
```bash
cdvalueHelper -c sz3.config
```
Step 2: Compress by h5repack-shared with generated cd_values:
```bash
h5repack-shared -f UD=[put cd_values here] data.h5 data.sz3.h5
```

#### Decompression
```bash
h5repack-shared -f NONE data.sz3.h5 data.sz3_decompressed.h5
```

### Python (h5py)
To use H5Z-SZ3 with h5py in Python:

1. Ensure h5py is linked against the same HDF5 library used to build the SZ3 filter (typically system HDF5).
2. Set `HDF5_PLUGIN_PATH` to the directory containing `libhdf5sz3.so` or `libhdf5sz3.dylib` **before importing h5py**.
3. Use h5py's compression interface with SZ3 filter ID 32024 and appropriate compression options.

```python
import os
os.environ["HDF5_PLUGIN_PATH"] = "/path/to/h5z-sz3/lib"
import h5py
import numpy as np

from cdvalueHelper import SZ3
config = SZ3('ALGO_INTERP_LORENZO', absolute=1e-3)

with h5py.File('data.h5', 'w') as f:
    f.create_dataset('dataset', data=np.random.rand(100, 100), 
                     compression=32024, compression_opts=tuple(config.cd_values))
```

**Note on HDF5 Runtime Compatibility**: If h5py was installed with its bundled HDF5 (common on macOS), it may not load plugins built against system HDF5 due to library conflicts. To resolve:
- Rebuild h5py against system HDF5: `export HDF5_DIR=/path/to/system/hdf5 && pip install --force-reinstall --no-binary h5py h5py`
- Ensure HDF5 versions match between h5py and the plugin.

### C/C++
See examples `sz3ToHDF5.cpp` and `dsz3FromHDF5.cpp` for how to use the H5Z-SZ3 filter in your C/C++ projects.
