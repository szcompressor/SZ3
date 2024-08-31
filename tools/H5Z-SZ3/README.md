
# H5Z-SZ3 Filter Integration

The H5Z-SZ3 filter integrates the SZ3 compression library with HDF5, providing an efficient way to compress and decompress data within HDF5 files.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
  - [HDF5 Executables](#hdf5-executables)
  - [Compression Methods](#compression-methods)
  - [Decompression](#decompression)
- [Integration in C/C++](#integration-in-cc)

## Installation

### Step 1: Build and Install SZ3
Compile and install SZ3 with the H5Z-SZ3 filter enabled:
```bash
add -DBUILD_H5Z_FILTER=true to cmake command to enable H5Z-SZ3 filter in SZ3
```

### Step 2: Configure Environment
Add the directory containing the H5Z-SZ3 library (`libhdf5sz3.so` or `libhdf5sz3.dylib`) to your `HDF5_PLUGIN_PATH` and `LD_LIBRARY_PATH` environment variables:
```bash
export HDF5_PLUGIN_PATH=${HDF5_PLUGIN_PATH}:<PATH_TO_YOUR_H5Z-SZ3_LIB>/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_HOME/lib:$<PATH_TO_YOUR_H5Z-SZ3_LIB>
```


## Use H5Z-SZ3 in HDF5 Executables

Note, if the HDF5 in your system contains both `h5repack-shared` and `h5repack`, you need to use `h5repack-shared` for filters, as outlined in [H5Z-ZFP Issue #137](https://github.com/LLNL/H5Z-ZFP/issues/137).

### Compression

**Method 1: Compression with default settings (see SZ3/utils/Config.hpp for defaults)**
```bash
h5repack-shared -f UD=32024,0 data.h5 data.sz3.h5
```

**Method 2: Compression with customized settings in configuration file**

Step 1, generate H5Z-SZ3 parameters from SZ3 configuration file:
```bash
print_h5repack_args -c sz3.config
```
Step 2, compress by h5repack-shared with generated parameters:
```bash
h5repack-shared -f UD=32024,0,11,4077060608,16843009,0,0,4054449152,1348619730,16826431,257,2147483904,2147483648,16777216 data.h5 data.sz3.h5
```

Alternatively, use the script to generate parameters and compress at once:
```bash
h5repack.sh ~/code/sz3/tools/sz3/sz3.config data.h5 data.sz3.h5
```

### Decompression
```bash
h5repack-shared -f NONE data.sz3.h5 data.sz3_decompressed.h5
```

## Use H5Z-SZ3 in code (C/C++)
See examples `sz3ToHDF5.cpp` and `dsz3FromHDF5.cpp` for how to use the H5Z-SZ3 filter in your C/C++ projects.
