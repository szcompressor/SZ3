"""Cython declarations for SZ3 compression/decompression API."""

from pysz cimport pyConfig
from libc.stdint cimport int32_t, int64_t
from libc.stddef cimport size_t

cdef extern from "SZ3/api/sz.hpp":
    # Compression with pre-allocated buffer
    size_t SZ_compress[T](const pyConfig.Config &conf, 
                          const T *data, 
                          char *cmpData, 
                          size_t cmpCap) except +
    
    # Decompression with pre-allocated buffer
    # Note: Config is non-const because it gets loaded from compressed data
    void SZ_decompress[T](pyConfig.Config &conf, 
                         const char *cmpData, 
                         size_t cmpSize, 
                         T *&decData) except +

