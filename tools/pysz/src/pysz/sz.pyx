# distutils: language = c++
"""Python interface for SZ3 compression library."""

from pysz cimport sz as c_sz
from pysz cimport pyConfig
cimport cython
from cython.operator cimport dereference
import numpy as np
cimport numpy as cnp
from libc.stdint cimport int8_t, int32_t, int64_t
from libc.stddef cimport size_t
from libc.string cimport memcpy
from typing import Tuple

# Initialize NumPy C API
cnp.import_array()


cdef class sz:
    """SZ3 compression/decompression with zero-copy NumPy API."""
    
    # Supported dtypes
    _SUPPORTED_DTYPES = frozenset([
        np.dtype(np.float32),
        np.dtype(np.float64),
        np.dtype(np.int32),
        np.dtype(np.int64),
    ])
    
    @staticmethod
    def compress(cnp.ndarray data, config) -> Tuple[cnp.ndarray, float]:
        """
        Compress NumPy array using SZ3.
        
        Parameters
        ----------
        data : numpy.ndarray
            Input data (float32, float64, int32, or int64)
        config : pyConfig or str
            Configuration object or path to config file
        
        Returns
        -------
        compressed : numpy.ndarray
            Compressed data as uint8 array
        ratio : float
            Compression ratio (original_size / compressed_size)
        """
        # Validate dtype
        if data.dtype not in sz._SUPPORTED_DTYPES:
            raise TypeError(
                f"Unsupported dtype: {data.dtype}. "
                f"Supported: float32, float64, int32, int64"
            )
        
        # Handle config parameter
        cdef pyConfig.pyConfig conf
        if isinstance(config, str):
            # Config file path provided - convert shape to tuple
            shape_tuple = tuple(<Py_ssize_t>data.shape[i] for i in range(data.ndim))
            conf = pyConfig.pyConfig(*shape_tuple)
            conf.loadcfg(config)
        elif isinstance(config, pyConfig.pyConfig):
            # Config object provided - set dimensions if needed
            conf = config
            if conf.num_elements != data.size:
                shape_tuple = tuple(<Py_ssize_t>data.shape[i] for i in range(data.ndim))
                conf.setDims(*shape_tuple)
        else:
            raise TypeError(f"config must be pyConfig or str, got {type(config)}")
        
        # Ensure C-contiguous for zero-copy
        if not data.flags['C_CONTIGUOUS']:
            data = np.ascontiguousarray(data)
        
        # Get data pointer (zero-copy input!)
        cdef void* data_ptr = <void*> cnp.PyArray_DATA(data)
        cdef size_t original_size = data.nbytes
        
        # Allocate buffer for compressed data (2x original size to be safe)
        cdef size_t buffer_size = <size_t>(original_size * 2)
        cdef cnp.ndarray[cnp.uint8_t, ndim=1] compressed = np.empty(buffer_size, dtype=np.uint8)
        cdef char* compressed_ptr = <char*> cnp.PyArray_DATA(compressed)
        cdef size_t compressed_size = 0
        
        # Compress into pre-allocated buffer
        if data.dtype == np.float32:
            compressed_size = c_sz.SZ_compress[float](
                dereference(&conf.conf),
                <float*>data_ptr,
                compressed_ptr,
                buffer_size
            )
        elif data.dtype == np.float64:
            compressed_size = c_sz.SZ_compress[double](
                dereference(&conf.conf),
                <double*>data_ptr,
                compressed_ptr,
                buffer_size
            )
        elif data.dtype == np.int32:
            compressed_size = c_sz.SZ_compress[int32_t](
                dereference(&conf.conf),
                <int32_t*>data_ptr,
                compressed_ptr,
                buffer_size
            )
        elif data.dtype == np.int64:
            compressed_size = c_sz.SZ_compress[int64_t](
                dereference(&conf.conf),
                <int64_t*>data_ptr,
                compressed_ptr,
                buffer_size
            )
        
        if compressed_size == 0:
            raise RuntimeError("Compression failed")
        
        # Resize to actual compressed size (no copy, just view)
        compressed = compressed[:compressed_size]
        
        ratio = original_size / float(compressed_size)
        return compressed, ratio
    
    @staticmethod
    def decompress(cnp.ndarray compressed, dtype, shape, config) -> cnp.ndarray:
        """
        Decompress SZ3-compressed data.
        
        Parameters
        ----------
        compressed : numpy.ndarray
            Compressed data (uint8 array)
        dtype : numpy.dtype or type
            Data type of original data
        shape : tuple
            Shape of the original data
        config : pyConfig or str
            Configuration object or path to config file
        
        Returns
        -------
        decompressed : numpy.ndarray
            Decompressed data with the specified shape
        """
        # Validate compressed data
        if compressed.dtype != np.uint8:
            raise TypeError(f"Compressed data must be uint8, got {compressed.dtype}")
        
        # Validate and normalize dtype
        dtype = np.dtype(dtype)
        if dtype not in sz._SUPPORTED_DTYPES:
            raise TypeError(
                f"Unsupported dtype: {dtype}. "
                f"Supported: float32, float64, int32, int64"
            )
        
        # Handle config parameter
        cdef pyConfig.pyConfig conf
        if isinstance(config, str):
            # Config file path provided
            conf = pyConfig.pyConfig(*shape)
            conf.loadcfg(config)
        elif isinstance(config, pyConfig.pyConfig):
            # Config object provided - set dimensions if needed
            conf = config
            expected_size = int(np.prod(shape))
            if conf.num_elements != expected_size:
                conf.setDims(*shape)
        else:
            raise TypeError(f"config must be pyConfig or str, got {type(config)}")
        
        # Ensure compressed data is contiguous
        if not compressed.flags['C_CONTIGUOUS']:
            compressed = np.ascontiguousarray(compressed)
        
        cdef char* compressed_ptr = <char*> cnp.PyArray_DATA(compressed)
        cdef size_t compressed_size = compressed.size
        cdef size_t num_elements = conf.num_elements
        
        # Pre-allocate NumPy array for decompressed data
        cdef cnp.ndarray result = np.empty(num_elements, dtype=dtype)
        
        # Get pointer to NumPy array data
        # Pass this to SZ_decompress - it will decompress directly into our buffer!
        cdef float* float_ptr
        cdef double* double_ptr
        cdef int32_t* int32_ptr
        cdef int64_t* int64_ptr
        
        if dtype == np.float32:
            float_ptr = <float*> cnp.PyArray_DATA(result)
            c_sz.SZ_decompress[float](
                conf.conf,
                compressed_ptr,
                compressed_size,
                float_ptr
            )
        elif dtype == np.float64:
            double_ptr = <double*> cnp.PyArray_DATA(result)
            c_sz.SZ_decompress[double](
                conf.conf,
                compressed_ptr,
                compressed_size,
                double_ptr
            )
        elif dtype == np.int32:
            int32_ptr = <int32_t*> cnp.PyArray_DATA(result)
            c_sz.SZ_decompress[int32_t](
                conf.conf,
                compressed_ptr,
                compressed_size,
                int32_ptr
            )
        elif dtype == np.int64:
            int64_ptr = <int64_t*> cnp.PyArray_DATA(result)
            c_sz.SZ_decompress[int64_t](
                conf.conf,
                compressed_ptr,
                compressed_size,
                int64_ptr
            )
        
        # Data is now in our NumPy array! (zero-copy decompression)
        # Reshape to original dimensions
        return result.reshape(shape)
    
    @staticmethod
    def verify(cnp.ndarray src_data, cnp.ndarray dec_data) -> Tuple[float, float, float]:
        """
        Compare decompressed data with original data.
        
        Parameters
        ----------
        src_data : numpy.ndarray
            Original data before compression
        dec_data : numpy.ndarray
            Decompressed data to verify
        
        Returns
        -------
        max_diff : float
            Maximum absolute difference
        psnr : float
            Peak Signal-to-Noise Ratio in dB
        nrmse : float
            Normalized Root Mean Square Error
        """
        # Check shapes match by comparing size and ndim
        if src_data.ndim != dec_data.ndim or src_data.size != dec_data.size:
            raise ValueError("Shape mismatch between src_data and dec_data")
        
        if src_data.dtype != dec_data.dtype:
            raise ValueError("Dtype mismatch between src_data and dec_data")
        
        # Calculate data range and difference
        cdef double data_range = np.max(src_data) - np.min(src_data)
        cdef cnp.ndarray diff = src_data - dec_data
        cdef double max_diff = np.max(np.abs(diff))
        
        # Calculate MSE and derived metrics
        cdef double mse = np.mean(diff ** 2)
        cdef double nrmse = np.sqrt(mse) / data_range if data_range > 0 else 0.0
        cdef double psnr = 20 * np.log10(data_range) - 10 * np.log10(mse) if mse > 0 else float('inf')
        
        return max_diff, psnr, nrmse
