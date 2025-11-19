"""Cython declarations for SZ3 compression/decompression API."""

from libcpp cimport bool
from libc.stdint cimport uint8_t, uint32_t, int32_t, int64_t
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stddef cimport size_t

cdef extern from "SZ3/utils/Config.hpp" namespace "SZ3":
    cdef cppclass Config:
        # Constructor
        Config() except +
        
        # Methods
        size_t setDims[Iter](Iter begin, Iter end)
        void loadcfg(const string &cfgpath) 
        void load_ini(const string &ini_content)
        string save_ini() const
        size_t save(unsigned char *&c) const
        void load(const unsigned char *&c)
        void print()
        size_t size_est() const

        # Basic configuration
        uint32_t sz3MagicNumber
        uint32_t sz3DataVer
        char N
        vector[size_t] dims
        size_t num
        
        # Compression settings
        uint8_t cmprAlgo
        uint8_t errorBoundMode
        double absErrorBound
        double relErrorBound
        double psnrErrorBound
        double l2normErrorBound
        bool openmp


cdef class szConfig:
    cdef Config conf


cdef extern from "SZ3/api/sz.hpp":
    # Compression with pre-allocated buffer
    size_t SZ_compress[T](const Config &conf,
                          const T *data, 
                          char *cmpData, 
                          size_t cmpCap) except +
    
    # Decompression with pre-allocated buffer
    # Note: Config is non-const because it gets loaded from compressed data
    void SZ_decompress[T](Config &conf,
                         const char *cmpData, 
                         size_t cmpSize, 
                         T *&decData) except +

