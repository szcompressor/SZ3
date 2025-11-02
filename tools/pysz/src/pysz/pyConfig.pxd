"""Cython declarations for SZ3 Config class."""

from libcpp cimport bool
from libc.stdint cimport uint8_t, uint32_t
from libcpp.vector cimport vector
from libcpp.string cimport string

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
        
        # Predictor settings
        bool lorenzo
        bool lorenzo2
        bool regression
        bool regression2
        bool openmp
        
        # Encoding settings
        int quantbinCnt
        int blockSize
        uint8_t predDim
        uint8_t dataType
        
        # Interpolation settings
        uint8_t interpAlgo
        uint8_t interpDirection
        int interpAnchorStride
        double interpAlpha
        double interpBeta


cdef class pyConfig:
    cdef Config conf
