# distutils: language = c++
"""Python wrapper for SZ3 Config class."""

from pysz cimport pyConfig
cimport cython


cdef class pyConfig:
    """
    Configuration class for SZ3 compression.
    
    Enum classes: pyConfig.EB, pyConfig.ALGO, pyConfig.INTERP_ALGO
    
    Parameters: Dimension sizes as individual ints, tuple, or list
        pyConfig(data.shape), pyConfig(100, 200, 300), etc.
    """

    def __init__(self, *args):
        """Initialize config with optional dimensions."""
        self.conf = Config()
        if args:
            self.setDims(*args)
    
    def setDims(self, *args):
        """Set dimensions. Accepts tuple/list or individual integers."""
        cdef vector[size_t] dims
        
        # Handle both setDims(100, 200, 300) and setDims((100, 200, 300))
        if len(args) == 1 and hasattr(args[0], '__iter__'):
            dims_iter = args[0]
        else:
            dims_iter = args
        
        if not dims_iter:
            raise ValueError("At least one dimension required")
        
        for arg in dims_iter:
            if not isinstance(arg, int) or arg <= 0:
                raise ValueError(f"Dimension must be positive integer, got {arg}")
            dims.push_back(<size_t>arg)
        
        self.conf.setDims(dims.begin(), dims.end())

    def loadcfg(self, cfgpath: str):
        """Load configuration from INI file."""
        cdef string cfgpathStr = cfgpath.encode('utf-8')
        try:
            self.conf.loadcfg(cfgpathStr)
        except RuntimeError as e:
            raise RuntimeError(f"Failed to load '{cfgpath}': {e}")
    
    # Read-only properties
    @property
    def dims(self):
        """Get dimensions (read-only)."""
        return tuple(self.conf.dims)
    
    @property
    def num_elements(self):
        """Get total number of elements (read-only)."""
        return self.conf.num
    
    @property
    def ndim(self):
        """Get number of dimensions (read-only)."""
        return self.conf.N
    
    # Error bounds (double)
    @property
    def absErrorBound(self):
        return self.conf.absErrorBound
    
    @absErrorBound.setter
    def absErrorBound(self, double value):
        self.conf.absErrorBound = value
    
    @property
    def relErrorBound(self):
        return self.conf.relErrorBound
    
    @relErrorBound.setter
    def relErrorBound(self, double value):
        self.conf.relErrorBound = value
    
    @property
    def psnrErrorBound(self):
        return self.conf.psnrErrorBound
    
    @psnrErrorBound.setter
    def psnrErrorBound(self, double value):
        self.conf.psnrErrorBound = value
    
    @property
    def l2normErrorBound(self):
        return self.conf.l2normErrorBound
    
    @l2normErrorBound.setter
    def l2normErrorBound(self, double value):
        self.conf.l2normErrorBound = value
    
    # Enum properties (accept enum or int)
    @property
    def errorBoundMode(self):
        """Get/set error bound mode. Use pyConfig.EB enum."""
        return self.conf.errorBoundMode
    
    @errorBoundMode.setter
    def errorBoundMode(self, value):
        if hasattr(value, 'value'):  # Is an enum
            self.conf.errorBoundMode = value.value
        else:
            raise TypeError(f"Use pyConfig.EB enum, got {type(value).__name__}")
    
    @property
    def cmprAlgo(self):
        """Get/set compression algorithm. Use pyConfig.ALGO enum."""
        return self.conf.cmprAlgo
    
    @cmprAlgo.setter
    def cmprAlgo(self, value):
        if hasattr(value, 'value'):  # Is an enum
            self.conf.cmprAlgo = value.value
        else:
            raise TypeError(f"Use pyConfig.ALGO enum, got {type(value).__name__}")
    
    @property
    def interpAlgo(self):
        """Get/set interpolation algorithm. Use pyConfig.INTERP_ALGO enum."""
        return self.conf.interpAlgo
    
    @interpAlgo.setter
    def interpAlgo(self, value):
        if hasattr(value, 'value'):  # Is an enum
            self.conf.interpAlgo = value.value
        else:
            raise TypeError(f"Use pyConfig.INTERP_ALGO enum, got {type(value).__name__}")
    
    # Boolean predictors
    @property
    def lorenzo(self):
        return self.conf.lorenzo
    
    @lorenzo.setter
    def lorenzo(self, bint value):
        self.conf.lorenzo = value
    
    @property
    def lorenzo2(self):
        return self.conf.lorenzo2
    
    @lorenzo2.setter
    def lorenzo2(self, bint value):
        self.conf.lorenzo2 = value
    
    @property
    def regression(self):
        return self.conf.regression
    
    @regression.setter
    def regression(self, bint value):
        self.conf.regression = value
    
    @property
    def regression2(self):
        return self.conf.regression2
    
    @regression2.setter
    def regression2(self, bint value):
        self.conf.regression2 = value
    
    @property
    def openmp(self):
        return self.conf.openmp
    
    @openmp.setter
    def openmp(self, bint value):
        self.conf.openmp = value
    
    # Integer settings
    @property
    def quantbinCnt(self):
        return self.conf.quantbinCnt
    
    @quantbinCnt.setter
    def quantbinCnt(self, int value):
        self.conf.quantbinCnt = value
    
    @property
    def blockSize(self):
        return self.conf.blockSize
    
    @blockSize.setter
    def blockSize(self, int value):
        self.conf.blockSize = value
    
    @property
    def interpDirection(self):
        return self.conf.interpDirection
    
    @interpDirection.setter
    def interpDirection(self, int value):
        self.conf.interpDirection = value
    
    @property
    def interpAnchorStride(self):
        return self.conf.interpAnchorStride
    
    @interpAnchorStride.setter
    def interpAnchorStride(self, int value):
        self.conf.interpAnchorStride = value
    
    # Double interpolation settings
    @property
    def interpAlpha(self):
        return self.conf.interpAlpha
    
    @interpAlpha.setter
    def interpAlpha(self, double value):
        self.conf.interpAlpha = value
    
    @property
    def interpBeta(self):
        return self.conf.interpBeta
    
    @interpBeta.setter
    def interpBeta(self, double value):
        self.conf.interpBeta = value
    
    def __repr__(self):
        return f"pyConfig(dims={self.dims}, num_elements={self.num_elements})"
