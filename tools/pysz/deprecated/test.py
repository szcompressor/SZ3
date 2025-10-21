import numpy as np
from cpysz import SZ
import sys

# prepare your data in numpy array format
data = np.random.uniform(low=0, high=1, size=(100, 100))

# init SZ (both SZ2 and SZ3 are supported)
# Please change the path to the SZ dynamic library file in your system
lib_extention = {
    "darwin": "libSZ3c.dylib",
    "windows": "SZ3c.dll",
}.get(sys.platform, "libSZ3c.so")

sz = SZ("../../install/lib/{}".format(lib_extention))

# compress, both input and output data are numpy array
data_cmpr, cmpr_ratio = sz.compress(data, 0, 1e-3, 0, 0)
print("compression ratio = {:5G}".format(cmpr_ratio))

# decompress, both input and output data are numpy array
data_dec = sz.decompress(data_cmpr, data.shape, data.dtype)

# verify
sz.verify(data, data_dec)
