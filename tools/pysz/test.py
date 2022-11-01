import numpy as np
from pathlib import Path
from pysz import SZ
import platform

# prepare your data
HOME = str(Path.home())
data = np.fromfile(HOME + '/data/hurricane-100x500x500/Uf48.bin.dat', dtype=np.float32)
data = np.reshape(data, (100, 500, 500))

lib_extention = "so" if platform.system() == 'Linux' else "dylib"
# init SZ with the c dynamic library. Please change the path based on your system
sz3 = SZ("../../build/tools/sz3c/libSZ3c.{}".format(lib_extention))
# sz2 = SZ("../../../sz2/build/sz/libSZ.{}".format(lib_extention))
sz = sz3

# compress, both input and output data are numpy array
data_cmpr, cmpr_ratio = sz.compress(data, 0, 1e-3, 0, 0)
print("compression ratio = {:5G}".format(cmpr_ratio))

# decompress, both input and output data are numpy array
data_dec = sz.decompress(data_cmpr, data.shape, data.dtype)

# verify
sz.verify(data, data_dec)
