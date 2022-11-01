import numpy as np
from pathlib import Path
from pysz3 import compress, decompress, verify

HOME = str(Path.home())
data = np.fromfile(HOME + '/data/hurricane-100x500x500/Uf48.bin.dat', dtype=np.float32)
data = np.reshape(data, (100, 500, 500))

data_cmpr = compress(data, 0, 1e-3, 0, 0)

data_dec = decompress(data_cmpr, data.shape)

verify(data, data_dec)
