from ctypes import *
from pathlib import Path
import os
import ctypes
import numpy as np
import random

HOME = str(Path.home())

sotest = cdll.LoadLibrary("../../build/tools/sz3c/libsz3c.dylib")

# sotest.SZ_fast_compress_args.argtypes = [c_float, c_float]
sotest.SZ_compress_args.restype = c_char_p

outSize = c_size_t()

# INPUT = c_float * 1000
# data = INPUT()
# for i in range(1000):
#     data[i] = random.random()
#     # data[i] = 1
data = np.fromfile(HOME + '/data/hurricane-100x500x500/Uf48.bin.dat',
                   dtype=np.float32)
r3 = c_size_t(100)
r2 = c_size_t(500)
r1 = c_size_t(500)
data_p = data.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
addr = sotest.SZ_compress_args(0, data_p, byref(outSize),
                               0, c_double(1e-6), c_double(0), c_double(0),
                               c_size_t(0), c_size_t(0), r3, r2, r1)

# print(type(addr))
# print(addr)

print('CR={:.3}'.format(r3.value * r2.value * r2.value * 4.0 / outSize.value))
