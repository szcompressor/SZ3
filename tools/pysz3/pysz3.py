import ctypes
from ctypes.util import find_library
import numpy as np

sz = ctypes.cdll.LoadLibrary("../../build/tools/sz3c/libsz3c.dylib")
libc = ctypes.CDLL(ctypes.util.find_library('c'))
libc.free.argtypes = (ctypes.c_void_p,)


def verify(src_data, dec_data):
    """
    Compare the decompressed data with original data
    :param src_data: original data, numpy array
    :param dec_data: decompressed data, numpy array
    :return: max_diff, psnr, nrmse
    """
    data_range = np.max(src_data) - np.min(src_data)
    diff = src_data - dec_data
    max_diff = np.max(abs(diff))
    print("abs err={:.8G}".format(max_diff))
    mse = np.mean(diff ** 2)
    nrmse = np.sqrt(mse) / data_range
    psnr = 20 * np.log10(data_range) - 10 * np.log10(mse)
    return max_diff, psnr, nrmse


def decompress(data_cmpr, original_shape):
    """
    Decompress
    :param data_cmpr: compressed data, numpy array format, dtype should be np.int8
    :param original_shape: the shape of original data
    :return: decompressed data,numpy array format
    """
    sz.SZ_decompress.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_ubyte), ctypes.c_size_t,
                                 ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
                                 ctypes.c_size_t)
    sz.SZ_decompress.restype = ctypes.POINTER(ctypes.c_float)

    r5, r4, r3, r2, r1 = [0] * (5 - len(original_shape)) + list(original_shape)
    data_dec_c = sz.SZ_decompress(0, data_cmpr.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)), np.size(data_cmpr),
                                  r5, r4, r3, r2, r1)

    data_dec = np.array(data_dec_c[:np.prod(original_shape)]).reshape(original_shape)
    libc.free(data_dec_c)
    return data_dec


def compress(data, eb_mode, eb_abs, eb_rel, eb_pwr):
    """
    Compress
    :param data: original data, numpy array format
    :param eb_mode:# error bound mode, integer (0: ABS, 1:REL, 2:ABS_AND_REL, 3:ABS_OR_REL, 4:PSNR, 5:NORM, 10:PW_REL)
    :param eb_abs: optional, abs error bound, double
    :param eb_rel: optional, rel error bound, double
    :param eb_pwr: optional, pwr error bound, double
    :return: compressed data, numpy array format, dtype=np.int8
    """
    sz.SZ_compress_args.argtypes = (ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t),
                                    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                    ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
                                    ctypes.c_size_t)
    sz.SZ_compress_args.restype = ctypes.POINTER(ctypes.c_ubyte)

    cmpr_size = ctypes.c_size_t()
    r5, r4, r3, r2, r1 = [0] * (5 - len(data.shape)) + list(data.shape)
    data_cmpr_c = sz.SZ_compress_args(0, data.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), ctypes.byref(cmpr_size),
                                      eb_mode, eb_abs, eb_rel, eb_pwr,
                                      r5, r4, r3, r2, r1)
    print('CR={:.3}'.format(np.size(data) * 4.0 / cmpr_size.value))

    data_cmpr = np.array(data_cmpr_c[:cmpr_size.value], dtype=np.int8)
    libc.free(data_cmpr_c)
    return data_cmpr
