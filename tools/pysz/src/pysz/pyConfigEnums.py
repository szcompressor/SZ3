"""Enum classes for SZ3 configuration matching C++ enums."""

from enum import IntEnum


class EB(IntEnum):
    """Error bound modes"""
    ABS = 0
    REL = 1
    PSNR = 2
    L2NORM = 3
    ABS_AND_REL = 4
    ABS_OR_REL = 5


class ALGO(IntEnum):
    """Compression algorithms"""
    LORENZO_REG = 0
    INTERP_LORENZO = 1
    INTERP = 2
    NOPRED = 3
    LOSSLESS = 4


class INTERP_ALGO(IntEnum):
    """Interpolation algorithms"""
    LINEAR = 0
    CUBIC = 1
