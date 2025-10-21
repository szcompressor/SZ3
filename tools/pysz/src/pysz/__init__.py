"""Python bindings for SZ3 error-bounded lossy compression."""

from pysz.sz import sz
from pysz.pyConfig import pyConfig as _pyConfigBase
from pysz.pyConfigEnums import EB, ALGO, INTERP_ALGO


class pyConfig(_pyConfigBase):
    """SZ3 configuration with integrated enums."""
    EB = EB
    ALGO = ALGO
    INTERP_ALGO = INTERP_ALGO


__all__ = ["sz", "pyConfig"]
