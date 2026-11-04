"""reflux from Python (M10a): thermodynamics and flashes on NumPy arrays,
through reflux's C API, capi/reflux.h. Units are SI throughout: K, Pa, mol, J.

    import reflux
    db = reflux.Database()
    fl = db.fluid("PR", ["methane", "n-butane"])

So far it has component data and fluids; properties and flashes come in
M10a's later steps (docs/m10a-plan.md).
"""

import ctypes

from ._native import (ABI_VERSION, BadArgument, BadHandle, DataError, NotFound, RefluxError,
                      check, lib)
from .database import Database, data_root
from .fluid import Fluid

__all__ = ["ABI_VERSION", "BadArgument", "BadHandle", "DataError", "Database", "Fluid",
           "NotFound", "RefluxError", "abi_version", "collect", "data_root", "version"]


def abi_version():
    """The library's ABI version, which this package checked when it loaded it."""
    return lib.reflux_abi_version()


def version():
    """reflux's version, v.mod's, such as "0.1.0"."""
    buf = ctypes.create_string_buffer(64)
    check(lib.reflux_version(buf, len(buf)), None)
    return buf.value.decode("utf-8")


def collect():
    """Run a full garbage collection in the library, and return its heap's size in
    bytes: for watching memory over many calls."""
    heap = ctypes.c_int64()
    check(lib.reflux_collect(ctypes.byref(heap)), None)
    return heap.value


__version__ = version()
