"""reflux from Python (M10a): thermodynamics and flashes on NumPy arrays,
through reflux's C API, capi/reflux.h. Units are SI throughout: K, Pa, mol, J.

    import reflux
    db = reflux.Database()
    fl = db.fluid("PR", ["methane", "n-butane"])
    pr = fl.props(T=300.0, P=[1e5, 10e5, 50e5], z=[0.3, 0.7])
    pr.h, pr.status

It has component data, fluids, properties, vapour pressures, stability
and every flash: PT, PH, PS, TV, UV and the vapour-fraction flashes, with
bubble and dew points (docs/m10a-plan.md).
"""

import ctypes

from ._native import (ABI_VERSION, BadArgument, BadHandle, DataError, NotFound, RefluxError,
                      check, lib)
from .database import Database, data_root
from .fluid import Fluid
from .points import (FlashResult, Phase, PhaseView, Properties, Saturation, Stability,
                     Status)

__all__ = ["ABI_VERSION", "BadArgument", "BadHandle", "DataError", "Database", "FlashResult",
           "Fluid", "NotFound", "Phase", "PhaseView", "Properties", "RefluxError", "Saturation",
           "Stability", "Status", "abi_version", "collect", "data_root", "version"]


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
