"""The C API as ctypes sees it: each function capi/reflux.h declares, the
library that exports them, and the codes they return.

Integers at the boundary are fixed widths, ctypes.c_int32 and c_int64, never
ctypes.c_int: V's `int` is 64-bit, C's 32. Strings come back through buffers
the caller allocates, so nothing V allocated ever reaches Python.
"""

import ctypes
import os
import pathlib
import sys

ABI_VERSION = 1

# A call's result (reflux.h).
OK = 0
BAD_ARGUMENT = -1
BAD_HANDLE = -2

# Every function in reflux.h: its result's type and its arguments' types.
FUNCTIONS = {
    "reflux_abi_version": (ctypes.c_int32, []),
    "reflux_version": (ctypes.c_int32, [ctypes.c_char_p, ctypes.c_int32]),
    "reflux_free": (ctypes.c_int32, [ctypes.c_int64, ctypes.c_char_p, ctypes.c_int32]),
    "reflux_collect": (ctypes.c_int32, [ctypes.POINTER(ctypes.c_int64)]),
}

# Room for a reason. A longer one is cut at a character boundary.
MESSAGE_SIZE = 1024


class RefluxError(Exception):
    """A failure reflux reports: its code, and its reason."""

    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


class BadArgument(RefluxError, ValueError):
    """An argument the library can't use."""


class BadHandle(RefluxError):
    """An id that isn't a live handle: never given out, or freed already."""


_ERRORS = {BAD_ARGUMENT: BadArgument, BAD_HANDLE: BadHandle}


def check(code, buf):
    """Raise the error a call's code stands for, with the reason in buf."""
    if code == OK:
        return
    message = buf.value.decode("utf-8", "replace") if buf is not None else ""
    raise _ERRORS.get(code, RefluxError)(code, message or f"reflux failed with code {code}")


def library_path():
    """Where the library is: $REFLUX_LIB, or else _lib/ beside this file."""
    env = os.environ.get("REFLUX_LIB")
    if env:
        return pathlib.Path(env)
    name = "reflux.dll" if sys.platform == "win32" else "reflux.so"
    return pathlib.Path(__file__).parent / "_lib" / name


def load(path):
    """Load the library and declare its functions.

    The library starts V's runtime on its first call, from whichever thread
    makes it, and that start mustn't race (reflux.h). So the first call is
    made here, by the import, which Python's import lock serialises.
    """
    if not path.exists():
        raise ImportError(f"reflux's library isn't built: there's no {path}. "
                          "Build it with `v run tools/python.vsh build`.")
    lib = ctypes.CDLL(str(path))
    for name, (restype, argtypes) in FUNCTIONS.items():
        f = getattr(lib, name)
        f.restype = restype
        f.argtypes = argtypes
    abi = lib.reflux_abi_version()
    if abi != ABI_VERSION:
        raise ImportError(f"{path} has ABI version {abi}, and this package needs {ABI_VERSION}")
    return lib


lib = load(library_path())


def message_buffer():
    """A buffer for a call's reason."""
    return ctypes.create_string_buffer(MESSAGE_SIZE)


def free(handle):
    """Drop a handle: an object's own close calls it."""
    buf = message_buffer()
    check(lib.reflux_free(handle, buf, MESSAGE_SIZE), buf)
