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
NOT_FOUND = -3
DATA_ERROR = -4
TOO_SMALL = -5

# Equations of state and α functions, as reflux.h numbers them.
EQUATIONS = {"VDW": 0, "RK": 1, "SRK": 2, "PR": 3}
ALPHAS = {"equation": 0, "mathias_copeman": 1}

# Volume roots, for reflux_props.
ROOTS = {"liquid": 0, "vapor": 1, "stable": 2}

# reflux_props' properties, in the order it writes them (REFLUX_PROP_*).
PROPS = ("z", "v", "h", "s", "g", "u", "a", "cp", "cv", "w", "mu_jt")

# Flashes, for reflux_flash: what a point's two values are (REFLUX_FLASH_PT...).
FLASHES = {"pt": 0, "ph": 1, "ps": 2, "tv": 3, "uv": 4, "tvf": 5, "pvf": 6}

# A flash's values per point, and a phase's, in the order reflux_flash
# writes them (REFLUX_FLASH_BETA..., REFLUX_PHASE_FRACTION...).
FLASH_VALUES = ("beta", "t", "p", "h", "s", "v")
PHASE_VALUES = ("fraction", "v", "h", "s")

_i32, _i64, _f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
_doubles = ctypes.POINTER(_f64)
_ints = ctypes.POINTER(_i32)
_text = ctypes.c_char_p
_texts = ctypes.POINTER(ctypes.c_char_p)
_out = [_text, _i32, ctypes.POINTER(_i32)]  # text out: buf, len, needed
_err = [_text, _i32]  # a reason: err, err_len

# Every function in reflux.h: its result's type and its arguments' types.
FUNCTIONS = {
    "reflux_abi_version": (_i32, []),
    "reflux_version": (_i32, [_text, _i32]),
    "reflux_free": (_i32, [_i64, *_err]),
    "reflux_collect": (_i32, [ctypes.POINTER(_i64)]),
    "reflux_data_root": (_i32, _out),
    "reflux_database_open": (_i32, [_text, _texts, _i32, ctypes.POINTER(_i64), *_err]),
    "reflux_database_index": (_i32, [_i64, *_out, *_err]),
    "reflux_compound": (_i32, [_i64, _text, *_out, *_err]),
    "reflux_fluid_new": (_i32, [_i64, _i32, _texts, _i32, _i32, _i32, _texts, _texts,
                                ctypes.POINTER(_f64), _i32, ctypes.POINTER(_i64), *_err]),
    "reflux_fluid_setup": (_i32, [_i64, *_out, *_err]),
    "reflux_props": (_i32, [_i64, _i32, _i64, _doubles, _i64, _doubles, _i64, _doubles, _i64,
                            _doubles, _doubles, _ints, _ints, *_err]),
    "reflux_psat": (_i32, [_i64, _i32, _i64, _doubles, _doubles, _ints, *_err]),
    "reflux_stability": (_i32, [_i64, _i64, _doubles, _i64, _doubles, _i64, _doubles, _i64, _ints,
                                _doubles, _doubles, _ints, *_err]),
    "reflux_flash": (_i32, [_i64, _i32, _i64, _doubles, _i64, _doubles, _i64, _doubles, _i64,
                            _doubles, _doubles, _doubles, _ints, _ints, *_err]),
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
    """An id that isn't a live handle of its kind: never given out, freed
    already, or a fluid's where a database's belongs."""


class NotFound(RefluxError, LookupError):
    """A compound the database doesn't have, or a name that matches more than
    one. The message names the closest."""


class DataError(RefluxError):
    """Component data that can't be read, or can't make what was asked: a
    directory that isn't a data root, a malformed file, or a compound without
    the constants an equation of state needs."""


_ERRORS = {BAD_ARGUMENT: BadArgument, BAD_HANDLE: BadHandle, NOT_FOUND: NotFound,
           DATA_ERROR: DataError}


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


def free_quietly(handle):
    """Drop a handle, ignoring any failure: for a finalizer, which can't
    raise."""
    lib.reflux_free(handle, None, 0)


# Room for text coming out, at first. More is asked for when it's needed.
TEXT_SIZE = 1 << 16


def text_out(f, *args, reason=True):
    """Call f(*args, buf, len, needed[, err, err_len]), a function that writes
    text, and return the text. It's called again with more room if the first
    buffer was too small (reflux.h)."""
    size = TEXT_SIZE
    while True:
        out = ctypes.create_string_buffer(size)
        needed = ctypes.c_int32()
        if reason:
            buf = message_buffer()
            code = f(*args, out, size, ctypes.byref(needed), buf, MESSAGE_SIZE)
        else:
            buf = None
            code = f(*args, out, size, ctypes.byref(needed))
        if code == TOO_SMALL and needed.value > size:
            size = needed.value
            continue
        check(code, buf)
        return out.value.decode("utf-8")


def text_in(s, what):
    """A str as the C API's text: UTF-8 with no NUL inside, which C would take
    for its end."""
    if not isinstance(s, str):
        raise TypeError(f"{what} should be a str, not {type(s).__name__}")
    if "\0" in s:
        raise ValueError(f"{what} has a NUL character")
    return s.encode("utf-8")


def texts_in(items, what):
    """A list of str as the C API's array of texts, and its length."""
    encoded = [text_in(s, f"{what}[{i}]") for i, s in enumerate(items)]
    return (ctypes.c_char_p * len(encoded))(*encoded), len(encoded)
