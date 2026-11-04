"""The library itself (M10a step 1): its versions, the functions it exports,
which must be exactly capi/reflux.h's (D3), and the conventions every call
keeps (D2): a reason cut to fit the caller's buffer, a code even with nowhere
to write it, and an id that isn't a live handle an error."""

import ctypes
import re
import subprocess
import sys

import pytest

import reflux
from reflux import _native
from support import ROOT, header_defines, header_functions, v_mod_version


def exported(path):
    """The names in a shared library's table of exports."""
    if sys.platform == "win32":
        out = subprocess.run(["objdump", "-p", str(path)], capture_output=True, text=True,
                             check=True).stdout
        table = out.split("[Ordinal/Name Pointer] Table")[1].split("\n\n")[0]
        return {line.split()[-1] for line in table.splitlines() if "+base[" in line}
    out = subprocess.run(["nm", "-D", "--defined-only", str(path)], capture_output=True, text=True,
                         check=True).stdout
    return {line.split()[-1] for line in out.splitlines() if line.strip()}


def test_the_versions_agree():
    assert reflux.abi_version() == reflux.ABI_VERSION == header_defines()["ABI_VERSION"]
    pyproject = (ROOT / "python" / "pyproject.toml").read_text()
    assert reflux.__version__ == v_mod_version() == re.search(r'^version = "([^"]+)"', pyproject,
                                                              re.M)[1]


def test_a_library_of_another_abi_is_refused(monkeypatch):
    # The package asks the library its ABI version before anything else, so a
    # library built from another reflux.h is refused at import, not misread.
    monkeypatch.setattr(_native, "ABI_VERSION", _native.ABI_VERSION + 1)
    want = f"has ABI version {reflux.ABI_VERSION}, and this package needs {reflux.ABI_VERSION + 1}"
    with pytest.raises(ImportError, match=want):
        _native.load(_native.library_path())


def test_the_codes_are_the_headers():
    # Every number the package uses is the header's, and every one of the
    # header's has its place in the package.
    want = {"ABI_VERSION": _native.ABI_VERSION}
    for name in ("OK", "BAD_ARGUMENT", "BAD_HANDLE", "NOT_FOUND", "DATA_ERROR", "TOO_SMALL"):
        want[name] = getattr(_native, name)
    for name, code in _native.EQUATIONS.items():
        want[name] = code
    for name, code in _native.ALPHAS.items():
        want[f"ALPHA_{name.upper()}"] = code
    for status in reflux.Status:
        want[f"POINT_{status.name}"] = status.value
    for phase in reflux.Phase:
        want["NO_PHASE" if phase.name == "NONE" else phase.name] = phase.value
    for name, code in _native.ROOTS.items():
        want[f"ROOT_{name.upper()}"] = code
    for k, name in enumerate(_native.PROPS):
        want[f"PROP_{name.upper()}"] = k
    want["PROPS"] = len(_native.PROPS)
    assert header_defines() == want


def test_the_library_exports_exactly_the_header():
    # Not the collector's internals, which gcc would export on Windows (M1's
    # spike found 1,600 of them), and not a function the header lacks.
    assert exported(_native.library_path()) == set(header_functions())
    assert set(_native.FUNCTIONS) == set(header_functions())


def test_a_reason_is_cut_to_fit():
    lib = _native.lib
    full = _native.message_buffer()
    assert lib.reflux_free(10**15, full, len(full)) == _native.BAD_HANDLE
    reason = full.value
    assert reason == b"1000000000000000 is not a handle"
    for n in range(1, len(reason) + 3):
        buf = ctypes.create_string_buffer(b"\xff" * (n + 1), n + 1)
        assert lib.reflux_free(10**15, buf, n) == _native.BAD_HANDLE
        kept = min(n - 1, len(reason))
        assert buf.raw[:kept] == reason[:kept] and buf.raw[kept] == 0
        assert buf.raw[n] == 0xFF, f"wrote past a buffer of {n}"
    # Nowhere to write it: the code all the same.
    assert lib.reflux_free(10**15, None, 0) == _native.BAD_HANDLE
    assert lib.reflux_free(10**15, None, 100) == _native.BAD_HANDLE


def test_an_id_that_isnt_a_handle_is_an_error():
    for handle in (0, -1, 2**63 - 1, -(2**63)):
        try:
            _native.free(handle)
        except reflux.BadHandle as e:
            assert str(e) == f"{handle} is not a handle" and e.code == _native.BAD_HANDLE
        else:
            raise AssertionError(f"freeing {handle} raised nothing")


def test_the_version_needs_room():
    lib = _native.lib
    v = reflux.__version__.encode()
    for n in (0, 1, len(v)):
        buf = ctypes.create_string_buffer(b"\xff" * 16, 16)
        assert lib.reflux_version(buf, n) == _native.BAD_ARGUMENT
        assert buf.raw == b"\xff" * 16
    assert lib.reflux_version(None, 64) == _native.BAD_ARGUMENT
    buf = ctypes.create_string_buffer(len(v) + 1)
    assert lib.reflux_version(buf, len(v) + 1) == _native.OK and buf.value == v


def test_collect_gives_the_heap():
    assert reflux.collect() > 0
    assert _native.lib.reflux_collect(None) == _native.OK
