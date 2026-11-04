"""Component data (M10a step 2): a database opened on a data root, its
compounds, and fluids built from them (docs/m10a-plan.md, D8 and D9)."""

import ctypes
import os
import tomllib
import weakref

from . import _native
from .fluid import Fluid

__all__ = ["Database", "data_root"]


def data_root():
    """The data root a Database opens when given none: $REFLUX_DATA if that's
    set, or else the data/ directory of the source tree the library was built
    from. A data root holds components/, a file per compound, and
    vendor/chemsep/, ChemSep's kᵢⱼ banks and Mathias–Copeman constants."""
    return _native.text_out(_native.lib.reflux_data_root, reason=False)


class Database:
    """A data root's compounds, with any directories of your own read after
    them (a compound of yours replaces the root's of the same CAS number), and
    the root's kᵢⱼ banks and Mathias–Copeman constants.

        db = reflux.Database()
        fl = db.fluid("PR", ["methane", "n-butane"])

    Close it with close(), or use it in a `with` block; otherwise it's freed
    when it's garbage. Fluids built from it are independent of it.
    """

    def __init__(self, root=None, extra=()):
        self.root = os.fspath(root) if root is not None else data_root()
        self.extra = [os.fspath(d) for d in extra]
        dirs, n = _native.texts_in(self.extra, "extra")
        handle = ctypes.c_int64()
        buf = _native.message_buffer()
        _native.check(_native.lib.reflux_database_open(
            _native.text_in(self.root, "root"), dirs, n, ctypes.byref(handle), buf,
            _native.MESSAGE_SIZE), buf)
        self._handle = handle.value
        self._finalizer = weakref.finalize(self, _native.free_quietly, self._handle)

    @property
    def handle(self):
        """The library's handle, while the database is open."""
        if not self._finalizer.alive:
            raise ValueError("the database is closed")
        return self._handle

    def close(self):
        """Free the database in the library. Closing it again does nothing."""
        self._finalizer()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __repr__(self):
        state = "" if self._finalizer.alive else ", closed"
        return f"Database({self.root!r}{state})"

    def compound(self, query):
        """A compound's component file, as a dict (docs/data-format.md), found by
        its CAS number, name, alias or InChIKey, ignoring case. NotFound if
        there's none, or more than one; its message names the closest."""
        return tomllib.loads(_native.text_out(_native.lib.reflux_compound, self.handle,
                                              _native.text_in(query, "query")))

    def compounds(self):
        """Every compound the database holds, sorted by CAS number: a dict each,
        with its cas, name, and inchikey and aliases where it has them."""
        index = tomllib.loads(_native.text_out(_native.lib.reflux_database_index, self.handle))
        return [{"cas": cas, **entry} for cas, entry in index.get("compounds", {}).items()]

    def fluid(self, equation, names, *, alpha="equation", kij=True, overrides=None):
        """A fluid of the named compounds for a cubic equation of state.

        equation: "PR", "SRK", "RK" or "vdW", in any case.
        names: the compounds, each by CAS number, name, alias or InChIKey.
        alpha: "equation", the equation's own α (Soave's for SRK and PR), or
            "mathias_copeman", ChemSep's constants, for PR. A compound mc-PR.lib
            lacks keeps the equation's own, and Fluid.setup says so.
        kij: True takes ChemSep's kᵢⱼ for PR or SRK (a pair's first listing),
            False none.
        overrides: kᵢⱼ set by hand, {("methane", "n-butane"): 0.02, ...}.
        """
        code = _native.EQUATIONS.get(str(equation).upper())
        if code is None:
            raise ValueError(f"no equation {equation!r}: PR, SRK, RK or vdW")
        if alpha not in _native.ALPHAS:
            raise ValueError(f"no alpha {alpha!r}: {' or '.join(map(repr, _native.ALPHAS))}")
        if isinstance(names, str):
            raise TypeError("names should be a list of compounds, not one str")
        pairs = list((overrides or {}).items())
        firsts, n = _native.texts_in([a for (a, _), _ in pairs], "an override's first")
        seconds, _ = _native.texts_in([b for (_, b), _ in pairs], "an override's second")
        ks = (ctypes.c_double * n)(*[float(k) for _, k in pairs])
        comps, nc = _native.texts_in(list(names), "names")
        handle = ctypes.c_int64()
        buf = _native.message_buffer()
        _native.check(_native.lib.reflux_fluid_new(
            self.handle, code, comps, nc, _native.ALPHAS[alpha], 1 if kij else 0, firsts,
            seconds, ks, n, ctypes.byref(handle), buf, _native.MESSAGE_SIZE), buf)
        return Fluid._own(handle.value)
