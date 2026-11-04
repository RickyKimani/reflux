"""Fluids (M10a step 2): a cubic equation of state over compounds from a
Database, which later steps flash and evaluate on NumPy arrays."""

import tomllib
import weakref

import numpy as np

from . import _native

__all__ = ["Fluid"]


def _frozen(values):
    a = np.array(values, dtype=np.float64)
    a.flags.writeable = False
    return a


class Fluid:
    """A cubic equation of state over compounds from a Database: make one with
    Database.fluid. Its compounds' constants and kᵢⱼ are read-only arrays, in
    SI, and setup says how it was built.

    Close it with close(), or use it in a `with` block; otherwise it's freed
    when it's garbage.
    """

    def __init__(self):
        raise TypeError("make a Fluid with Database.fluid")

    @classmethod
    def _own(cls, handle):
        """A Fluid that owns a handle reflux_fluid_new gave."""
        self = cls.__new__(cls)
        self._handle = handle
        self._finalizer = weakref.finalize(self, _native.free_quietly, handle)
        setup = tomllib.loads(_native.text_out(_native.lib.reflux_fluid_setup, handle))
        comps = setup["components"]
        self.setup = setup
        self.equation = setup["equation"]
        self.names = [c["name"] for c in comps]
        self.cas = [c["cas"] for c in comps]
        self.tc = _frozen([c["tc"] for c in comps])
        self.pc = _frozen([c["pc"] for c in comps])
        self.omega = _frozen([c["omega"] for c in comps])
        self.molar_mass = _frozen([c["mw"] for c in comps])
        kij = np.zeros((len(comps), len(comps)))
        for p in setup.get("pairs", []):
            kij[p["i"], p["j"]] = kij[p["j"], p["i"]] = p["k"]
        kij.flags.writeable = False
        self.kij = kij
        return self

    @property
    def nc(self):
        """The number of compounds."""
        return len(self.names)

    @property
    def handle(self):
        """The library's handle, while the fluid is open."""
        if not self._finalizer.alive:
            raise ValueError("the fluid is closed")
        return self._handle

    def close(self):
        """Free the fluid in the library. Closing it again does nothing."""
        self._finalizer()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __repr__(self):
        state = "" if self._finalizer.alive else ", closed"
        return f"Fluid({self.equation!r}, {self.names!r}{state})"
