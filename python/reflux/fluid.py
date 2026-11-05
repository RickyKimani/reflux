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
    def absolute(self):
        """Whether every compound has its ideal gas's ΔH°f and S°, so that H, S,
        G, U and A are on the formation basis. Without, they sit on an
        arbitrary zero, and only differences at one composition mean
        anything."""
        return all("dhf" in c and "s0" in c for c in self.setup["components"])

    def index(self, compound):
        """A compound's position: given as one, or by its name or CAS number."""
        if isinstance(compound, (int, np.integer)):
            if not -self.nc <= compound < self.nc:
                raise IndexError(f"no compound {compound}: the fluid has {self.nc}")
            return int(compound) % self.nc
        for j, (name, cas) in enumerate(zip(self.names, self.cas)):
            if compound in (cas, name) or str(compound).casefold() == name.casefold():
                return j
        raise KeyError(f"{compound!r} isn't in the fluid: {', '.join(self.names)}")

    def props(self, T, P, z=None, *, root="stable"):
        """A phase's properties at many points: a Properties of arrays, per mole
        and in SI.

        T (K), P (Pa) and z broadcast together by NumPy's rules, z's last axis
        being the compounds' mole amounts (normalised by their sum), so one
        composition can serve every point. A pure fluid needs no z. Arrays
        already of the points' shape, C-contiguous float64, aren't copied.
        root is "liquid" (the smallest volume), "vapor" (the largest) or
        "stable" (the lower Gibbs energy); where there's one root, it's that.
        """
        from .points import props
        return props(self, T, P, z, root)

    def stability(self, T, P, z=None):
        """Whether each feed is stable as one phase at T (K) and P (Pa): a
        Stability, with tm and the phase that would appear. T, P and z
        broadcast as for props."""
        from .points import stability
        return stability(self, T, P, z)

    def flash_pt(self, T, P, z=None):
        """The PT flash at many points: a FlashResult, its arrays of the
        points' shape. T (K), P (Pa) and z broadcast as for props."""
        from .points import flash
        return flash(self, "pt", T, P, z)

    def psat(self, T, compound=None):
        """Vapour pressures (Pa) at temperatures T (K): a Saturation, with each
        compound's along a last axis, or one compound's (an index, name or
        CAS number). It's the fluid's equation and α for the compound alone;
        at or above its Tc there's none (Status.SUPERCRITICAL)."""
        from .points import psat
        return psat(self, T, compound)

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
