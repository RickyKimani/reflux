"""Many points at once (M10a step 3): inputs broadcast by NumPy's rules and
handed to the library without copies where they can be, and results whose
every point has a status (D6, D7)."""

import ctypes
import enum
import math

import numpy as np

from . import _native

__all__ = ["Phase", "Properties", "Saturation", "Status"]


class Status(enum.IntEnum):
    """A point's status: OK, or why it has no answer (REFLUX_POINT_* in
    reflux.h). A point without an answer has NaN for every value."""

    OK = 0
    BAD_INPUT = 1  # a value it can't use: not finite, T or P ≤ 0, z negative or all 0
    NO_CONVERGENCE = 2  # an iteration ran out of steps
    NO_VOLUME_ROOT = 3  # no volume root at a state it needed
    MORE_PHASES = 4  # a flash's phase is itself unstable
    PHASES_SWAPPED = 5  # a vapour-fraction flash's phases swapped
    OUTSIDE_ENVELOPE = 6  # no such point on the phase envelope
    SUPERCRITICAL = 7  # at or above Tc: no vapour pressure
    NOT_FINITE = 8  # an answer that isn't finite: beyond the model


class Phase(enum.IntEnum):
    """A phase's label, by its phase identification parameter."""

    NONE = -1
    LIQUID = 0
    VAPOR = 1


def doubles(a):
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def ints(a):
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))


def flat(a, shape, n):
    """An input as the library takes it: one value for every point, or n,
    contiguous. A float64 array of the points' shape that's contiguous already
    is passed as it is, not copied."""
    if a.size == 1:
        return np.ascontiguousarray(a.reshape(1)), 1
    if a.shape != shape:
        a = np.broadcast_to(a, shape)
    return np.ascontiguousarray(a).reshape(n), n


def flat_rows(z, shape, n, nc):
    """Compositions as the library takes them: one row for every point, or
    n rows of nc."""
    if math.prod(z.shape[:-1]) == 1:
        return np.ascontiguousarray(z.reshape(1, nc)), 1
    if z.shape[:-1] != shape:
        z = np.broadcast_to(z, shape + (nc,))
    return np.ascontiguousarray(z).reshape(n, nc), n


def floats(x):
    return np.asarray(x, dtype=np.float64)


def call(f, *args):
    """Call a many-point function, whose err holds the first failed point's
    reason, and return that reason."""
    buf = _native.message_buffer()
    _native.check(f(*args, buf, _native.MESSAGE_SIZE), buf)
    return buf.value.decode("utf-8", "replace")


def flat_index(index, shape):
    """A point's index, an int or a tuple with an entry per axis, negative
    from the end, as its position in the flattened points."""
    index = index if isinstance(index, tuple) else (index,)
    if len(index) != len(shape) or not all(-s <= i < s for i, s in zip(index, shape)):
        raise IndexError(f"no point {index} among points of shape {shape}")
    return int(np.ravel_multi_index(tuple(i % s for i, s in zip(index, shape)), shape))


class Properties:
    """A phase's properties at many points, each an array of the points'
    shape, per mole and in SI:

    - z, the compressibility factor; v, the molar volume (m³/mol);
    - h, s, g, u and a: enthalpy, entropy, and Gibbs, internal and Helmholtz
      energy (J/mol, J/(mol·K)), on the ideal gas's formation basis
      (Fluid.absolute says whether every compound has one);
    - cp and cv (J/(mol·K)); w, the speed of sound (m/s); mu_jt, the
      Joule–Thomson coefficient (K/Pa);
    - ln_phi, ln φ, with one more axis, a compound's;
    - phase, the phase's label (Phase), and status (Status), with ok its
      mask. A point without an answer has NaN and Phase.NONE.

    reason(index) says why a point has no answer.
    """

    def __init__(self, fluid, root, shape, inputs, out, ln_phi, phase, status, why):
        self._fluid, self._root, self._inputs, self._why = fluid, root, inputs, why
        for k, name in enumerate(_native.PROPS):
            setattr(self, name, out[k].reshape(shape))
        self.ln_phi = ln_phi.reshape(shape + (ln_phi.shape[-1],))
        self.phase = phase.reshape(shape)
        self.status = status.reshape(shape)

    @property
    def ok(self):
        """Where there's an answer."""
        return self.status == Status.OK

    def reason(self, index):
        """Why the point at index has no answer, or "" if it has one: the
        library's reason for that point alone (D7)."""
        i = flat_index(index, self.status.shape)
        (t, t_len), (p, p_len), (z, z_rows) = self._inputs
        return props(self._fluid, t[i if t_len > 1 else 0], p[i if p_len > 1 else 0],
                     z[i if z_rows > 1 else 0], self._root)._why

    def __repr__(self):
        return f"Properties({self.status.shape}, {int(self.ok.sum())} of {self.status.size} ok)"


def props(fluid, T, P, z, root):
    """Fluid.props' work: broadcast, call, shape."""
    if root not in _native.ROOTS:
        raise ValueError(f"no root {root!r}: {', '.join(map(repr, _native.ROOTS))}")
    nc = fluid.nc
    if z is None:
        if nc != 1:
            raise TypeError(f"z is needed: the fluid has {nc} compounds")
        z = [1.0]
    T, P, z = floats(T), floats(P), floats(z)
    if z.ndim == 0 or z.shape[-1] != nc:
        raise ValueError(f"z's last axis should have the fluid's {nc} compounds, and its "
                         f"shape is {z.shape}")
    shape = np.broadcast_shapes(T.shape, P.shape, z.shape[:-1])
    n = math.prod(shape)
    (t, t_len), (p, p_len) = flat(T, shape, n), flat(P, shape, n)
    zz, z_rows = flat_rows(z, shape, n, nc)
    out = np.empty((len(_native.PROPS), n))
    ln_phi = np.empty((n, nc))
    phase = np.empty(n, dtype=np.int32)
    status = np.empty(n, dtype=np.int32)
    why = call(_native.lib.reflux_props, fluid.handle, _native.ROOTS[root], n, doubles(t), t_len,
               doubles(p), p_len, doubles(zz), z_rows, doubles(out), doubles(ln_phi), ints(phase),
               ints(status))
    return Properties(fluid, root, shape, ((t, t_len), (p, p_len), (zz, z_rows)), out, ln_phi,
                      phase, status, why)


class Saturation:
    """Vapour pressures at many temperatures: p (Pa), with status (Status)
    and ok, its mask. A point without one has NaN. With every compound's,
    the last axis is the compound's. reason(index) says why a point has
    none."""

    def __init__(self, fluid, compound, t, t_shape, p, status):
        # compound is the one compound's index, or None for every compound's,
        # along a last axis.
        self._fluid, self._compound, self._t, self._t_shape = fluid, compound, t, t_shape
        self.p, self.status = p, status

    @property
    def ok(self):
        """Where there's an answer."""
        return self.status == Status.OK

    def reason(self, index):
        """Why the point at index has no vapour pressure, or "" if it has one:
        the library's reason for that compound at that temperature alone."""
        flat_index(index, self.status.shape)
        index = index if isinstance(index, tuple) else (index,)
        if self._compound is None:
            where, j = index[:-1], index[-1] % self._fluid.nc
        else:
            where, j = index, self._compound
        i = flat_index(where, self._t_shape)
        return psat_one(self._fluid, j, self._t[i:i + 1])[2]

    def __repr__(self):
        return f"Saturation({self.status.shape}, {int(self.ok.sum())} of {self.status.size} ok)"


def psat_one(fluid, j, t):
    """One compound's vapour pressure at the temperatures t, contiguous."""
    out = np.empty(t.size)
    status = np.empty(t.size, dtype=np.int32)
    why = call(_native.lib.reflux_psat, fluid.handle, j, t.size, doubles(t), doubles(out),
               ints(status))
    return out, status, why


def psat(fluid, T, compound):
    """Fluid.psat's work."""
    T = floats(T)
    t = np.ascontiguousarray(T).reshape(-1)
    if compound is not None:
        j = fluid.index(compound)
        out, status, _ = psat_one(fluid, j, t)
        return Saturation(fluid, j, t, T.shape, out.reshape(T.shape), status.reshape(T.shape))
    ps, statuses = zip(*(psat_one(fluid, j, t)[:2] for j in range(fluid.nc)))
    shape = T.shape + (fluid.nc,)
    return Saturation(fluid, None, t, T.shape,
                      np.stack(ps, axis=-1).reshape(shape),
                      np.stack(statuses, axis=-1).reshape(shape))
