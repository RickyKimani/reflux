"""Two phase diagrams from reflux, on NumPy arrays: M10a's exit script
(docs/m10a-plan.md, exit criterion 8).

- A T–x–y diagram of methane and n-butane at 20 bar, below the binary's
  critical locus, so that every composition has a bubble and a dew point:
  101 of each, with tie lines from the bubble points to their first bubbles.
- A phase-fraction map of a natural gas: PT flashes over 81 temperatures and
  80 pressures, up to 120 bar, past its cricondenbar, the vapour fraction
  coloured and the two-phase region outlined.

Both are Peng–Robinson with ChemSep's kᵢⱼ. With the library built (`v run
tools/python.vsh build`) and the package importable:

    python python/examples/phase_diagrams.py [figure.png]

It writes the figure, phase_diagrams.png unless named, and needs matplotlib.
python/tests/test_examples.py runs it in CI, and checks its numbers against
the V API's bit for bit.
"""

import sys

import numpy as np

import reflux

# The T–x–y. Methane's mole fraction is i/100, exactly as V computes it, so
# the V reference (tools/capi_reference.v) flashes the same compositions.
TXY_NAMES = ["methane", "n-butane"]
TXY_P = 20e5  # Pa
X = np.arange(101) / 100

# The map. Every T and P is exact in binary, for the same reason.
GAS_NAMES = ["methane", "ethane", "propane", "n-butane", "n-pentane"]
GAS_Z = np.array([0.8, 0.1, 0.05, 0.03, 0.02])
MAP_T = 150.0 + 2.5 * np.arange(81)  # K, 150 to 350
MAP_P = 1.5e5 * np.arange(1, 81)  # Pa, 1.5 to 120 bar


def txy(db):
    """The T–x–y: the bubble and dew points at TXY_P of each composition in
    X, as two FlashResults. A bubble point's vapour is its first bubble, and a
    dew point's liquid its first drop."""
    fluid = db.fluid("PR", TXY_NAMES)
    z = np.stack([X, 1 - X], axis=-1)
    return fluid.bubble_t(TXY_P, z), fluid.dew_t(TXY_P, z)


def phase_map(db):
    """The gas's PT flash at every (P, T) of the grid: one FlashResult of
    shape (80, 81), a row per pressure."""
    fluid = db.fluid("PR", GAS_NAMES)
    return fluid.flash_pt(MAP_T, MAP_P[:, None], GAS_Z)


def plot(bubble, dew, gas, path):
    """Draw both diagrams side by side into path."""
    import matplotlib
    matplotlib.use("Agg")  # a file, never a window
    import matplotlib.pyplot as plt

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")

    # T against x on the bubble curve and against y on the dew curve: a dew
    # point's vapour is its feed, so both curves are drawn over X.
    left.plot(X, bubble.t, color="C0", label="bubble point, T(x)")
    left.plot(X, dew.t, color="C1", label="dew point, T(y)")
    for i in range(10, 100, 10):
        y = bubble.vapor.x[i, 0]
        left.plot([X[i], y], [bubble.t[i]] * 2, color="0.6", lw=0.8, zorder=0)
    left.plot([], [], color="0.6", lw=0.8, label="tie line, to the first bubble")
    left.set(xlabel="methane mole fraction", ylabel="T (K)", xlim=(0, 1),
             title=f"methane and n-butane at {TXY_P / 1e5:g} bar")
    left.legend(frameon=False, loc="upper right")

    mesh = right.pcolormesh(MAP_T, MAP_P / 1e5, gas.beta, cmap="viridis", vmin=0, vmax=1,
                            shading="nearest")
    right.contour(MAP_T, MAP_P / 1e5, gas.phases, levels=[1.5], colors="white", linewidths=1)
    fig.colorbar(mesh, ax=right, label="vapour fraction β")
    right.set(xlabel="T (K)", ylabel="P (bar)",
              title="a natural gas, 80% methane: two phases inside the line")

    fig.suptitle("Peng–Robinson, ChemSep's kᵢⱼ, by reflux")
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main(argv):
    path = argv[1] if len(argv) > 1 else "phase_diagrams.png"
    with reflux.Database() as db:
        bubble, dew = txy(db)
        gas = phase_map(db)
    for name, r in (("bubble points", bubble), ("dew points", dew), ("map", gas)):
        failed = np.argwhere(~r.ok)
        print(f"{name}: {r.ok.sum()} of {r.ok.size} have answers", end="")
        print(f"; the first without, {tuple(failed[0])}: {r.reason(tuple(failed[0]))}"
              if failed.size else "")
    print(f"map: {(gas.phases == 2).sum()} of {gas.phases.size} points split")
    plot(bubble, dew, gas, path)
    print(f"wrote {path}")


if __name__ == "__main__":
    main(sys.argv)
