"""Write spike/seam/seamflash/: flash's PT flash and stability test with
thermo.Model replaced by a phase-model interface, Phases (M6's D3). The
copy is mechanical, so that its only difference from flash is the seam.

    python spike/seam/make.py

Run from the repository root. Not part of the simulator.
"""
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "spike" / "seam" / "seamflash"
OUT.mkdir(exist_ok=True)

REPLACE = [
    ("module flash", "module seamflash"),
    ("thermo.solve_tp(sp.m, sp.t, sp.p, x, .stable)", "sp.m.phase(sp.t, sp.p, x, .stable)"),
    ("thermo.solve_tp(f.m, f.t, f.p, x, .stable)", "f.m.phase(f.t, f.p, x, .stable)"),
    ("thermo.solve_tp(m, t, p, feed_x, .stable)", "m.phase(t, p, feed_x, .stable)"),
    ("thermo.solve_tp(m, t, p, feed, .stable)", "m.phase(t, p, feed, .stable)"),
    ("thermo.ln_phi_derivs(sp.m, sy)", "sp.m.derivs(sy)"),
    ("thermo.ln_phi_derivs(sp.m, sx)", "sp.m.derivs(sx)"),
    ("thermo.ln_phi_derivs(f.m, cur.st)", "f.m.derivs(cur.st)"),
    ("thermo.phase_label(m, st)", "m.label(st)"),
    ("thermo.Model", "Phases"),
]

for name in ["pt.v", "stability.v", "wilson.v", "error.v"]:
    text = (ROOT / "flash" / name).read_text(encoding="utf-8")
    for a, b in REPLACE:
        text = text.replace(a, b)
    left = [w for w in ["thermo.solve_tp", "thermo.ln_phi_derivs", "thermo.phase_label", "thermo.Model"] if w in text]
    assert not left, (name, left)
    (OUT / name).write_bytes(text.encode("utf-8"))

(OUT / "phases.v").write_bytes('''module seamflash

import thermo

// Phases is M6's D3 in its least form: what the PT flash and the stability
// test ask of a model, a phase at (T, P, n) with its ln φ, their derivatives
// and its label. Not part of the simulator.
pub interface Phases {
	components() []thermo.Component
	phase(t f64, p f64, n []f64, root thermo.Root) !thermo.State
	derivs(st thermo.State) thermo.LnPhiDerivs
	label(st thermo.State) thermo.PhaseLabel
}

// Eos is a Helmholtz model as a Phases, by its volume roots: exactly what
// flash calls on a thermo.Model, one call further away.
pub struct Eos {
pub:
	m thermo.Model
}

pub fn (e Eos) components() []thermo.Component {
	return e.m.components()
}

pub fn (e Eos) phase(t f64, p f64, n []f64, root thermo.Root) !thermo.State {
	return thermo.solve_tp(e.m, t, p, n, root)
}

pub fn (e Eos) derivs(st thermo.State) thermo.LnPhiDerivs {
	return thermo.ln_phi_derivs(e.m, st)
}

pub fn (e Eos) label(st thermo.State) thermo.PhaseLabel {
	return thermo.phase_label(e.m, st)
}
'''.encode("utf-8"))
print("wrote", OUT.relative_to(ROOT))
