# /// script
# requires-python = ">=3.10"
# dependencies = ["thermo==0.6.1", "CoolProp==8.0.0"]
# ///
"""Generate the M1 reference fixtures for reflux.

Run from the repository root:

    uv run tools/gen_m1_fixtures.py

uv reads the dependency block above, so there is no virtual environment to
manage. The library versions are pinned, so re-running gives the same numbers.

It writes two files to testdata/:

m1_thermo.json
    The same models as reflux (PR and SRK with van der Waals mixing and the
    exact critical constants), computed by Caleb Bell's `thermo`. reflux must
    agree to about 1e-9: this checks the implementation, not the model.

m1_coolprop_butane.json
    n-butane from CoolProp's reference equation of state (Buecker & Wagner),
    plus the deviations of PR from it at the time of generation, computed with
    `thermo`'s PR. reflux's PR must not deviate by more. This is a regression
    guard on the model's accuracy, not a claim that PR is accurate.

The V tests only read these files. They never run this script.
"""

import datetime
import json
import math
import pathlib
import platform

import CoolProp
import CoolProp.CoolProp as CP
import thermo
from thermo import PR, PRMIX, SRK, SRKMIX

R = 8.314462618  # J/(mol K); the same value as reflux and thermo
OUT = pathlib.Path(__file__).resolve().parent.parent / "testdata"

# Smith, Van Ness & Abbott, Table B.1 (critical constants, acentric factor),
# and molar masses from IUPAC atomic weights. Units are SI, as in reflux.
COMPONENTS = [
    {"name": "methane", "tc": 190.6, "pc": 45.99e5, "omega": 0.012, "mw": 0.016043},
    {"name": "n-butane", "tc": 425.1, "pc": 37.96e5, "omega": 0.200, "mw": 0.058123},
    {"name": "carbon dioxide", "tc": 304.2, "pc": 73.83e5, "omega": 0.224, "mw": 0.044010},
]

KIJ_TERNARY = [[0.0, 0.02, 0.10], [0.02, 0.0, 0.13], [0.10, 0.13, 0.0]]


def provenance():
    return {
        "generated": datetime.date.today().isoformat(),
        "generator": "tools/gen_m1_fixtures.py",
        "python": platform.python_version(),
        "thermo": thermo.__version__,
        "coolprop": CoolProp.__version__,
    }


def phase_values(eos, suffix, mixture):
    """The values reflux checks, for one root ('l' or 'g') of a thermo EOS."""
    ln_phi = getattr(eos, "lnphis_" + suffix) if mixture else [getattr(eos, "lnphi_" + suffix)]
    return {
        "phase": "liquid" if suffix == "l" else "vapor",
        "z": getattr(eos, "Z_" + suffix),
        "v": getattr(eos, "V_" + suffix),
        "ln_phi": list(ln_phi),
        "h_dep": getattr(eos, "H_dep_" + suffix),
        "s_dep": getattr(eos, "S_dep_" + suffix),
        "cp_dep": getattr(eos, "Cp_dep_" + suffix),
    }


def thermo_case(kind, comp_idx, zs, t, p):
    comps = [COMPONENTS[i] for i in comp_idx]
    mixture = len(comps) > 1
    if mixture:
        cls = PRMIX if kind == "pr" else SRKMIX
        kij = [[KIJ_TERNARY[i][j] for j in comp_idx] for i in comp_idx]
        eos = cls(
            Tcs=[c["tc"] for c in comps],
            Pcs=[c["pc"] for c in comps],
            omegas=[c["omega"] for c in comps],
            zs=zs,
            kijs=kij,
            T=t,
            P=p,
        )
    else:
        cls = PR if kind == "pr" else SRK
        c = comps[0]
        kij = [[0.0]]
        eos = cls(Tc=c["tc"], Pc=c["pc"], omega=c["omega"], T=t, P=p)
    phases = [phase_values(eos, s, mixture) for s in ("l", "g") if hasattr(eos, "Z_" + s)]
    return {
        "eos": kind,
        "components": comp_idx,
        "kij": [x for row in kij for x in row],
        "n": zs,
        "t": t,
        "p": p,
        "phases": phases,
    }


def thermo_fixture():
    cases = []
    for kind in ("pr", "srk"):
        for t, p in [(300.0, 2e5), (350.0, 9.4573e5), (400.0, 20e5), (450.0, 50e5)]:
            cases.append(thermo_case(kind, [1], [1.0], t, p))
        for t, p in [(300.0, 10e5), (250.0, 60e5), (400.0, 20e5)]:
            cases.append(thermo_case(kind, [0, 1, 2], [0.2, 0.5, 0.3], t, p))
    return {
        "source": provenance(),
        "note": "Same models as reflux, computed by thermo. Residuals are departures "
        "at the same T and P. Per mole of mixture (n sums to 1).",
        "components": COMPONENTS,
        "cases": cases,
    }


def pr_root(eos, phase):
    """The requested root, or the only one if the state has a single root."""
    want, other = ("l", "g") if phase == "liquid" else ("g", "l")
    return want if hasattr(eos, "Z_" + want) else other


def coolprop_fixture():
    butane = COMPONENTS[1]
    ref = CP.AbstractState("HEOS", "n-Butane")
    points = []
    for t, p, phase in [
        (350.0, 2e5, "vapor"),
        (400.0, 5e5, "vapor"),
        (450.0, 50e5, "vapor"),
        (430.0, 100e5, "liquid"),
        (300.0, 30e5, "liquid"),
        (250.0, 10e5, "liquid"),
    ]:
        ref.update(CP.PT_INPUTS, p, t)
        rho = ref.rhomolar()
        z = ref.compressibility_factor()
        cp0 = ref.cp0molar()
        point = {
            "t": t,
            "p": p,
            "phase": phase,
            "v": 1.0 / rho,
            "z": z,
            # Enthalpy departure is the same on either basis.
            "h_dep": ref.hmolar_residual(),
            # CoolProp's smolar_residual is S − S_ig at the same T and density;
            # adding R ln Z puts it on reflux's (T, P) basis.
            "s_dep": ref.smolar_residual() + R * math.log(z),
            "cp0": cp0,
            "cp": ref.cpmolar(),
            "cv": ref.cvmolar(),
            "w": ref.speed_sound(),
        }

        # PR with the same component constants reflux uses, and CoolProp's
        # ideal-gas Cp, so that only the residual model differs.
        eos = PR(Tc=butane["tc"], Pc=butane["pc"], omega=butane["omega"], T=t, P=p)
        s = pr_root(eos, phase)
        v = getattr(eos, "V_" + s)
        cp = cp0 + getattr(eos, "Cp_dep_" + s)
        cv = cp - getattr(eos, "Cp_minus_Cv_" + s)
        w = math.sqrt(-v * v / butane["mw"] * cp / cv * getattr(eos, "dP_dV_" + s))
        point["pr_deviation"] = {
            "v_rel": (v - point["v"]) / point["v"],
            "h_dep": getattr(eos, "H_dep_" + s) - point["h_dep"],
            "s_dep": getattr(eos, "S_dep_" + s) - point["s_dep"],
            "cp_rel": (cp - point["cp"]) / point["cp"],
            "w_rel": (w - point["w"]) / point["w"],
        }
        points.append(point)
    return {
        "source": provenance(),
        "note": "n-butane from CoolProp's reference EOS (HEOS backend). pr_deviation is "
        "PR (constants below, CoolProp's ideal-gas Cp) minus the reference, computed "
        "with thermo. v, h_dep and s_dep are molar; *_rel are relative.",
        "component": butane,
        "points": points,
    }


def write(name, data):
    path = OUT / name
    path.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print("wrote", path.relative_to(OUT.parent))


if __name__ == "__main__":
    write("m1_thermo.json", thermo_fixture())
    write("m1_coolprop_butane.json", coolprop_fixture())
