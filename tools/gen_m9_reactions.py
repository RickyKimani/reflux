# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["mpmath==1.3.0", "chemicals==1.5.2", "thermo==0.6.1", "rdkit==2024.9.6"]
# ///
"""Generate M9's fixtures for formulas and reactions (step 4).

Run from the repository root:

    uv run --python 3.13 tools/gen_m9_reactions.py

It writes two files to testdata/:

m9_formulas.json
    Every component file's composition as two independent programs read it:
    RDKit from its SMILES string (each atom's element and hydrogens, and the
    net charge), and chemicals from its own database by the compound's CAS
    number (the formula PubChem gives). Each is written in Hill's order. Also
    chemicals' atomic weights, element by element, for the standard atomic
    weights reflux holds.

m9_reactions.json
    Standard reaction properties from ChemSep's own data in the component
    files: ΔH°, ΔS°, ΔG° and ln K at 1 bar, and d ln K/dT, for a dozen
    reactions at temperatures from 298.15 to 1500 K. Each species' ideal-gas
    Cp correlation is evaluated here from its coefficients, written
    independently of reflux, and integrated by mpmath at 50 digits:

        H°(T) = ΔH°f + ∫ Cp dT,   S°(T) = S° + ∫ Cp/T dT,   from 298.15 K,
        ΔG° = ΔH° − TΔS°,   ln K = −ΔG°/RT,   d ln K/dT = ΔH°/RT²,

    with R = N_A·k_B exactly. Beside them, Σν ΔG°f from the files, and ln K
    from thermo's own data: chemicals' ΔH°f and S° (at 1 bar) and thermo's
    choice of ideal-gas Cp, which reflux's ChemSep data are measured against.
    And equilibrium constants that others fit, with ln K's distance from
    reflux's at 700 to 1300 K: those the owner's DWSIM case types in, on
    pressures in atm, and those later papers use with Xu and Froment's
    kinetics, on pressures in bar (Xu and Froment's own paper gives none).

The output is the same, byte for byte, for the same versions.
"""

import collections
import json
import pathlib
import tomllib

import chemicals
import mpmath as mp
import rdkit
import thermo
from chemicals import elements, identifiers, reaction
from rdkit import Chem, RDLogger
from thermo import HeatCapacityGas

RDLogger.DisableLog("rdApp.*")
mp.mp.dps = 50

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "components"
OUT = ROOT / "testdata"

R = mp.mpf("8.31446261815324")  # N_A·k_B, exact: 6.02214076e23 × 1.380649e-23
T0 = mp.mpf("298.15")
ATM = mp.mpf(101325) / mp.mpf(100000)  # 1 atm in bar

TEMPERATURES = ["298.15", "300", "400", "500", "600", "700", "800", "900", "1000", "1100", "1200", "1300",
                "1500"]

# Each reaction's equation, as reflux reads it (by formula), and its species
# by CAS number in the equation's order with their coefficients.
REACTIONS = [
    ("steam reforming", "CH4 + H2O = CO + 3 H2", [("74-82-8", -1), ("7732-18-5", -1), ("630-08-0", 1),
                                                   ("1333-74-0", 3)]),
    ("water-gas shift", "CO + H2O = CO2 + H2", [("630-08-0", -1), ("7732-18-5", -1), ("124-38-9", 1),
                                                ("1333-74-0", 1)]),
    ("overall reforming", "CH4 + 2 H2O = CO2 + 4 H2", [("74-82-8", -1), ("7732-18-5", -2), ("124-38-9", 1),
                                                       ("1333-74-0", 4)]),
    ("methanol from CO2", "CO2 + 3 H2 = CH3OH + H2O", [("124-38-9", -1), ("1333-74-0", -3), ("67-56-1", 1),
                                                       ("7732-18-5", 1)]),
    ("methanol from CO", "CO + 2 H2 = CH3OH", [("630-08-0", -1), ("1333-74-0", -2), ("67-56-1", 1)]),
    ("ammonia", "N2 + 3 H2 = 2 NH3", [("7727-37-9", -1), ("1333-74-0", -3), ("7664-41-7", 2)]),
    ("methane combustion", "CH4 + 2 O2 = CO2 + 2 H2O", [("74-82-8", -1), ("7782-44-7", -2), ("124-38-9", 1),
                                                        ("7732-18-5", 2)]),
    ("ethylene hydration", "C2H4 + H2O = C2H5OH", [("74-85-1", -1), ("7732-18-5", -1), ("64-17-5", 1)]),
    ("SO2 oxidation", "2 SO2 + O2 = 2 SO3", [("7446-09-5", -2), ("7782-44-7", -1), ("7446-11-9", 2)]),
    ("butane isomerisation", "n-butane = isobutane", [("106-97-8", -1), ("75-28-5", 1)]),
    ("propane dehydrogenation", "C3H8 = C3H6 + H2", [("74-98-6", -1), ("115-07-1", 1), ("1333-74-0", 1)]),
    ("ethylbenzene dehydrogenation", "Ethylbenzene = Styrene + H2", [("100-41-4", -1), ("100-42-5", 1),
                                                                     ("1333-74-0", 1)]),
]

# Equilibrium constants others fit, ln K = a + b/T: the owner's DWSIM case's,
# as it types them into its rates (spike/dwsim_pfr/run.py), on partial
# pressures in atm; and those later papers use with Xu and Froment's
# kinetics, on pressures in bar, their paper giving K1, K2 and K3 units but
# no values (K3 is K1·K2).
CORRELATIONS = [
    ("DWSIM", "steam reforming", "30.42", "-27106", "atm"),
    ("DWSIM", "water-gas shift", "-3.798", "4160", "atm"),
    ("DWSIM", "overall reforming", "34.218", "-31266", "atm"),
    ("literature", "steam reforming", "30.114", "-26830", "bar"),
    ("literature", "water-gas shift", "-4.036", "4400", "bar"),
    ("literature", "overall reforming", "26.078", "-22430", "bar"),
]
CORRELATION_TEMPERATURES = ["700", "800", "900", "1000", "1100", "1200", "1300"]


def hill(counts):
    """A formula in Hill's order: C, then H, then the rest by symbol; with no
    carbon, every element by symbol."""
    counts = {k: v for k, v in counts.items() if v}
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else []) + sorted(k for k in counts if k not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(k + (str(counts[k]) if counts[k] > 1 else "") for k in order)


def rdkit_formula(smiles):
    """RDKit's reading of a SMILES string: its Hill formula, its net charge,
    and whether RDKit sanitized it. A string RDKit won't sanitize is read
    without, its hydrogens from RDKit's valence model alone."""
    m = Chem.MolFromSmiles(smiles)
    sanitized = m is not None
    if m is None:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return None, None, False
        m.UpdatePropertyCache(strict=False)
    counts = collections.Counter()
    charge = 0
    for a in m.GetAtoms():
        counts[a.GetSymbol()] += 1
        counts["H"] += a.GetTotalNumHs()
        charge += a.GetFormalCharge()
    return hill(counts), charge, sanitized


def chemicals_formula(cas):
    """chemicals' formula for a CAS number, in Hill's order, or None."""
    try:
        f = identifiers.search_chemical(cas).formula
    except Exception:
        return None
    try:
        return hill(elements.nested_formula_parser(f))
    except Exception:
        return None


def formulas():
    compounds = []
    for path in sorted(DATA.glob("*-*.toml")):
        d = tomllib.loads(path.read_text(encoding="utf-8"))
        ids = d["ids"]
        smiles = ids.get("smiles")
        rd, charge, sanitized = rdkit_formula(smiles) if smiles else (None, None, False)
        compounds.append({
            "cas": ids["cas"], "name": d["name"],
            "rdkit": rd, "rdkit_charge": charge, "rdkit_sanitized": sanitized,
            "chemicals": chemicals_formula(ids["cas"]),
        })
    weights = {}
    for z in range(1, 93):
        e = elements.periodic_table[z]
        weights[e.symbol] = e.MW
    return {
        "source": {
            "by": f"tools/gen_m9_reactions.py, with rdkit {rdkit.__version__} and chemicals {chemicals.__version__}",
            "rdkit": "each component file's SMILES read by RDKit: every atom's element, its hydrogens "
                     "(GetTotalNumHs) and formal charge, summed; sanitized where RDKit allows, and read by "
                     "its valence model alone where it doesn't (rdkit_sanitized false)",
            "chemicals": "chemicals.identifiers.search_chemical(cas).formula, from chemicals' database "
                         "(PubChem's formulas), or null where it has no such CAS number",
            "weights": "chemicals.elements.periodic_table[z].MW, g/mol, for z from 1 to 92",
        },
        "weights": weights,
        "compounds": compounds,
    }


def cp_function(eq, k):
    """ChemSep's ideal-gas Cp correlation in J/(mol·K), from its J/(kmol·K)
    coefficients: the forms its files use."""
    A, B, C, D, E = (mp.mpf(float(x)) for x in (list(k) + [0] * 5)[:5])
    if eq == 1:
        return lambda t: A / 1000
    if eq == 100:
        return lambda t: (A + B * t + C * t**2 + D * t**3 + E * t**4) / 1000
    if eq == 16:
        return lambda t: (A + mp.exp(B / t + C + D * t + E * t**2)) / 1000

    def aly_lee(t):
        x, y = C / t, E / t
        return (A + B * (x / mp.sinh(x)) ** 2 + D * (y / mp.cosh(y)) ** 2) / 1000

    if eq == 107:
        return aly_lee
    raise ValueError(eq)


def chemsep(cas):
    d = tomllib.loads((DATA / f"{cas}.toml").read_text(encoding="utf-8"))
    c = d["constants"]
    cp = d["correlations"]["ideal_gas_cp"]
    assert cp["units"] == "J/kmol/K", cas
    return {
        "name": d["name"],
        "dhf": mp.mpf(float(c["dhf"])), "s0": mp.mpf(float(c["s0"])), "dgf": mp.mpf(float(c["dgf"])),
        "cp": cp_function(cp["eq"], cp["coeffs"]),
    }


def standard(species, t):
    """ΔH°, ΔS° of a reaction at T from ChemSep's data, by mpmath."""
    dh = ds = mp.mpf(0)
    for s, nu in species:
        h = s["dhf"] + (mp.quad(s["cp"], [T0, t]) if t != T0 else 0)
        sv = s["s0"] + (mp.quad(lambda x: s["cp"](x) / x, [T0, t]) if t != T0 else 0)
        dh += nu * h
        ds += nu * sv
    return dh, ds


def thermo_ln_k(cases, t):
    """ln K at T from thermo's data: chemicals' ΔH°f and S° at 1 bar and
    thermo's ideal-gas Cp, integrated by thermo."""
    dh = ds = 0.0
    for cas, nu in cases:
        hf, s0 = reaction.Hfg(cas), reaction.S0g(cas)
        if hf is None or s0 is None:
            return None
        cp = HeatCapacityGas(CASRN=cas)
        dh += nu * (hf + cp.T_dependent_property_integral(298.15, t))
        ds += nu * (s0 + cp.T_dependent_property_integral_over_T(298.15, t))
    return -(dh - t * ds) / (8.31446261815324 * t)


def reactions():
    out = []
    by_name = {}
    for name, equation, cases in REACTIONS:
        species = [(chemsep(cas), nu) for cas, nu in cases]
        points = []
        for ts in TEMPERATURES:
            t = mp.mpf(ts)
            dh, ds = standard(species, t)
            dg = dh - t * ds
            ln_k = -dg / (R * t)
            th = thermo_ln_k(cases, float(t))
            points.append({
                "t": float(t), "dh": float(dh), "ds": float(ds), "dg": float(dg), "ln_k": float(ln_k),
                "dln_k_dt": float(dh / (R * t * t)),
                "thermo_ln_k": th,
            })
        dgf = sum(nu * s["dgf"] for s, nu in species)
        dh0, ds0 = standard(species, T0)
        by_name[name] = species
        out.append({
            "name": name, "equation": equation,
            "species": [{"cas": cas, "nu": nu} for cas, nu in cases],
            "thermo_cp": {cas: HeatCapacityGas(CASRN=cas).method for cas, _ in cases},
            "points": points,
            "dgf_sum": float(dgf),
            "gap": float(dgf - (dh0 - T0 * ds0)),
        })
    correlations = []
    for source, name, a, b, unit in CORRELATIONS:
        species = by_name[name]
        dn = sum(nu for _, nu in species)
        diffs = []
        for ts in CORRELATION_TEMPERATURES:
            t = mp.mpf(ts)
            dh, ds = standard(species, t)
            ln_k = -(dh - t * ds) / (R * t)
            # K on pressures in bar is K on pressures in atm times (1.01325)^Δn.
            file_ln_k = mp.mpf(a) + mp.mpf(b) / t + (dn * mp.log(ATM) if unit == "atm" else 0)
            diffs.append(float(ln_k - file_ln_k))
        correlations.append({"source": source, "reaction": name, "a": float(a), "b": float(b), "unit": unit,
                             "dn": dn,
                             "temperatures": [float(t) for t in CORRELATION_TEMPERATURES], "diffs": diffs})
    return {
        "source": {
            "by": f"tools/gen_m9_reactions.py, with mpmath {mp.__version__} at 50 digits, thermo "
                  f"{thermo.__version__}, chemicals {chemicals.__version__}",
            "values": "from ChemSep's data in data/components/ (ΔH°f, S°, ΔG°f, and the ideal-gas Cp "
                      "correlation, evaluated and integrated by mpmath at 50 digits), in J/mol, J/(mol·K), "
                      "and ln K at 1 bar; R = 8.31446261815324 J/(mol·K)",
            "gap": "Σν ΔG°f from the files, less ΔG° at 298.15 K from their ΔH°f and S°, J/mol",
            "thermo_ln_k": "ln K at 1 bar from chemicals' ΔH°f and S° (Hfg, S0g) and thermo's ideal-gas Cp "
                           "(HeatCapacityGas's default method, named in thermo_cp), integrated by thermo; "
                           "null where chemicals lacks a value",
            "correlations": "equilibrium constants others fit, ln K = a + b/T: the owner's DWSIM case's on "
                            "partial pressures in atm (spike/dwsim_pfr/run.py), and the literature's for Xu and "
                            "Froment's kinetics on pressures in bar, the overall reaction's the sum of the other "
                            "two; diffs are reflux's ln K at 1 bar less the correlation's, converted to bar by "
                            "Δn·ln(1.01325) where it's in atm",
        },
        "reactions": out,
        "correlations": correlations,
    }


def write(name, data):
    path = OUT / name
    path.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {path.relative_to(ROOT)}")


def main():
    write("m9_formulas.json", formulas())
    write("m9_reactions.json", reactions())


if __name__ == "__main__":
    main()
