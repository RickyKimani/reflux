# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5"]
# ///
"""Write M6's DWSIM references (step 12, D17, criterion 11):
testdata/m6_dwsim.json.

DWSIM 9.0.5's activity packages, by script, every one at its defaults (an
ideal-gas vapour, Poynting on, ChemSep's liquid densities, Lee-Kesler
enthalpies), on the compounds of its ChemSep database:

- each package's gamma by state, methanol-water at 300, 330 and 360 K and
  x1 = 0.1 to 0.9, from its liquid phi with its own Poynting term taken out,
  gamma = phi P / (Psat exp(V_L |P - Psat|/RT)), V_L = M/rho_L(T) by its
  own Psat and density: DWSIM's Poynting term has |P - Psat| where the
  textbook's has P - Psat, so above Psat it's the textbook's and below it
  it raises phi by as much as it should lower it (found here: methanol at
  360 K and 1 atm). gamma_signed is the textbook's, for the record;
- each compound's vapour-pressure and density constants as DWSIM holds them;
- UNIQUAC's two records for methanol-water, Wilson's parameters and
  UNIFAC-LL's groups, as its packages hold them, and its gamma worked out
  again here from its own code (DanWBR/dwsim, PropertyPackages/Models/
  UNIQUAC.vb, WilsonModel.vb and UNIFAC.vb with PropertyPackages/
  UNIFACLL.vb): UNIQUAC with tau_ij of each ordered pair from that pair's
  own record; Wilson with R = 1.9872 and V_j/V_i at 298.15 K; UNIFAC-LL as
  its VLE UNIFAC with each compound's group fraction e_ki = nu_k Q_k over
  the Q of its Dortmund group numbers read in the VLE subgroup table;
- NRTL and UNIQUAC PT flashes of methanol-water at 1 atm, with its loops at
  1e-10, and its Psat and density at each T;
- its three-phase flash of water-1-butanol by NRTL, at 0.5 and 0.8 water,
  by each of its flash approaches and with its three-phase stability test
  asked for each compound;
- a pure liquid's enthalpy less its ideal gas's at 330 K and 1 atm, by each
  of its enthalpy modes.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m6_dwsim.py

On Windows with DWSIM 9 installed. About two minutes. Its numbers reproduce
byte for byte on the same day: the three-phase flash's, which move from run
to run, are kept to 10 significant figures by NestedLoops and 3 by the
others, whose answers aren't converged, and with no count of iterations.
"""
import datetime
import json
import math
import os
import pathlib
import sys

import pythonnet

pythonnet.load("netfx")
import clr  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
HOME = os.path.join(os.environ["LOCALAPPDATA"], "DWSIM")
sys.path.append(HOME)
os.chdir(HOME)  # DWSIM finds its databases relative to the working directory
for name in ["CapeOpen", "DWSIM.Automation", "DWSIM.Interfaces", "DWSIM.GlobalSettings", "DWSIM.SharedClasses",
             "DWSIM.Thermodynamics", "DWSIM.UnitOperations", "DWSIM.Inspector", "System.Buffers"]:
    clr.AddReference(name)

from System import Array, Double, Enum  # noqa: E402
from DWSIM.Automation import Automation3  # noqa: E402
from DWSIM.Interfaces.Enums import FlashCalculationType, FlashSetting  # noqa: E402
from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType  # noqa: E402
from DWSIM.Thermodynamics.PropertyPackages import State  # noqa: E402

P = 101325.0
R_SI = 8.314462618  # DWSIM's own R in its Poynting term
TOL = "1E-10"
vec = lambda xs: Array[Double](list(xs))  # noqa: E731
auto = Automation3()
VERSION = str(auto.GetVersion())
PACKAGES = ["NRTL", "UNIQUAC", "Wilson", "UNIFAC", "Modified UNIFAC (Dortmund)", "UNIFAC-LL"]
BUTANOL = next(k for k in auto.AvailableCompounds.Keys if k.lower() in ("1-butanol", "n-butanol"))


def figures(value, n):
    """Every float in a flash's result to n significant figures: DWSIM's
    three-phase flash moves from run to run, its converged answers by about
    1e-13 (its stability test's random trials) and InsideOut's and
    GibbsMinimization's unconverged ones by up to 1e-4, and the file must
    reproduce byte for byte."""
    if isinstance(value, float):
        return float(f"{value:.{n}g}")
    if isinstance(value, list):
        return [figures(v, n) for v in value]
    if isinstance(value, dict):
        return {k: figures(v, n) for k, v in value.items()}
    return value


def set_enum(obj, attr, name):
    """Set an enum property by its member's name, on the property's own type."""
    setattr(obj, attr, Enum.Parse(getattr(obj, attr).GetType(), name))


class Package:
    """A flowsheet and one property package, its PT loops at tol, or at
    DWSIM's defaults for none (tools/gen_m6_bench.py times those)."""

    def __init__(self, compounds, package, z, tol=TOL):
        sim = auto.CreateFlowsheet()
        for c in compounds:
            comp = sim.AvailableCompounds[c]
            assert comp.OriginalDB == "ChemSep", (c, comp.OriginalDB)
            sim.SelectedCompounds.Add(c, comp)
        self.sim, self.compounds = sim, compounds
        self.pp = sim.CreateAndAddPropertyPackage(package).__implementation__
        s = sim.AddObject(ObjectType.MaterialStream, 50, 50, "feed").GetAsObject()
        s.SetTemperature(330.0)
        s.SetPressure(P)
        s.SetMolarFlow(1.0)
        s.SetOverallComposition(vec(z))
        self.pp.CurrentMaterialStream = s
        if tol is not None:
            for key in (FlashSetting.PTFlash_External_Loop_Tolerance, FlashSetting.PTFlash_Internal_Loop_Tolerance):
                self.pp.FlashSettings[key] = tol

    def comp(self, i):
        return self.sim.SelectedCompounds[self.compounds[i]]

    def psat(self, i, t):
        return self.pp.AUX_PVAPi(i, t)

    def rho(self, i, t):
        return self.pp.AUX_LIQDENSi(self.comp(i), t)

    def gamma(self, x, t, signed=False):
        """gamma from the liquid's phi, its Poynting term taken out on its own
        Psat and density: as DWSIM writes it, with |P - Psat|, or signed."""
        phi = list(self.pp.DW_CalcFugCoeff(vec(x), t, P, State.Liquid))
        out = []
        for i in range(len(x)):
            ps = self.psat(i, t)
            vl = self.comp(i).Molar_Weight / 1000.0 / self.rho(i, t)
            dp = P - ps if signed else abs(P - ps)
            out.append(phi[i] * P / ps / math.exp(vl * dp / (R_SI * t)))
        return out

    def flash(self, t, z):
        r = self.pp.CalculateEquilibrium(FlashCalculationType.PressureTemperature, P, t, vec(z), None, 0.0)
        l2 = r.GetLiquidPhase2MoleFraction()
        return {"beta": r.GetVaporPhaseMoleFraction(), "y": list(r.GetVaporPhaseMoleFractions()),
                "liquid1": r.GetLiquidPhase1MoleFraction(), "x1": list(r.GetLiquidPhase1MoleFractions()),
                "liquid2": l2, "x2": list(r.GetLiquidPhase2MoleFractions()) if l2 > 0 else [],
                "iterations": int(r.IterationsTaken)}


def compound_data(pk):
    out = []
    for i, name in enumerate(pk.compounds):
        c = pk.comp(i)
        out.append({
            "name": name, "cas": c.CAS_Number, "mw": c.Molar_Weight / 1000.0,
            "psat_eq": int(c.VaporPressureEquation),
            "psat": [c.Vapor_Pressure_Constant_A, c.Vapor_Pressure_Constant_B, c.Vapor_Pressure_Constant_C,
                     c.Vapor_Pressure_Constant_D, c.Vapor_Pressure_Constant_E],
            "density_eq": int(c.LiquidDensityEquation),
            "density": [c.Liquid_Density_Const_A, c.Liquid_Density_Const_B, c.Liquid_Density_Const_C,
                        c.Liquid_Density_Const_D, c.Liquid_Density_Const_E],
            "uniquac_r": c.UNIQUAC_R, "uniquac_q": c.UNIQUAC_Q,
            "rho_298": pk.rho(i, 298.15),
        })
    return out


def uniquac_records(pk):
    out = {}
    ip = pk.pp.m_uni.InteractionParameters
    for a in pk.compounds:
        for b in pk.compounds:
            if a != b and ip.ContainsKey(a) and ip[a].ContainsKey(b):
                d = ip[a][b]
                out[f"{a}|{b}"] = {"A12": d.A12, "A21": d.A21, "B12": d.B12, "B21": d.B21, "C12": d.C12,
                                   "C21": d.C21, "comment": str(d.comment)}
    return out


def uniquac_dwsim(x, t, r, q, recs, names):
    """DWSIM's GAMMA_MR: each ordered pair (i, j) takes a12, a21 from its own
    record [i][j], else [j][i] swapped; tau_ij = exp(-a12/(1.98721 T)),
    tau_ji = exp(-a21/(1.98721 T)); S_i = sum_j theta_j tau_ji(i)(j),
    sum1_i = sum_j theta_j tau_ij(i)(j)/S_j; z = 10."""
    n = len(x)
    a12 = [[0.0] * n for _ in range(n)]
    a21 = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if f"{names[i]}|{names[j]}" in recs:
                d = recs[f"{names[i]}|{names[j]}"]
                a12[i][j], a21[i][j] = d["A12"], d["A21"]
            elif f"{names[j]}|{names[i]}" in recs:
                d = recs[f"{names[j]}|{names[i]}"]
                a12[i][j], a21[i][j] = d["A21"], d["A12"]
    tij = [[math.exp(-a12[i][j] / (1.98721 * t)) for j in range(n)] for i in range(n)]
    tji = [[math.exp(-a21[i][j] / (1.98721 * t)) for j in range(n)] for i in range(n)]
    rr = sum(xi * v for xi, v in zip(x, r))
    qq = sum(xi * v for xi, v in zip(x, q))
    fi = [xi * v / rr for xi, v in zip(x, r)]
    te = [xi * v / qq for xi, v in zip(x, q)]
    s = [sum(te[j] * tji[i][j] for j in range(n)) for i in range(n)]
    s1 = [sum(te[j] * tij[i][j] / s[j] for j in range(n)) for i in range(n)]
    out = []
    for i in range(n):
        lngc = 1 - r[i] / rr + math.log(r[i] / rr) - 5.0 * q[i] * (1 - fi[i] / te[i] + math.log(fi[i] / te[i]))
        lngr = q[i] * (1 - math.log(s[i]) - s1[i])
        out.append(math.exp(lngc + lngr))
    return out


def wilson_bips(pk):
    """Wilson's parameters as DWSIM's package gets them, cal/mol, i with j."""
    m = pk.pp.WilsonM
    out = {}
    for i, a in enumerate(pk.compounds):
        for j, b in enumerate(pk.compounds):
            if i != j:
                out[f"{a}|{b}"] = float(m.GetBIPs(pk.comp(i).CAS_Number, pk.comp(j).CAS_Number))
    return out


def wilson_dwsim(x, t, bips, v, names):
    """DWSIM's Wilson: Lambda_ij = (V_j/V_i) exp(-BIP_ij/(1.9872 T))."""
    n = len(x)
    lam = [[1.0 if i == j else v[j] / v[i] * math.exp(-bips[f"{names[i]}|{names[j]}"] / (1.9872 * t))
            for j in range(n)] for i in range(n)]
    out = []
    for i in range(n):
        s = sum(x[j] * lam[i][j] for j in range(n))
        t2 = sum(x[k] * lam[k][i] / sum(x[j] * lam[k][j] for j in range(n)) for k in range(n))
        out.append(math.exp(-math.log(s) + 1 - t2))
    return out


def unifac_ll_dwsim(x, t, r, q, e, a):
    """DWSIM's UNIFAC for molecules of one group each, molecule i's group i
    with fraction e_i: theta_k = x_k q_k e_k / sum x q, s_k = sum_m theta_m
    psi_mk, ln gamma_i^R = q_i (1 - sum_k theta_k beta_ik/s_k + e_i
    ln(beta_ii/s_i)), beta_ik = e_i psi_ik, psi_mk = exp(-a_mk/T); the
    combinatorial term with z = 10."""
    n = len(x)
    psi = [[1.0 if m == k else math.exp(-a[m][k] / t) for k in range(n)] for m in range(n)]
    rr = sum(x[i] * r[i] for i in range(n))
    qq = sum(x[i] * q[i] for i in range(n))
    theta = [x[k] * q[k] * e[k] / qq for k in range(n)]
    s = [sum(theta[m] * psi[m][k] for m in range(n)) for k in range(n)]
    out = []
    for i in range(n):
        beta = [e[i] * psi[i][k] for k in range(n)]
        tmp = sum(theta[k] * beta[k] / s[k] for k in range(n)) - e[i] * math.log(beta[i] / s[i])
        j, l = r[i] / rr, q[i] / qq
        c = 1 - j + math.log(j) - 5 * q[i] * (1 - j / l + math.log(j / l))
        out.append(math.exp(c + q[i] * (1 - tmp)))
    return out


def unifac_ll_data(pk):
    """UNIFAC-LL's groups as DWSIM holds them: each compound's UNIFAC and
    Dortmund group numbers, its UNIFAC subgroup's R and Q, the Q its Dortmund
    number has in the VLE subgroup table, its group fraction (RET_VEKI), and
    the main groups' a."""
    groups = pk.pp.m_uni.UnifGroups.Groups

    def group(k):
        for key in groups.Keys:
            if str(key) == str(k):
                return groups[key]
        raise KeyError(k)

    out = {"compounds": []}
    veki = pk.pp.RET_VEKI()
    for i in range(len(pk.compounds)):
        cp = pk.comp(i)
        uv = {str(k): cp.UNIFACGroups[k] for k in cp.UNIFACGroups.Keys}
        mf = {str(k): cp.MODFACGroups[k] for k in cp.MODFACGroups.Keys}
        assert len(uv) == 1 and len(mf) == 1
        g = group(next(iter(uv)))
        out["compounds"].append({"unifac": uv, "modfac": mf, "main": str(g.PrimaryGroup), "R": g.R, "Q": g.Q,
                                 "modfac_as_unifac_Q": group(next(iter(mf))).Q,
                                 "e": [veki[i][k] for k in veki[i].Keys][0]})
    ip = pk.pp.m_uni.UnifGroups.InteracParam
    mains = [c["main"] for c in out["compounds"]]
    a = [[0.0] * len(mains) for _ in mains]
    for m, gm in enumerate(mains):
        for k, gk in enumerate(mains):
            if m != k:
                for ka in ip.Keys:
                    if str(ka) == gm:
                        for kb in ip[ka].Keys:
                            if str(kb) == gk:
                                a[m][k] = ip[ka][kb]
    out["a"] = a
    return out


def main():
    mw = ["Methanol", "Water"]
    states = [(t, [x1, 1 - x1]) for t in (300.0, 330.0, 360.0) for x1 in (0.1, 0.3, 0.5, 0.7, 0.9)]
    out = {"source": {"generator": "tools/gen_m6_dwsim.py", "dwsim": VERSION,
                      "generated": datetime.date.today().isoformat()},
           "note": ("DWSIM 9.0.5's activity packages at their defaults on its ChemSep compounds. gamma has DWSIM's "
                    "Poynting term taken out on its own Psat and density; the _formula values are DWSIM's code "
                    "worked out again here on the parameters it holds."),
           "P": P, "methanol_water": {}, "packages": {}}
    nrtl = Package(mw, "NRTL", [0.5, 0.5])
    out["methanol_water"]["compounds"] = compound_data(nrtl)
    for package in PACKAGES:
        pk = Package(mw, package, [0.5, 0.5])
        entry = {"class": pk.pp.GetType().FullName,
                 "gamma": [{"T": t, "x": x, "gamma": pk.gamma(x, t), "gamma_signed": pk.gamma(x, t, True)}
                           for t, x in states],
                 "psat": [{"T": t, "psat": [pk.psat(i, t) for i in range(2)]} for t in (300.0, 330.0, 360.0)],
                 "rho": [{"T": t, "rho": [pk.rho(i, t) for i in range(2)]} for t in (300.0, 330.0, 360.0)]}
        if package == "UNIQUAC":
            recs = uniquac_records(pk)
            r = [pk.comp(i).UNIQUAC_R for i in range(2)]
            q = [pk.comp(i).UNIQUAC_Q for i in range(2)]
            first = {k: v for k, v in recs.items() if k == "Methanol|Water"}
            entry["records"] = recs
            for g in entry["gamma"]:
                g["formula"] = uniquac_dwsim(g["x"], g["T"], r, q, recs, mw)
                g["formula_one_record"] = uniquac_dwsim(g["x"], g["T"], r, q, first, mw)
        if package == "UNIFAC-LL":
            ll = unifac_ll_data(pk)
            entry["groups"] = ll
            r = [c["R"] for c in ll["compounds"]]
            q = [c["Q"] for c in ll["compounds"]]
            e = [c["e"] for c in ll["compounds"]]
            for g in entry["gamma"]:
                g["formula"] = unifac_ll_dwsim(g["x"], g["T"], r, q, e, ll["a"])
                g["formula_e1"] = unifac_ll_dwsim(g["x"], g["T"], r, q, [1.0] * len(e), ll["a"])
        if package == "Wilson":
            bips = wilson_bips(pk)
            v = [pk.comp(i).Molar_Weight / pk.rho(i, 298.15) for i in range(2)]
            entry["bips"] = bips
            entry["volumes_298"] = v
            for g in entry["gamma"]:
                g["formula"] = wilson_dwsim(g["x"], g["T"], bips, v, mw)
        out["packages"][package] = entry
        worst = max(abs(math.log(a / b)) for g in entry["gamma"] if "formula" in g for a, b in zip(g["gamma"], g["formula"])) \
            if "formula" in entry["gamma"][0] else None
        print(package, "gamma at 330 K, x1 0.3:", entry["gamma"][6]["gamma"], "formula gap", worst)
    # NRTL PT flashes of methanol-water, two-phase at 1 atm.
    out["flashes"] = []
    uniquac = Package(mw, "UNIQUAC", [0.5, 0.5])
    for package, pk in [("NRTL", nrtl), ("UNIQUAC", uniquac)]:
        for t, z in [(350.0, [0.5, 0.5]), (357.0, [0.3, 0.7]), (346.0, [0.7, 0.3])]:
            f = pk.flash(t, z)
            out["flashes"].append({"package": package, "T": t, "z": z, **f, "gamma_z": pk.gamma(z, t),
                                   "psat": [pk.psat(i, t) for i in range(2)],
                                   "rho": [pk.rho(i, t) for i in range(2)]})
            print(package, "flash", t, z, f["beta"], f["y"], f["x1"])
    # The three-phase flash of water-1-butanol by NRTL, by each approach.
    wb = ["Water", BUTANOL]
    out["three_phase"] = {"compounds": wb, "runs": []}
    for approach in ["NestedLoops", "InsideOut", "GibbsMinimization"]:
        for stab in [False, True]:
            pk = Package(wb, "NRTL", [0.5, 0.5])
            set_enum(pk.pp, "FlashCalculationApproach", approach)
            if stab:
                pk.pp.FlashSettings[FlashSetting.ThreePhaseFlashStabTestSeverity] = "2"
                pk.pp.FlashSettings[FlashSetting.ThreePhaseFlashStabTestCompIds] = ",".join(wb)
            pk.pp.ForceNewFlashAlgorithmInstance = True
            # Water 0.5 is outside ChemSep's NRTL split (0.6008 to 0.9944
            # water at 298.15 K), so one liquid is right there; 0.8 is inside.
            for zw in (0.5, 0.8):
                for t in (298.15, 330.0, 360.0):
                    try:
                        f = pk.flash(t, [zw, 1.0 - zw])
                    except Exception as e:  # noqa: BLE001
                        f = {"error": str(e).splitlines()[0][:200]}
                    out["three_phase"]["runs"].append({"approach": approach, "stability_test": stab, "T": t,
                                                       "z": [zw, 1.0 - zw],
                                                       **figures({k: v for k, v in f.items() if k != "iterations"},
                                                                 10 if approach == "NestedLoops" else 3)})
                    print("three phases", approach, stab, zw, t, {k: v for k, v in f.items() if k != "iterations"})
    # A pure liquid's enthalpy less its ideal gas's, J/mol, by each mode: its
    # H as a liquid and as a vapour, the ideal gas, at 330 K and 1 atm.
    out["enthalpy"] = []
    for mode in ["LeeKesler", "Ideal", "Excess", "ExpData", "LiqCp_Excess"]:
        pk = Package(mw, "NRTL", [0.5, 0.5])
        set_enum(pk.pp, "EnthalpyEntropyCpCvCalculationMode", mode)
        row = {"mode": mode, "T": 330.0, "liquid_less_gas": [],
               "hvap": [pk.pp.RET_HVAPM(vec([1.0 if j == i else 0.0 for j in range(2)]), 330.0)
                        * pk.comp(i).Molar_Weight for i in range(2)]}
        for i in range(2):
            pure = [1.0 if j == i else 0.0 for j in range(2)]
            m = pk.comp(i).Molar_Weight
            try:
                hl = pk.pp.DW_CalcEnthalpy(vec(pure), 330.0, P, State.Liquid) * m
                hv = pk.pp.DW_CalcEnthalpy(vec(pure), 330.0, P, State.Vapor) * m
                row["liquid_less_gas"].append(hl - hv)
            except Exception as e:  # noqa: BLE001
                row["liquid_less_gas"].append(None)
                row["error"] = str(e).splitlines()[0][:200]
        out["enthalpy"].append(row)
        print("enthalpy", mode, row["liquid_less_gas"])
    (ROOT / "testdata" / "m6_dwsim.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8", newline="\n")
    sys.stdout.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
