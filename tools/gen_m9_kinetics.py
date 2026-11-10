# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["mpmath==1.3.0"]
# ///
"""Generate M9's fixture for rate laws (step 5).

Run from the repository root:

    uv run --python 3.13 tools/gen_m9_kinetics.py

It writes testdata/m9_kinetics.json: two published kinetic models, each
rate evaluated by mpmath at 50 digits from the paper's own formulas,
written here independently of reflux, at states spread over each model's
range, with each rate's derivative with respect to T and to each partial
pressure (partial pressures held fixed for the one, T for the others).

Xu and Froment (AIChE J. 35, 88, 1989): methane steam reforming, the
water-gas shift and their sum on Ni/MgAl2O4, three LHHW rates over one
adsorption term, squared (their Eq. 3). Rates in kmol/(kg·h), pressures in
bar. Their paper names the equilibrium constants K1, K2 and K3 and gives
their units, but no values: the correlations here are the ones later
papers use with these kinetics.

DWSIM: the owner's reforming case's rates as DWSIM has them typed, from
spike/dwsim_pfr/run.py: Xu and Froment's constants rounded, R typed as
8.314, partial pressures in atm, each reaction's own denominator and its own
equilibrium constants. Rates in kmol/(kg·h).

Vanden Bussche and Froment (J. Catal. 161, 1, 1996): methanol from CO2
and the reverse water-gas shift on Cu/ZnO/Al2O3, over one site's term,
cubed for methanol (their Eq. 3). Rates in mol/(kg·s), pressures in bar,
and the equilibrium constants Graaf, Stamhuis and Beenackers' (Chem. Eng.
Sci. 41, 2883, 1986), as Vanden Bussche and Froment print them.

The constants are below, each with the table it comes from. Each paper
gives its constants twice, centred on reference temperatures and as
pre-exponential factors, and this script checks that the two agree within
their rounding before it writes anything. The rates are written in SI
here, mol/(kg·s) and Pa, as reflux's are.

The output is the same, byte for byte, for the same versions.
"""

import json
import pathlib

import mpmath as mp

mp.mp.dps = 50

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m9_kinetics.json"

R = mp.mpf("8.31446261815324")
BAR = mp.mpf(100000)
KMOL_PER_H = mp.mpf(1000) / 3600  # kmol/h in mol/s


def M(s):
    return mp.mpf(s)


# Xu and Froment's constants: A from their Table 6 ("preexponential factors
# for final model", at the catalyst's reference activity: Table 5's text
# multiplies the rate coefficients by 2.246 for fresh catalyst), E and ΔH
# from Table 5. Rate constants k = A·exp(−E/RT), kmol·bar^0.5/(kg·h) for k1
# and k3 and kmol/(kg·h·bar) for k2; adsorption constants K = A·exp(−ΔH/RT),
# 1/bar, K_H2O's a pure number. The equilibrium constants are not the
# paper's: later papers use ln K1 = 30.114 − 26830/T (bar²) and ln K2 =
# −4.036 + 4400/T with these kinetics, with K3 = K1·K2; one prints 30.11 and
# −4.063.
XF = {
    "k1": (M("4.225e15"), M("240.1e3")),
    "k2": (M("1.955e6"), M("67.13e3")),
    "k3": (M("1.020e15"), M("243.9e3")),
    "K_CO": (M("8.23e-5"), M("-70.65e3")),
    "K_H2": (M("6.12e-9"), M("-82.90e3")),
    "K_CH4": (M("6.65e-4"), M("-38.28e3")),
    "K_H2O": (M("1.77e5"), M("88.68e3")),
    "K1": (M("30.114"), M("-26830")),  # ln K1 = a + b/T, bar²
    "K2": (M("-4.036"), M("4400")),  # ln K2 = a + b/T
}

# Xu and Froment's Table 5: the same constants centred, k(T) = k_Tr·exp(−E/R·
# (1/T − 1/Tr)), with Tr = 648 K for the rate constants, K_CO and K_H2, and
# 823 K for K_CH4 and K_H2O.
XF_CENTRED = {
    "k1": (M("1.842e-4"), 648),
    "k2": (M("7.558"), 648),
    "k3": (M("2.193e-5"), 648),
    "K_CO": (M("40.91"), 648),
    "K_H2": (M("0.02960"), 648),
    "K_CH4": (M("0.1791"), 823),
    "K_H2O": (M("0.4152"), 823),
}


def xf_rates(t, p):
    """Xu and Froment's three rates, mol/(kg·s), at T in K and partial
    pressures in Pa (CH4, H2O, H2, CO, CO2)."""
    pb = {k: v / BAR for k, v in p.items()}  # bar

    def arr(name):
        a, e = XF[name]
        return a * mp.exp(-e / (R * t))

    k1, k2, k3 = arr("k1"), arr("k2"), arr("k3")
    kco, kh2, kch4, kh2o = arr("K_CO"), arr("K_H2"), arr("K_CH4"), arr("K_H2O")
    K1 = mp.exp(XF["K1"][0] + XF["K1"][1] / t)
    K2 = mp.exp(XF["K2"][0] + XF["K2"][1] / t)
    K3 = K1 * K2
    den = 1 + kco * pb["CO"] + kh2 * pb["H2"] + kch4 * pb["CH4"] + kh2o * pb["H2O"] / pb["H2"]
    r1 = k1 / pb["H2"] ** M("2.5") * (pb["CH4"] * pb["H2O"] - pb["H2"] ** 3 * pb["CO"] / K1) / den**2
    r2 = k2 / pb["H2"] * (pb["CO"] * pb["H2O"] - pb["H2"] * pb["CO2"] / K2) / den**2
    r3 = k3 / pb["H2"] ** M("3.5") * (pb["CH4"] * pb["H2O"] ** 2 - pb["H2"] ** 4 * pb["CO2"] / K3) / den**2
    return [r * KMOL_PER_H for r in (r1, r2, r3)]


# Vanden Bussche and Froment's constants, from their Table 2, each
# A·exp(B/RT): the methanol rate's lumped constant k'5a·K'2·K3·K4·K_H2 in
# mol/(kg·s·bar²); the RWGS rate's k'1 in mol/(kg·s·bar); K_H2O/(K8·K9·K_H2),
# a pure number; √K_H2 in bar^−0.5; K_H2O in 1/bar. Graaf's equilibrium
# constants as their text prints them, log10 K1* = 3066/T − 10.592 for
# methanol from CO2 (bar^−2, the pressures being in bar) and
# log10(1/K3*) = −2073/T + 2.029, 1/K3* being the RWGS's K.
VBF = {
    "kM": (M("1.07"), M("36696")),
    "kR": (M("1.22e10"), M("-94765")),
    "Ka": (M("3453.38"), M("0")),
    "Kb": (M("0.499"), M("17197")),
    "Kc": (M("6.62e-11"), M("124119")),
    "KM": (M("-10.592"), M("3066")),
    "KR": (M("2.029"), M("-2073")),
}

# Vanden Bussche and Froment's Table 1: the same constants as fitted,
# κ = A*·exp(−B/R·(1/T_av − 1/T)), T_av = 501.57 K.
VBF_CENTRED = {
    "kM": M("7070.34"),
    "kR": M("1.65"),
    "Ka": M("3453.38"),
    "Kb": M("30.82"),
    "Kc": M("558.17"),
}

# How far a paper's two forms could disagree from their rounding alone: half
# a unit in the last figure of each of the A, the centred value and the
# energy (as a factor exp(δE/(R·Tr)) − 1). Vanden Bussche and Froment's
# agree within it. Xu and Froment's don't, for k1, k2, k3, K_CO and K_H2:
# their Table 6 puts k1 and k3 1.55% and 1.72% above Table 5 at 648 K, as if
# from activation energies of 240.18 and 243.99 kJ/mol rather than the 240.1
# and 243.9 printed, and no other R reconciles them. Each pair's ratio is
# recorded beside its rounding, and anything over 2%, which a misread figure
# would be, stops the script.


def half_unit(s):
    """Half a unit in the last written figure of s, relative to s."""
    mant = s.lower().split("e")[0].replace("-", "")
    places = len(mant.split(".")[1]) if "." in mant else 0
    return (M(1) / 2) * M(10) ** (-places) / M(mant)


def check_two_forms():
    """Each paper's centred table against its pre-exponential one, at each
    reference temperature: their ratio, and how far their rounding alone
    could take it from 1."""
    out = {}
    a_text = {"k1": "4.225", "k2": "1.955", "k3": "1.020", "K_CO": "8.23", "K_H2": "6.12", "K_CH4": "6.65",
              "K_H2O": "1.77"}
    c_text = {"k1": "1.842", "k2": "7.558", "k3": "2.193", "K_CO": "40.91", "K_H2": "0.02960", "K_CH4": "0.1791",
              "K_H2O": "0.4152"}
    e_text = {"k1": "240.1", "k2": "67.13", "k3": "243.9", "K_CO": "70.65", "K_H2": "82.90", "K_CH4": "38.28",
              "K_H2O": "88.68"}
    for name, (centred, tr) in XF_CENTRED.items():
        a, e = XF[name]
        from_a = a * mp.exp(-e / (R * tr))
        e_unit = M(10) ** (-len(e_text[name].split(".")[1])) / 2 * 1000
        rounding = (mp.exp(e_unit / (R * tr)) - 1) + half_unit(a_text[name]) + half_unit(c_text[name])
        ratio = from_a / centred
        if abs(ratio - 1) > M("0.02"):
            raise SystemExit(f"Xu and Froment's {name}: Table 6 gives {from_a} at {tr} K, Table 5 {centred}")
        out[f"Xu and Froment {name}"] = {"ratio": float(ratio), "rounding": float(rounding)}
    vbf_a = {"kM": "1.07", "kR": "1.22", "Ka": "3453.38", "Kb": "0.499", "Kc": "6.62"}
    vbf_c = {"kM": "7070.34", "kR": "1.65", "Ka": "3453.38", "Kb": "30.82", "Kc": "558.17"}
    tav = M("501.57")
    for name, centred in VBF_CENTRED.items():
        a, b = VBF[name]
        from_a = a * mp.exp(b / (R * tav))
        rounding = half_unit(vbf_a[name]) + half_unit(vbf_c[name]) + (mp.exp(M(1) / (R * tav)) - 1)
        ratio = from_a / centred
        if abs(ratio - 1) > rounding:
            raise SystemExit(f"Vanden Bussche and Froment's {name}: Table 2 gives {from_a}, Table 1 {centred}")
        out[f"Vanden Bussche and Froment {name}"] = {"ratio": float(ratio), "rounding": float(rounding)}
    return out

def vbf_rates(t, p):
    """Vanden Bussche and Froment's two rates, mol/(kg·s), at T in K and
    partial pressures in Pa (CO2, H2, CH3OH, H2O, CO)."""
    pb = {k: v / BAR for k, v in p.items()}

    def ab(name):
        a, b = VBF[name]
        return a * mp.exp(b / (R * t))

    km, kr, ka, kb, kc = ab("kM"), ab("kR"), ab("Ka"), ab("Kb"), ab("Kc")
    KM = mp.power(10, VBF["KM"][0] + VBF["KM"][1] / t)
    KR = mp.power(10, VBF["KR"][0] + VBF["KR"][1] / t)
    theta = 1 + ka * pb["H2O"] / pb["H2"] + kb * mp.sqrt(pb["H2"]) + kc * pb["H2O"]
    rm = km * pb["CO2"] * pb["H2"] * (1 - pb["H2O"] * pb["CH3OH"] / (pb["H2"] ** 3 * pb["CO2"]) / KM) / theta**3
    rr = kr * pb["CO2"] * (1 - pb["H2O"] * pb["CO"] / (pb["CO2"] * pb["H2"]) / KR) / theta
    return [rm, rr]


def states(species, ts, ps, fractions):
    out = []
    for t in ts:
        for ptot in ps:
            for y in fractions:
                out.append((M(t), {s: M(ptot) * M(f) for s, f in zip(species, y)}))
    return out


def model(name, species, rates, sts):
    cases = []
    for t, p in sts:
        r = rates(t, p)
        drdt = [mp.diff(lambda x, i=i: rates(x, p)[i], t) for i in range(len(r))]
        drdp = []
        for s in species:
            def at(x, s=s, i=0):
                q = dict(p)
                q[s] = x
                return rates(t, q)
            drdp.append([mp.diff(lambda x, s=s, i=i: at(x, s)[i], p[s]) for i in range(len(r))])
        cases.append({
            "t": float(t),
            "p": [float(p[s]) for s in species],
            "rates": [float(x) for x in r],
            "drdt": [float(x) for x in drdt],
            "drdp": [[float(x) for x in row] for row in drdp],
            "scales": [float(max(abs(x) for x in terms)) for terms in scale_terms(name, t, p)],
        })
    return {"name": name, "species": species, "cases": cases}


ATM = mp.mpf(101325)


def dwsim_rates(t, p, forward_only=False):
    """The owner's DWSIM case's three rates as DWSIM has them typed
    (spike/dwsim_pfr/run.py): partial pressures in atm, R typed as 8.314,
    rates in kmol/(kg·h), each reaction's own denominator, and its own
    equilibrium constants. A reaction's R1, R2, P1 and P2 are its reactants
    and products in the order of its compound list. Returned as steam
    reforming, the shift and the overall reaction, in mol/(kg·s)."""
    a = {k: v / ATM for k, v in p.items()}
    e = lambda x: mp.exp(x)
    r8 = M("8.314")
    big = M("1e300") if forward_only else M(0)
    k_sr = e(M("30.42") - M(27106) / t + big)
    k_wgs = e(M("-3.798") + M(4160) / t + big)
    k_ov = e(M("34.218") - M(31266) / t + big)
    # Steam reforming: R1 CH4, R2 H2O, P1 H2, P2 CO.
    num = M("4.22e15") * e(M(-240100) / r8 / t) / a["H2"] ** M("2.5") * (a["CH4"] * a["H2O"] - a["H2"] ** 3 * a["CO"] / k_sr)
    den = (1 + M("8.23e-5") * e(M(70650) / r8 / t) * a["CO"] + M("6.65e-4") * e(M(38280) / r8 / t) * a["CH4"]
           + M("1.77e5") * e(M(-88680) / r8 / t) * a["H2O"] / a["H2"] + M("6.12e-9") * e(M(82900) / r8 / t) * a["H2"]) ** 2
    sr = num / den
    # The shift: R1 H2O, R2 CO, P1 H2, P2 CO2.
    num = M("1.96e6") * e(M(-67130) / r8 / t) / a["H2"] * (a["H2O"] * a["CO"] - a["H2"] * a["CO2"] / k_wgs)
    den = (1 + M("1.77e5") * e(M(-88680) / r8 / t) * a["H2O"] / a["H2"] + M("6.12e-9") * e(M(82900) / r8 / t) * a["H2"]
           + M("8.23e-5") * e(M(70650) / r8 / t) * a["CO"]) ** 2
    wgs = num / den
    # The overall reaction: R1 CH4, R2 H2O, P1 H2, P2 CO2.
    num = M("1.02e15") * e(M(-243900) / r8 / t) / a["H2"] ** M("3.5") * (a["CH4"] * a["H2O"] ** 2 - a["H2"] ** 4 * a["CO2"] / k_ov)
    den = (1 + M("6.65e-4") * e(M(38280) / r8 / t) * a["CH4"] + M("1.77e5") * e(M(-88680) / r8 / t) * a["H2O"] / a["H2"]
           + M("6.12e-9") * e(M(82900) / r8 / t) * a["H2"]) ** 2
    ov = num / den
    return [r * KMOL_PER_H for r in (sr, wgs, ov)]


def scale_terms(name, t, p):
    """Each rate's forward and backward terms, the magnitudes its rounding is
    measured against where the two nearly cancel."""
    if name == "DWSIM":
        return [[a, b] for a, b in zip(dwsim_rates(t, p), dwsim_rates(t, p, True))]
    if name == "Xu and Froment":
        sp = dict(p)
        fwd = xf_rates(t, sp)
        # The forward term alone: the rate with the backward one dropped.
        XF2 = dict(XF)
        XF["K1"] = (M("1e300"), M(0))
        XF["K2"] = (M("1e300"), M(0))
        forward = xf_rates(t, sp)
        XF.update(XF2)
        return [[a, b] for a, b in zip(fwd, forward)]
    fwd = vbf_rates(t, p)
    VBF2 = dict(VBF)
    VBF["KM"] = (M("1e300"), M(0))
    VBF["KR"] = (M("1e300"), M(0))
    forward = vbf_rates(t, p)
    VBF.update(VBF2)
    return [[a, b] for a, b in zip(fwd, forward)]


def main():
    agreement = check_two_forms()
    xf_species = ["CH4", "H2O", "H2", "CO", "CO2"]
    xf_states = states(xf_species, ["700", "800", "900", "1000", "1100"], ["1e5", "15e5", "30e5"], [
        ["0.30", "0.60", "0.02", "0.05", "0.03"],
        ["0.20", "0.45", "0.20", "0.10", "0.05"],
        ["0.05", "0.30", "0.45", "0.12", "0.08"],
    ])
    vbf_species = ["CO2", "H2", "CH3OH", "H2O", "CO"]
    vbf_states = states(vbf_species, ["473.15", "503.15", "533.15", "563.15"], ["50e5", "80e5"], [
        ["0.10", "0.70", "0.01", "0.01", "0.05"],
        ["0.08", "0.65", "0.04", "0.03", "0.08"],
        ["0.20", "0.60", "0.002", "0.004", "0.01"],
    ])
    data = {
        "source": {
            "by": f"tools/gen_m9_kinetics.py, with mpmath {mp.__version__} at 50 digits",
            "models": "each rate from its paper's own formulas, written independently of reflux; the constants "
                      "are in the script with where each comes from",
            "constants": "Xu and Froment's A from their Table 6 and E and dH from Table 5, at the catalyst's "
                         "reference activity; their K1 and K2 the correlations later papers use, since theirs gives "
                         "none. Vanden Bussche and Froment's from their Table 2, and Graaf's equilibrium constants "
                         "as they print them. Each paper's centred table agrees with its pre-exponential one within "
                         "their rounding but Xu and Froment's k1, k2, k3, K_CO and K_H2: agreement is each pair's ratio, "
                         "pre-exponential over centred at the reference temperature, and how far their rounding "
                         "alone could take it from 1",
            "units": "T in K, partial pressures in Pa in each model's species order, rates in mol/(kg·s); drdt "
                     "is each rate's derivative with respect to T at fixed partial pressures, drdp[j][i] rate "
                     "i's with respect to species j's partial pressure; scales are each rate's largest term, "
                     "its forward one or itself",
        },
        "agreement": agreement,
        "models": [
            model("Xu and Froment", xf_species, xf_rates, xf_states),
            model("Vanden Bussche and Froment", vbf_species, vbf_rates, vbf_states),
            model("DWSIM", xf_species, dwsim_rates, xf_states),
        ],
    }
    OUT.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}: {sum(len(m['cases']) for m in data['models'])} states")


main()
