# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "iapws==1.5.4"]
# ///
"""M5a plan spike: DWSIM 9.0.5's packages for pure water, driven from Python
with no GUI, against IAPWS-IF97 and IAPWS-95 as the iapws package has them.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/water/dwsim.py

It needs DWSIM 9 in %LOCALAPPDATA%\\DWSIM. It:
  1. lists the property packages whose names suggest water or CoolProp;
  2. for each of DWSIM's "Steam Tables (IAPWS-IF97)" and "CoolProp" packages,
     flashes pure water at (T, P) states in each IF97 region, by PH and PS
     back, and at a vapour fraction, and compares h, s, v and the phase with
     IF97's and IAPWS-95's;
  3. times its PT and PH flashes.

DWSIM reports specific enthalpy in kJ/kg and entropy in kJ/(kg K).
"""
import math
import os
import sys
import time

from iapws import IAPWS95, IAPWS97

# (T K, P Pa): IF97's verification states by region, and a few more.
STATES = [
    (300.0, 3e6, "region 1, Table 5"), (300.0, 80e6, "region 1, Table 5"),
    (500.0, 3e6, "region 1, Table 5"),
    (300.0, 3.5e3, "region 2, Table 15"), (700.0, 3.5e3, "region 2, Table 15"),
    (700.0, 30e6, "region 2, Table 15"),
    (650.0, 25.5837018e6, "region 3, Table 33 (500 kg/m3)"),
    (650.0, 22.2930643e6, "region 3, Table 33 (200 kg/m3)"),
    (750.0, 78.3095639e6, "region 3, Table 33 (500 kg/m3)"),
    (1500.0, 0.5e6, "region 5, Table 42"), (1500.0, 30e6, "region 5, Table 42"),
    (2000.0, 30e6, "region 5, Table 42"),
    (373.15, 101325.0, "1 atm, 100 degC"), (647.0, 22.0e6, "near critical"),
]


def load_dwsim():
    import pythonnet

    pythonnet.load("netfx")
    import clr

    home = os.path.join(os.environ["LOCALAPPDATA"], "DWSIM")
    sys.path.append(home)
    os.chdir(home)
    for name in ["CapeOpen", "DWSIM.Automation", "DWSIM.Interfaces", "DWSIM.GlobalSettings",
                 "DWSIM.SharedClasses", "DWSIM.Thermodynamics", "DWSIM.UnitOperations",
                 "DWSIM.Inspector", "System.Buffers"]:
        clr.AddReference(name)


def ref95(t, p):
    w = IAPWS95(T=t, P=p / 1e6)
    return dict(h=w.h, s=w.s, v=1 / w.rho, phase=w.phase)


def ref97(t, p):
    w = IAPWS97(T=t, P=p / 1e6)
    return dict(h=w.h, s=w.s, v=w.v, region=w.region)


def main():
    load_dwsim()
    from System import Array, Double
    from System.Diagnostics import FileVersionInfo
    from DWSIM.Automation import Automation3
    from DWSIM.Interfaces.Enums import FlashCalculationType
    from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType

    auto = Automation3()
    print(auto.GetVersion())
    home = os.environ["LOCALAPPDATA"] + "\\DWSIM\\"
    for dll in ["CoolProp.dll", "DWSIM.Thermodynamics.CoolPropInterface.dll"]:
        fv = FileVersionInfo.GetVersionInfo(home + dll)
        print(f"  {dll}: file version {fv.FileVersion}, product {fv.ProductVersion}")
    sim = auto.CreateFlowsheet()
    names = [str(k) for k in sim.AvailablePropertyPackages.Keys]
    print("property packages:", len(names))
    for n in names:
        if any(w in n.lower() for w in ["steam", "iapws", "water", "coolprop", "seawater"]):
            print("  ", n)
    sim.SelectedCompounds.Add("Water", sim.AvailableCompounds["Water"])
    mw = sim.SelectedCompounds["Water"].Molar_Weight
    print(f"  Water's molar mass in DWSIM: {mw} g/mol")
    vec = lambda xs: Array[Double](list(xs))
    s = sim.AddObject(ObjectType.MaterialStream, 50, 50, "w").GetAsObject()
    s.SetTemperature(300.0)
    s.SetPressure(1e5)
    s.SetMolarFlow(1.0)
    s.SetOverallComposition(vec([1.0]))

    for ppname in [n for n in names if "IF97" in n or n.strip().lower() == "coolprop"
                   or n.startswith("CoolProp (")]:
        if "Incompressible" in ppname:
            continue
        print(f"\n=== {ppname}")
        try:
            pp = sim.CreateAndAddPropertyPackage(ppname).__implementation__
        except Exception as e:  # noqa: BLE001
            print("  can't make it:", e)
            continue
        pp.CurrentMaterialStream = s
        print(f"  class {pp.GetType().FullName}")

        def pt(t, p):
            r = pp.CalculateEquilibrium(FlashCalculationType.PressureTemperature, p, t, vec([1.0]),
                                        None, 0.0)
            return r

        offs_h = offs_s = None
        for t, p, what in STATES:
            try:
                r = pt(t, p)
                h, sv = r.CalculatedEnthalpy, r.CalculatedEntropy
                beta = r.GetVaporPhaseMoleFraction()
            except Exception as e:  # noqa: BLE001
                print(f"  {what:32} T {t} K P {p:.6g} Pa: FAILED {str(e).splitlines()[0][:120]}")
                continue
            a95, a97 = ref95(t, p), ref97(t, p)
            if offs_h is None:
                offs_h, offs_s = h - a95["h"], sv - a95["s"]
            print(f"  {what:32} T {t:7.2f} K P {p:12.6g} Pa: beta {beta:.3g} "
                  f"h {h:12.6f} kJ/kg (IF97 {a97['h']:12.6f}, d {h - a97['h']: .2e}; "
                  f"IAPWS-95 {a95['h']:12.6f}, d {h - a95['h']: .2e}) "
                  f"s d97 {sv - a97['s']: .2e} d95 {sv - a95['s']: .2e}")

        # PH and PS back, and a two-phase PH.
        for t, p in [(500.0, 3e6), (700.0, 30e6), (650.0, 25.5837018e6)]:
            h = pt(t, p).CalculatedEnthalpy
            try:
                r = pp.CalculateEquilibrium(FlashCalculationType.PressureEnthalpy, p, h, vec([1.0]),
                                            None, 0.0)
                print(f"  PH back at T {t} K P {p:.6g}: T {r.CalculatedTemperature!r} "
                      f"(off {r.CalculatedTemperature - t:.2e} K)")
            except Exception as e:  # noqa: BLE001
                print(f"  PH back at T {t} K P {p:.6g}: FAILED {str(e).splitlines()[0][:120]}")
        p = 1e5
        w = IAPWS97(P=0.1, x=0.5)
        try:
            r = pp.CalculateEquilibrium(FlashCalculationType.PressureEnthalpy, p, w.h, vec([1.0]),
                                        None, 0.0)
            print(f"  PH in the dome, 1 bar, IF97's h at x = 0.5 ({w.h:.6f} kJ/kg): "
                  f"T {r.CalculatedTemperature!r} (IF97 Ts {w.T!r}), beta "
                  f"{r.GetVaporPhaseMoleFraction()!r}")
        except Exception as e:  # noqa: BLE001
            print(f"  PH in the dome: FAILED {str(e).splitlines()[0][:160]}")
        try:
            r = pp.CalculateEquilibrium(FlashCalculationType.PressureVaporFraction, p, 0.0,
                                        vec([1.0]), None, 0.0)
            print(f"  bubble T at 1 bar: {r.CalculatedTemperature!r} K (IF97 {IAPWS97(P=0.1, x=0).T!r}, "
                  f"IAPWS-95 {IAPWS95(P=0.1, x=0).T!r})")
        except Exception as e:  # noqa: BLE001
            print(f"  bubble T: FAILED {str(e).splitlines()[0][:160]}")

        # Timings.
        n = 300
        t0 = time.perf_counter()
        for i in range(n):
            pt(400.0 + 0.5 * (i % 50), 5e6)
        t_pt = (time.perf_counter() - t0) / n
        h0 = pt(600.0, 5e6).CalculatedEnthalpy
        t0 = time.perf_counter()
        for i in range(n):
            pp.CalculateEquilibrium(FlashCalculationType.PressureEnthalpy, 5e6, h0 * (1 + 1e-3 * (i % 50)),
                                    vec([1.0]), None, 0.0)
        t_ph = (time.perf_counter() - t0) / n
        print(f"  timings: PT {t_pt * 1e3:.3f} ms, PH {t_ph * 1e3:.3f} ms ({n} each)")


if __name__ == "__main__":
    main()
