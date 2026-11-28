# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["chemicals==1.5.2"]
# ///
"""M5b step 0: ChemSep's water transport correlations (reflux's component
file, data/components/7732-18-5.toml) against IAPWS's (R12-08's viscosity
and R15-11's conductivity without their critical enhancements, R1-76's
surface tension, all by chemicals), on IAPWS-95's saturated liquid and on
the vapour at 1 kPa. It answers whether a mixture's water should take
ChemSep's correlations, as every other compound does, or IAPWS's.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/transport/water.py

Not part of the simulator.
"""
import math
import pathlib
import tomllib

from chemicals.dippr import EQ101, EQ102
from chemicals.iapws import iapws95_rho, iapws95_rhol_sat
from chemicals.interface import sigma_IAPWS
from chemicals.thermal_conductivity import k_IAPWS
from chemicals.viscosity import mu_IAPWS

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
WATER = tomllib.loads((ROOT / "data/components/7732-18-5.toml").read_text(encoding="utf-8"))["correlations"]


def eq16(t, a, b, c, d, e):
    return a + math.exp(b / t + c + d * t + e * t * t)


def chemsep(key, t):
    d = WATER[key]
    k = (d["coeffs"] + [0.0] * 5)[:5]
    return {101: lambda: EQ101(t, *k), 102: lambda: EQ102(t, *k[:4]), 16: lambda: eq16(t, *k)}[d["eq"]]()


def main():
    print("saturated liquid: T, ChemSep/IAPWS - 1 for mu, k, sigma")
    for t in [275.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 630.0]:
        rho = iapws95_rhol_sat(t)
        r = (chemsep("liquid_viscosity", t) / mu_IAPWS(t, rho) - 1,
             chemsep("liquid_thermal_conductivity", t) / k_IAPWS(t, rho) - 1,
             chemsep("surface_tension", t) / sigma_IAPWS(t) - 1)
        print(f"  {t:6.1f} K  mu {r[0]:+.4f}  k {r[1]:+.4f}  sigma {r[2]:+.4f}")
    print("vapour at 1 kPa: T, ChemSep/IAPWS - 1 for mu, k")
    for t in [300.0, 400.0, 600.0, 800.0, 1000.0, 1073.15]:
        rho = iapws95_rho(t, 1e3)
        r = (chemsep("vapor_viscosity", t) / mu_IAPWS(t, rho) - 1,
             chemsep("vapor_thermal_conductivity", t) / k_IAPWS(t, rho) - 1)
        print(f"  {t:6.1f} K  mu {r[0]:+.4f}  k {r[1]:+.4f}")


if __name__ == "__main__":
    main()
