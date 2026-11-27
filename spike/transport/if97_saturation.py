# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["iapws==1.5.4", "numpy==2.2.6", "scipy==1.14.1"]
# ///
"""M5b step 0: IF97 on region 4's saturation line. At ps(T) from region 4,
the Gibbs energies of region 1's liquid and region 2's vapour, and the T at
that pressure where they're equal: what the flash's fugacity equality would
take for saturation if IF97 went behind it as a phase model with nothing
more. The two lines differ by IF97's own inconsistency (R7-97's Table 44).

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/transport/if97_saturation.py

Not part of the simulator.
"""
from iapws.iapws97 import _PSat_T, _Region1, _Region2
from scipy.optimize import brentq

R = 0.461526  # kJ/(kg K), IF97's


def g(region, t, p):
    r = region(t, p)
    return r["h"] - t * r["s"]


def main():
    worst = [0.0, 0.0, 0.0]
    for t in [273.16, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 620.0, 623.15]:
        p = _PSat_T(t)
        dg = (g(_Region1, t, p) - g(_Region2, t, p)) / (R * t)
        teq = brentq(lambda x: g(_Region1, x, p) - g(_Region2, x, p), t - 1, t + 1, xtol=1e-12)
        dh = _Region1(teq, p)["h"] - _Region1(t, p)["h"]
        worst = [max(worst[0], abs(dg)), max(worst[1], abs(teq - t)), max(worst[2], abs(dh))]
        print(f"T {t:7.2f} K  ps {p:.6g} MPa  (g1 - g2)/RT {dg:+.3e}  equal-G T - Ts {1e3 * (teq - t):+.4f} mK  "
              f"liquid h moves {1e3 * dh:+.3f} J/kg")
    print(f"worst: (g1 - g2)/RT {worst[0]:.2e}, T {1e3 * worst[1]:.3f} mK, h {1e3 * worst[2]:.3f} J/kg")


if __name__ == "__main__":
    main()
