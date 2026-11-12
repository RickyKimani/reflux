# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["chemicals==1.5.2"]
# ///
"""Generate M9's fixture for the gas viscosity Ergun needs (step 9, D10 as
amended in step 0).

Run from the repository root:

    uv run --python 3.13 tools/gen_m9_viscosity.py

It writes testdata/m9_viscosity.json: gas viscosities at low pressure, each
compound's from its component file's vapour-viscosity correlation (ChemSep's
equation 102, evaluated by chemicals' EQ102), and a mixture's by Wilke's rule
(chemicals' Wilke, with the component files' molar masses). Pure gases from
300 to 1000 K, the owner's reforming gas at its feed and two points along
the reactor, a methanol synthesis gas, air, and hydrogen with carbon dioxide
across their compositions, where Wilke's rule and a mole-fraction mean
differ most.

The output is the same, byte for byte, for the same versions.
"""

import json
import pathlib
import tomllib

import chemicals
from chemicals.dippr import EQ102
from chemicals.viscosity import Wilke

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "components"
OUT = ROOT / "testdata" / "m9_viscosity.json"

CH4, H2O, H2, CO, CO2 = "74-82-8", "7732-18-5", "1333-74-0", "630-08-0", "124-38-9"
N2, O2, MEOH, NH3 = "7727-37-9", "7782-44-7", "67-56-1", "7664-41-7"


def record(cas):
    return tomllib.loads((DATA / f"{cas}.toml").read_text(encoding="utf-8"))


def pure(cas, t):
    c = record(cas)["correlations"]["vapor_viscosity"]
    assert c["eq"] == 102 and c["units"] == "Pa.s", (cas, c)
    a, b, cc, d = (c["coeffs"] + [0.0, 0.0, 0.0, 0.0])[:4]
    return EQ102(t, a, b, cc, d)


def mw(cas):
    return record(cas)["constants"]["mw"] * 1000  # g/mol, as chemicals takes it


def case(name, cas, y, t):
    total = sum(y)
    y = [v / total for v in y]
    mus = [pure(c, t) for c in cas]
    return {
        "name": name,
        "cas": cas,
        "y": y,
        "t": t,
        "pure": mus,
        "mixture": Wilke(y, mus, [mw(c) for c in cas]) if len(cas) > 1 else mus[0],
    }


def main():
    cases = []
    for cas in [CH4, H2O, H2, CO, CO2, N2, O2, MEOH, NH3]:
        for t in [300.0, 600.0, 1000.0]:
            cases.append(case(f"pure {record(cas)['name']} at {t:g} K", [cas], [1.0], t))
    reforming = [CH4, H2O, H2, CO, CO2]
    cases.append(case("the owner's reforming feed at 1000 K", reforming, [10, 10, 0.1, 0, 0], 1000.0))
    cases.append(case("reforming gas mid-reactor at 1000 K", reforming,
                      [1.5, 1.2, 4.6, 1.1, 0.25], 1000.0))
    cases.append(case("reforming outlet at 1000 K", reforming,
                      [0.67277862, 0.45147590, 6.56407797, 1.88369644, 0.22130272], 1000.0))
    cases.append(case("methanol synthesis gas at 513.15 K", [CO2, H2, MEOH, H2O, CO],
                      [0.2, 0.65, 0.04, 0.03, 0.08], 513.15))
    cases.append(case("air at 300 K", [N2, O2], [0.79, 0.21], 300.0))
    for y in [0.1, 0.3, 0.5, 0.7, 0.9]:
        cases.append(case(f"hydrogen and carbon dioxide, y(H2) = {y:g}, at 500 K", [H2, CO2],
                          [y, 1 - y], 500.0))
    doc = {
        "source": ("Gas viscosities at low pressure: each compound's from its component file's "
                   "vapour-viscosity correlation (ChemSep's equation 102) by chemicals "
                   + chemicals.__version__ + "' EQ102, and a mixture's by its Wilke, with the "
                   "component files' molar masses; mole fractions normalised from the amounts "
                   "given. tools/gen_m9_viscosity.py"),
        "cases": cases,
    }
    OUT.write_bytes((json.dumps(doc, indent=1) + "\n").encode("utf-8"))
    print(f"{len(cases)} cases")


if __name__ == "__main__":
    main()
