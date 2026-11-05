"""What the package's tests share: the repository's paths, the header's
declarations, and running code in a Python process of its own, so that a
crash is reported as one rather than taking the whole run with it (D11)."""

import json
import os
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
HEADER = (ROOT / "capi" / "reflux.h").read_text(encoding="utf-8")


def header_functions():
    """The functions capi/reflux.h declares."""
    return re.findall(r"^\w[\w \*]*?\b(reflux_\w+)\(", HEADER, re.M)


def header_defines():
    """The header's #define REFLUX_* values, without the prefix."""
    return {name: int(value) for name, value in
            re.findall(r"^#define REFLUX_(\w+) \(?(-?\d+)\)?", HEADER, re.M)}


def v_mod_version():
    """v.mod's version."""
    return re.search(r"version:\s*'([^']+)'", (ROOT / "v.mod").read_text())[1]


def isolated(code, timeout=300, env=None):
    """Run code in a Python process of its own, with this package importable
    and env's variables set, and return what it printed. A crash fails the
    test with its exit status."""
    env = {**os.environ, **(env or {}), "PYTHONPATH": str(ROOT / "python")}
    r = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True,
                       timeout=timeout)
    assert r.returncode == 0, f"exit {r.returncode}\n{r.stdout}\n{r.stderr}"
    return r.stdout


# M3's flash references (testdata/m3_thermo.json), from thermo 0.6.1, and the
# CAS numbers of their eleven compounds.
M3 = json.loads((ROOT / "testdata" / "m3_thermo.json").read_text(encoding="utf-8"))
M3_CAS = {"methane": "74-82-8", "ethane": "74-84-0", "propane": "74-98-6",
          "n-butane": "106-97-8", "n-pentane": "109-66-0", "n-hexane": "110-54-3",
          "n-heptane": "142-82-5", "n-decane": "124-18-5", "carbon dioxide": "124-38-9",
          "nitrogen": "7727-37-9", "isobutane": "75-28-5"}


def m3_compounds(directory):
    """Write M3's eleven compounds into directory as the caller's own: each
    the database's file with M3's Tc, Pc, ω and molar mass. Their ideal-gas Cp
    stays ChemSep's, so only what the model alone decides compares with the
    fixture: splits, volumes, boundaries."""
    for c in M3["components"]:
        cas = M3_CAS[c["name"]]
        text = (ROOT / "data" / "components" / f"{cas}.toml").read_text(encoding="utf-8")
        for key in ("tc", "pc", "omega", "mw"):
            text, n = re.subn(rf"^{key} = .*$", f"{key} = {c[key]!r}", text, flags=re.M)
            assert n == 1, (cas, key)
        (directory / f"{cas}.toml").write_text(text, encoding="utf-8")
    return directory


def m3_fluid(db, case):
    """A fixture case's fluid from M3's database: its equation, compounds and
    kij, every pair set by hand."""
    names = [M3["components"][i]["name"] for i in case["components"]]
    nc = len(names)
    k = case["kij"]
    overrides = {(names[i], names[j]): k[i * nc + j] for i in range(nc) for j in range(i + 1, nc)}
    return db.fluid(case.get("eos", "pr"), [M3_CAS[n] for n in names], kij=False,
                    overrides=overrides)
