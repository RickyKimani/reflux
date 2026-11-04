"""What the package's tests share: the repository's paths, the header's
declarations, and running code in a Python process of its own, so that a
crash is reported as one rather than taking the whole run with it (D11)."""

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
