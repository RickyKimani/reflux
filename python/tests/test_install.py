"""Installing the package (D10, exit criterion 9). Built and installed into a
fresh environment from a copy away from the source tree, it uses the library
inside itself and finds the data in the source tree the library was built
from: kᵢⱼ and Mathias–Copeman constants included. Installed editable, as the README says to, it's the source
tree's own package."""

import os
import shutil
import subprocess
import sys

from support import ROOT

CHECK = """
import reflux
print(reflux.__file__)
print(reflux.data_root())
print(reflux.Database().fluid("PR", ["methane", "n-butane"]).kij[0, 1])
mc = reflux.Database().fluid("PR", ["n-butane"], alpha="mathias_copeman")
# Not "the equation's own: mc-PR.lib has no …", its fallback.
print(mc.setup["components"][0]["alpha"].startswith("Mathias–Copeman, mc-PR.lib line"))
"""


def run(args, **kwargs):
    r = subprocess.run(args, capture_output=True, text=True, **kwargs)
    assert r.returncode == 0, f"{args}: exit {r.returncode}\n{r.stdout}\n{r.stderr}"
    return r.stdout


def test_the_package_installs(tmp_path):
    src = tmp_path / "src"
    shutil.copytree(ROOT / "python", src, ignore=shutil.ignore_patterns(
        "tests", "__pycache__", ".pytest_cache", "*.egg-info", "build"))
    venv = tmp_path / "venv"
    # Nothing from this test run leaks in: the import path, the data root, the
    # library's path or uv's notion of the environment.
    env = {k: v for k, v in os.environ.items()
           if k not in ("PYTHONPATH", "REFLUX_DATA", "REFLUX_LIB", "VIRTUAL_ENV")}
    run(["uv", "venv", "--quiet", "--python", sys.executable, str(venv)], env=env)
    python = venv / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")

    run(["uv", "pip", "install", "--quiet", "--python", str(python), str(src)], env=env)
    where, root, kij, mc = run([str(python), "-c", CHECK], cwd=tmp_path, env=env).split()
    assert where.startswith(str(venv)) and where.endswith("__init__.py")
    assert root == str(ROOT / "data") and kij == "0.0244" and mc == "True"

    run(["uv", "pip", "install", "--quiet", "--python", str(python), "-e",
         str(ROOT / "python")], env=env)
    where, root, kij, mc = run([str(python), "-c", CHECK], cwd=tmp_path, env=env).split()
    assert where == str(ROOT / "python" / "reflux" / "__init__.py")
    assert root == str(ROOT / "data") and kij == "0.0244" and mc == "True"
