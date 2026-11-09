# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["sympy==1.13.3", "mpmath==1.3.0", "scipy==1.14.1"]
# ///
"""Generate M9's fixtures for the expression language.

Run from the repository root, with V on the PATH:

    uv run --python 3.13 tools/gen_m9_fixtures.py

It writes testdata/m9_expressions.json:

cases
    Expressions in reflux's syntax: hand-written ones for precedence, every
    function and every operator's spelling, and seeded random ones. SymPy
    reads each (as Python syntax, ^ as **), and mpmath at 50 digits gives
    its value at four points, with the largest magnitude of any of its
    subexpressions there, the scale a rounding error is measured against,
    and its derivative with respect to each input there (step 3).
    Then reflux writes each back as plain text, Unicode and LaTeX
    (tools/expr_render.v), and this script decodes each of those itself,
    with a small reader for exactly what reflux writes, and holds the result
    to SymPy's value at every point: the renderings mean what the text
    does. reflux's renderings are recorded, so that its tests see any change
    to them, and a change means running this again.

    SymPy's own LaTeX readers can't be the judge: in sympy 1.13.3 the ANTLR
    one reads \\log_{10}{\\left(x\\right)} \\cdot y as log(x·y), and the Lark
    one returns ambiguous trees for a function followed by an operator.

constants
    The built-in constants: the 2019 SI's defining constants as the SI
    Brochure (9th edition, table 1) gives them, R = N_A·k_B and F = N_A·e
    computed exactly in decimal, standard gravity (3rd CGPM, 1901) and the
    standard pressure (IUPAC, 1 bar). Each is checked against
    scipy.constants (CODATA 2018) to the digits SciPy has.

The output is the same, byte for byte, for the same versions.
"""

import json
import math
import pathlib
import random
import re
import subprocess
import tempfile
from decimal import Decimal, getcontext

import mpmath as mp
import scipy.constants as sc
import sympy

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m9_expressions.json"
mp.mp.dps = 50
getcontext().prec = 60

VARS = ["a", "b", "x", "y"]
MEMBERS = {"p[CH4]": "p_CH4", "p[H2O]": "p_H2O"}
SYMBOLS = {name: sympy.Symbol(name, positive=True) for name in VARS + list(MEMBERS.values())}  # every point is in [0.5, 2]

HAND = [
    # Precedence and associativity, each as Python reads it.
    "a + b * x", "(a + b) * x", "a - b - x", "a - (b - x)", "a - b + x",
    "a / b / x", "a / (b / x)", "a / b * x", "a * b / x", "a / (b * x)",
    "-a^2", "(-a)^2", "a^-2", "a^-b", "a^b^x", "(a^b)^x", "-a * b", "a * -b",
    "a + -b", "- -a", "2^-1", "-2^2", "+a", "a^+2", "a^(b + x)", "-(a + b)",
    "a * (b + x) / (y - 2)", "1 / (1 + a * x + b * y)^2",
    # Every function.
    "exp(-a / b)", "ln(a) + log10(b)", "sqrt(a * b)", "abs(a - b)", "exp(ln(a))",
    "sqrt(abs(a - b) + 1)", "ln(exp(a) + exp(b))", "log10(1000 * a)",
    # Numbers in every form, and pi.
    "a^0.5", "a^(1/3)", "x^2.5", "4.22e15 * x^-2.5", "1E+3 * a", ".5 * a",
    "3.5e-2 + a", "1e-3 / a", "pi * a", "2.0 * pi * x",
    # Every spelling of an operator, and a line break.
    "a·b", "a × b", "a − b", "a +\n b * x",
    # A family's members, and a rate's shape.
    "p[CH4]^2 / p[H2O]",
    "a * exp(-b / x) * (p[CH4] * p[H2O] - x^3 * y / exp(30.42 - 27106 / (1000 * a)))"
    " / (1 + 0.1 * p[CH4] + y * p[H2O] / x)^2",
]


def to_python(text):
    """reflux's syntax as Python's: ^ is **, ln is log, members are names."""
    s = text.replace("^", "**").replace("·", "*").replace("×", "*").replace("−", "-").replace("\n", " ")
    for k, v in MEMBERS.items():
        s = s.replace(k, v)
    s = re.sub(r"\bln\(", "log(", s)
    s = re.sub(r"\blog10\(", "log10(", s)
    return s


LOCALS = dict(SYMBOLS, log=sympy.log, log10=lambda z: sympy.log(z, 10), exp=sympy.exp,
              sqrt=sympy.sqrt, abs=sympy.Abs, Abs=sympy.Abs, pi=sympy.pi)


def sympy_of(python_text):
    return sympy.parse_expr(python_text, local_dict=LOCALS, evaluate=False)


def at(e, point):
    """e's value at a point, by mpmath at 50 digits, the point's values being
    the doubles reflux reads."""
    v = e.evalf(50, subs={SYMBOLS[k]: mp.mpf(val) for k, val in point.items()})
    return v


def scale_of(e, point):
    """The largest magnitude of any subexpression at the point."""
    m = 0
    for sub in sympy.preorder_traversal(e):
        if sub.is_number or sub.free_symbols:
            v = at(sub, point)
            if v.is_real is False or not v.is_number:
                return None
            m = max(m, abs(float(v)))
    return m


def points(rng):
    out = []
    for _ in range(4):
        out.append({k: round(rng.uniform(0.5, 2.0), 6) for k in SYMBOLS})
    return out


class Gen:
    """Seeded random expressions, written with the fewest parentheses
    reflux's and Python's shared precedence allow, or with all of them."""

    def __init__(self, seed):
        self.rng = random.Random(seed)

    def tree(self, depth):
        r = self.rng
        if depth == 0 or r.random() < 0.25:
            if r.random() < 0.6:
                return ("leaf", r.choice(VARS + list(MEMBERS)))
            return ("leaf", r.choice(["2", "3", "0.5", "1.5", "2.5e-1", "4e2", "7", "0.37", "pi", "1E+1"]))
        k = r.random()
        if k < 0.45:
            return (r.choice(["+", "-", "*", "/"]), self.tree(depth - 1), self.tree(depth - 1))
        if k < 0.6:
            return ("^", ("leaf", r.choice(VARS)), ("leaf", r.choice(["2", "3", "0.5", "-1", "-2", "1.5", "b"])))
        if k < 0.72:
            return ("neg", self.tree(depth - 1))
        return ("call", r.choice(["exp", "ln", "log10", "sqrt", "abs"]), self.tree(depth - 1))

    PREC = {"+": 1, "-": 1, "*": 2, "/": 2, "neg": 3, "^": 4}

    def write(self, t, full):
        kind = t[0]
        if kind == "leaf":
            return t[1]
        if kind == "call":
            return f"{t[1]}({self.write(t[2], full)})"
        if kind == "neg":
            inner = self.write(t[1], full)
            if full or self.prec(t[1]) < 3 or t[1][0] == "neg":
                inner = f"({inner})"
            return "-" + inner
        a, b = self.write(t[1], full), self.write(t[2], full)
        p = self.PREC[kind]
        if full:
            return f"({a} {kind} {b})"
        if self.prec(t[1]) < p or (kind == "^" and self.prec(t[1]) <= p):
            a = f"({a})"
        if self.prec(t[2]) < p or (kind in "-/" and self.prec(t[2]) <= p):
            b = f"({b})"
        return f"{a}^{b}" if kind == "^" else f"{a} {kind} {b}"

    def prec(self, t):
        return self.PREC.get(t[0], 5) if t[0] not in ("leaf", "call") else 5


def good(e, pts):
    """Real, finite and moderate at every point."""
    for pt in pts:
        try:
            v = at(e, pt)
        except Exception:
            return False
        if not v.is_number or v.is_real is False or v.has(sympy.zoo, sympy.nan, sympy.oo):
            return False
        f = abs(float(v))
        if not math.isfinite(f) or f > 1e100:
            return False
        s = scale_of(e, pt)
        if s is None or not math.isfinite(s) or s > 1e100:
            return False
    return True


# Decoding what reflux writes, to check it means the same.

GREEK = {"\\pi": "pi"}


class LatexReader:
    """Reads exactly the LaTeX reflux writes, into Python syntax."""

    def __init__(self, s):
        self.toks = re.findall(r"\\left\(|\\right\)|\\left\||\\right\||\\[a-zA-Z]+|\d+\.?\d*|\.\d+|[a-zA-Z]|[{}()^_+\-]|\S", s)
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else None

    def take(self, want=None):
        t = self.peek()
        if want is not None and t != want:
            raise ValueError(f"expected {want}, got {t} at {self.i} in {self.toks}")
        self.i += 1
        return t

    def read(self):
        out = self.expr()
        if self.peek() is not None:
            raise ValueError(f"trailing {self.toks[self.i:]}")
        return out

    def expr(self):
        out = self.term()
        while self.peek() in ("+", "-"):
            op = self.take()
            out = f"({out}) {op} ({self.term()})"
        return out

    def term(self):
        out = self.factor()
        while self.peek() in ("\\cdot", "\\times"):
            self.take()
            out = f"({out}) * ({self.factor()})"
        return out

    def factor(self):
        if self.peek() == "-":
            self.take()
            return f"-({self.factor()})"
        return self.power()

    def power(self):
        base = self.atom()
        if self.peek() == "^":
            self.take()
            return f"({base})**({self.group()})"
        return base

    def group(self):
        self.take("{")
        out = self.expr()
        self.take("}")
        return out

    def atom(self):
        t = self.take()
        if t == "\\frac":
            return f"({self.group()}) / ({self.group()})"
        if t == "\\sqrt":
            return f"sqrt({self.group()})"
        if t == "\\left(":
            inner = self.expr()
            self.take("\\right)")
            return f"({inner})"
        if t == "\\left|":
            inner = self.expr()
            self.take("\\right|")
            return f"Abs({inner})"
        if t in ("\\exp", "\\ln"):
            self.take("\\left(")
            inner = self.expr()
            self.take("\\right)")
            return f"{'exp' if t == '\\exp' else 'log'}({inner})"
        if t == "\\log":
            self.take("_")
            self.take("{")
            self.take("10")
            self.take("}")
            self.take("\\left(")
            inner = self.expr()
            self.take("\\right)")
            return f"log10({inner})"
        if t in GREEK:
            return GREEK[t]
        if re.fullmatch(r"\d+\.?\d*|\.\d+", t):
            return t
        if re.fullmatch(r"[a-zA-Z]", t):
            if self.peek() == "_":
                self.take("_")
                self.take("{")
                self.take("\\mathrm")
                self.take("{")
                name = ""
                while self.peek() != "}":
                    tok = self.take()
                    if tok == "_":
                        self.take("{")
                        name += self.take()
                        self.take("}")
                    else:
                        name += tok
                self.take("}")
                self.take("}")
                return f"{t}_{name}"
            return t
        raise ValueError(f"unexpected {t}")


SUPER = dict(zip("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-"))
SUB = dict(zip("₀₁₂₃₄₅₆₇₈₉", "0123456789"))


class UnicodeReader:
    """Reads exactly the Unicode reflux writes, into Python syntax."""

    def __init__(self, s):
        s = "".join(SUB.get(c, c) for c in s)
        self.toks = re.findall(r"log10\(|exp\(|ln\(|sqrt|√|[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+|\d+\.?\d*|\.\d+|p_[A-Za-z0-9]+|[a-zA-Z]|π|[()|·/+\-−^]", s.replace("log₁₀(", "log10("))
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else None

    def take(self, want=None):
        t = self.peek()
        if want is not None and t != want:
            raise ValueError(f"expected {want}, got {t} in {self.toks}")
        self.i += 1
        return t

    def read(self):
        out = self.expr()
        if self.peek() is not None:
            raise ValueError(f"trailing {self.toks[self.i:]}")
        return out

    def expr(self):
        out = self.term()
        while self.peek() in ("+", "−", "-"):
            self.take()
            op = "+" if self.toks[self.i - 1] == "+" else "-"
            out = f"({out}) {op} ({self.term()})"
        return out

    def term(self):
        out = self.factor()
        while self.peek() in ("·", "/"):
            op = "*" if self.take() == "·" else "/"
            out = f"({out}) {op} ({self.factor()})"
        return out

    def factor(self):
        if self.peek() in ("−", "-"):
            self.take()
            return f"-({self.factor()})"
        return self.power()

    def power(self):
        base = self.atom()
        t = self.peek()
        if t is not None and t[0] in SUPER:
            self.take()
            return f"({base})**({''.join(SUPER[c] for c in t)})"
        if t == "^":
            self.take()
            return f"({base})**({self.factor()})"
        return base

    def atom(self):
        t = self.take()
        if t == "(":
            inner = self.expr()
            self.take(")")
            return f"({inner})"
        if t == "|":
            inner = self.expr()
            self.take("|")
            return f"Abs({inner})"
        if t in ("exp(", "ln(", "log10("):
            inner = self.expr()
            self.take(")")
            return f"{ {'exp(': 'exp', 'ln(': 'log', 'log10(': 'log10'}[t] }({inner})"
        if t == "√":
            return f"sqrt({self.atom()})"
        if t == "π":
            return "pi"
        if re.fullmatch(r"\d+\.?\d*|\.\d+", t) or re.fullmatch(r"p_[A-Za-z0-9]+|[a-zA-Z]", t):
            return t
        raise ValueError(f"unexpected {t}")


def render(texts):
    """reflux's renderings, from tools/expr_render.v."""
    with tempfile.TemporaryDirectory() as d:
        path = pathlib.Path(d) / "cases.json"
        path.write_text(json.dumps(texts), encoding="utf-8")
        r = subprocess.run(["v", "run", str(ROOT / "tools" / "expr_render.v"), str(path)], cwd=ROOT,
                           capture_output=True, text=True, encoding="utf-8", check=True)
        return json.loads(r.stdout)


def same(e, python_text, pts, what, text):
    """The rendering means the text: equal at every point to 1e-12 of the
    scale. A rendering writes a number differently (4.22e15 as
    4.22 \\times 10^{15}, a binary 4.22 times 10¹⁵), which moves it by a
    rounding; a wrong precedence moves it by much more."""
    other = sympy_of(python_text)
    for pt in pts:
        a, b = at(e, pt), at(other, pt)
        s = scale_of(e, pt)
        if abs(a - b) > mp.mpf("1e-12") * max(1, s):
            raise SystemExit(f"{what} of {text!r} means something else: {python_text!r}: {a} vs {b} at {pt}")


def constants():
    n_a = Decimal("6.02214076e23")
    k_b = Decimal("1.380649e-23")
    e = Decimal("1.602176634e-19")
    out = [
        ("R", n_a * k_b, "J/(mol*K)", "N_A·k_B, exact", sc.value("molar gas constant")),
        ("N_A", n_a, "1/mol", "the SI's defining constant", sc.value("Avogadro constant")),
        ("k_B", k_b, "J/K", "the SI's defining constant", sc.value("Boltzmann constant")),
        ("h", Decimal("6.62607015e-34"), "J*s", "the SI's defining constant", sc.value("Planck constant")),
        ("F", n_a * e, "C/mol", "N_A·e, exact", sc.value("Faraday constant")),
        ("c", Decimal("299792458"), "m/s", "the SI's defining constant", sc.value("speed of light in vacuum")),
        ("g_n", Decimal("9.80665"), "m/s^2", "3rd CGPM (1901), exact", sc.value("standard acceleration of gravity")),
        ("pi", Decimal(mp.nstr(mp.pi, 40)), "", "π", math.pi),
        ("p_std", Decimal("100000"), "Pa", "IUPAC standard pressure, 1 bar, exact", 1e5),
    ]
    rows = []
    for name, exact, unit, source, scipy_value in out:
        nearest = float(exact)
        # SciPy's CODATA 2018 values: the defining constants in full; R and F
        # rounded to the digits CODATA prints, so agree to 1e-9.
        if abs(nearest - scipy_value) > 1e-9 * abs(nearest):
            raise SystemExit(f"{name}: {nearest} vs scipy's {scipy_value}")
        rows.append({"name": name, "exact": str(exact), "value": nearest, "unit": unit, "source": source,
                     "scipy": scipy_value})
    return rows


def main():
    rng = random.Random(20261110)
    gen = Gen(20261110)
    texts = list(HAND)
    trees = []
    while len(trees) < 240:
        t = gen.tree(4)
        text = gen.write(t, full=gen.rng.random() < 0.3)
        if text in texts:
            continue
        e = sympy_of(to_python(text))
        if not e.free_symbols:
            continue
        trees.append(text)
        texts.append(text)
    cases = []
    kept = []
    for text in texts:
        e = sympy_of(to_python(text))
        pts = points(rng)
        if not good(e, pts):
            if text in HAND:
                raise SystemExit(f"hand case {text!r} isn't real and finite at {pts}")
            continue
        kept.append((text, e, pts))
    rendered = render([k[0] for k in kept])
    for (text, e, pts), r in zip(kept, rendered):
        if r["error"]:
            raise SystemExit(f"reflux refused {text!r}: {r['error']}")
        same(e, to_python(r["plain"]), pts, "plain", text)
        same(e, LatexReader(r["latex"]).read(), pts, "LaTeX", text)
        same(e, UnicodeReader(r["unicode"]).read(), pts, "Unicode", text)
        cases.append({
            "text": text,
            "points": [[pt[k] for k in SYMBOLS] for pt in pts],
            "values": [float(at(e, pt)) for pt in pts],
            "scales": [scale_of(e, pt) for pt in pts],
            "gradients": [[float(at(sympy.diff(e, sym), pt)) for sym in SYMBOLS.values()] for pt in pts],
            "plain": r["plain"], "unicode": r["unicode"], "latex": r["latex"],
        })
    data = {
        "source": {
            "by": "tools/gen_m9_fixtures.py, with sympy 1.13.3, mpmath 1.3.0 at 50 digits, scipy 1.14.1",
            "inputs": "points give " + ", ".join(SYMBOLS) + ", in that order: the context's variables a, b, x, "
                      "y, then p[CH4] and p[H2O]",
            "values": "SymPy's reading of each text as Python syntax (^ as **), at each point, by mpmath at 50 "
                      "digits, rounded to the nearest double; scales are each point's largest magnitude of any "
                      "subexpression",
            "gradients": "SymPy's derivative of each case with respect to each input, in the inputs' order, at "
                         "each point, by mpmath at 50 digits, rounded to the nearest double",
            "renderings": "reflux's, from tools/expr_render.v, each decoded here and found to mean the text, to "
                          "50 digits at every point",
        },
        "cases": cases,
        "constants": constants(),
    }
    OUT.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}: {len(cases)} cases ({len(HAND)} by hand), {len(data['constants'])} constants")


main()
