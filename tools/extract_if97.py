# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pymupdf==1.24.14"]
# ///
"""Transcribe IAPWS-IF97's coefficients from its releases into V
(docs/m5-plan.md, D10): water/if97/coefficients.v, coefficients_sr3.v and
coefficients_sr5.v.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/extract_if97.py R7-97.pdf SR3-03.pdf SR5-05.pdf

The three PDFs are IAPWS's, from iapws.org (they aren't in the repository):
  R7-97(2012)    IF97-Rev.pdf
  SR3-03(2014)   Supp-Tv(ph,ps)3-2014.pdf
  SR5-05(2016)   Supp-VPT3-2016.pdf

Each table's rows are read from the PDF's text, span by span, with each
span's size and baseline, so that a superscript exponent is told from the
line it sits on, and SR3-03's and SR5-05's digits, set in the Symbol font
at private-use code points (U+F031 for 1), are mapped back to ASCII. Every
table must have exactly the number of terms its release gives, or nothing
is written. The coefficients are written as printed, mantissa and power of
ten; a test holds every one to the `iapws` package's
(testdata/m5_if97.json), and the releases' verification tables hold the
equations.
"""
import pathlib
import re
import sys

import pymupdf

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "water" / "if97"

# Symbol-font private-use characters, U+F000 + their code in the font.
SYM = {0xB4: "×", 0x2D: "-"}
HEAD = re.compile(r"^Table (A?\d+(?:\.\d+)?)\.")
INT = re.compile(r"^-?\d+(?:\.\d{1,2})?$")  # R7-97's Table 25 has I = -1.5, -1.25, ...
MERGED = re.compile(r"^(\d+)\s+(-?\s?\d\.\d{3} \d.*)$")
VALUE = re.compile(r"^-?\s?\d\.\d{3} \d")  # a coefficient: 0.xxx xxx ...
GLUED_POS = re.compile(r"^-\d{3,4}$")
INTRUN = re.compile(r"^[-\d\s]+$")


def norm(t):
    return t.replace("–", "-").replace("−", "-").strip()


def clean(t):
    """A span's text with the Symbol font's private-use characters mapped
    back, and en dashes and minus signs as hyphens."""
    out = []
    for ch in t:
        c = ord(ch)
        if 0xF000 <= c <= 0xF0FF:
            c -= 0xF000
            out.append(SYM.get(c, chr(c)))
        else:
            out.append(ch)
    return "".join(out).replace("–", "-").replace("−", "-")


def rows(page):
    """The page's lines: full-size spans grouped by baseline, each smaller
    span attached to the nearest line as a superscript (^) if raised."""
    spans = []
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                t = clean(s["text"])
                if t.strip():
                    spans.append({"x": s["bbox"][0], "size": s["size"], "text": t, "base": s["origin"][1]})
    out = []
    for s in sorted((s for s in spans if s["size"] >= 9.4), key=lambda s: s["base"]):
        if out and abs(out[-1][0] - s["base"]) < 2.5:
            out[-1][1].append(s)
        else:
            out.append([s["base"], [s]])
    for s in (s for s in spans if s["size"] < 9.4):
        best = min(out, key=lambda r: abs(r[0] - s["base"]), default=None)
        if best is None or abs(best[0] - s["base"]) > 7:
            out.append([s["base"], [s]])
            continue
        s = dict(s)
        s["text"] = ("^" if s["base"] < best[0] - 0.5 else "_") + s["text"]
        best[1].append(s)
    out.sort(key=lambda r: r[0])
    return out


def unglue(t):
    """Integers printed with no space between them, as the supplements'
    tables print i, I and J: "1-10 -8", "-8-12", "2-12". A minus sign inside
    a run starts a new integer; "-1210", with none, splits where |I| <= 14,
    |J| <= 64 and J has no leading zero, if only one split does."""
    if not INTRUN.match(t) or not any(c.isdigit() for c in t):
        return [t]
    out = []
    for piece in t.split():
        parts = re.findall(r"-?\d+", piece)
        if len(parts) == 1 and GLUED_POS.match(piece):
            ok = [(piece[:k], piece[k:]) for k in range(2, len(piece))
                  if abs(int(piece[:k])) <= 14 and abs(int(piece[k:])) <= 64
                  and not (piece[k] == "0" and len(piece) - k > 1)]
            if len(ok) == 1:
                parts = list(ok[0])
        out.extend(parts)
    return out


def tokens(row):
    out = []
    for s in sorted(row, key=lambda s: s["x"]):
        t = s["text"].strip()
        m = MERGED.match(t)
        if m and not t.startswith("^"):
            out += [(s["x"], m.group(1)), (s["x"] + 1, m.group(2))]
        else:
            out += [(s["x"] + j * 0.1, u) for j, u in enumerate(unglue(t))]
    # a minus sign set apart from its integer: R7-97's Table 22 prints I as
    # "–" and "7", 6 pt apart
    merged = []
    for x, t in out:
        if merged and merged[-1][1] == "-" and re.fullmatch(r"\d+", t) and x - merged[-1][0] < 10:
            merged[-1] = (merged[-1][0], "-" + t)
        else:
            merged.append((x, t))
    return merged


def make(ints, sign, mant, exp):
    m = mant.split("×")[0].replace(" ", "")
    if m.startswith("-"):
        sign, m = "-", m[1:]
    v = f"{sign}{m}" + (f"e{exp.replace(' ', '')}" if exp not in (None, "0") else "")
    return {"ints": ints, "n": v}


def entries(toks):
    """A row's entries: integers (i, I, J), then a coefficient with its sign
    and power of ten. None if the row doesn't read as entries."""
    res = []
    ints, mant, exp, sign = [], None, None, ""
    k = 0
    while k < len(toks):
        _, t = toks[k]
        if mant is None:
            if t[0] in "^_":
                pass  # a footnote mark before the coefficient
            elif t == "-":
                # a sign only right before its coefficient: R7-97's Table 30
                # prints n₁'s missing I and J as dashes, "1 – – 0.106…"
                if k + 1 < len(toks) and VALUE.match(toks[k + 1][1]) and toks[k + 1][0] - toks[k][0] < 12:
                    sign = "-"
            elif VALUE.match(t):
                mant = t
            elif INT.match(t):
                ints.append(float(t) if "." in t else int(t))
            else:
                return None
        elif t.startswith("^"):
            exp = norm(t[1:])
        elif t in ("×", "× 10") or (t == "10" and toks[k - 1][1] == "×"):
            pass
        elif INT.match(t) or t == "-" or VALUE.match(t):
            res.append(make(ints, sign, mant, exp))
            ints, mant, exp, sign = [], None, None, ""
            continue  # this token starts the next entry
        else:
            return None
        k += 1
    if mant is not None:
        res.append(make(ints, sign, mant, exp))
    elif ints:
        return None
    return res


def tables(path, tag):
    doc = pymupdf.open(path)
    found, cur, label = {}, None, None
    for k, page in enumerate(doc):
        for _, row in rows(page):
            toks = tokens(row)
            if not toks:
                continue
            m = HEAD.match(toks[0][1]) or HEAD.match(" ".join(t for _, t in toks[:2]))
            if m:
                cur, label = f"{tag}:{m.group(1)}", None
                found.setdefault(cur, [])
                continue
            # a boundary equation's name heads its first row: T3ab(p)
            lab = []
            while toks and toks[0][0] < 170 and not INT.match(toks[0][1]):
                lab.append(toks.pop(0)[1])
            if lab and cur:
                name = "".join(u.lstrip("_^") for u in lab if u not in ("()", "( )"))
                label = name[:-1] if name.endswith("p") else name
            if not (cur and toks and INT.match(toks[0][1]) and toks[0][0] < 200
                    and any(VALUE.match(u) for _, u in toks)):
                continue
            if any(re.search(r"[A-Za-z]", t) and t[0] not in "^_" for _, t in toks):
                continue  # a line of text with units in it, not a table's row
            es = entries(toks)
            if cur not in WANTED:
                continue  # a verification table: the tests type those in
            if es is None:
                raise SystemExit(f"{cur}, page {k + 1}: a row that doesn't read: {[t for _, t in toks]}")
            for e in es:
                e["label"] = label
                found[cur].append(e)
    for es in found.values():
        es.sort(key=lambda e: (e["label"] or "", e["ints"][0] if e["ints"] else 0))
    return found


# (key, name, columns, terms): the columns each row has after i.
R7 = [
    ("R7-97:1", "b23", "", 5, "Table 1, B23's n₁ to n₅ (Eqs. 5 and 6)"),
    ("R7-97:2", "r1", "IJ", 34, "Table 2, region 1's γ (Eq. 7)"),
    ("R7-97:6", "t1ph", "IJ", 20, "Table 6, region 1's backward T(p, h) (Eq. 11)"),
    ("R7-97:8", "t1ps", "IJ", 20, "Table 8, region 1's backward T(p, s) (Eq. 13)"),
    ("R7-97:10", "r2o", "J", 9, "Table 10, region 2's ideal-gas part γ° (Eq. 16)"),
    ("R7-97:11", "r2", "IJ", 43, "Table 11, region 2's residual part γʳ (Eq. 17)"),
    ("R7-97:16", "r2m", "IJ", 13, "Table 16, the metastable vapour's residual part γʳ (Eq. 19)"),
    ("R7-97:19", "b2bc", "", 5, "Table 19, B2bc's n₁ to n₅ (Eqs. 20 and 21)"),
    ("R7-97:20", "t2aph", "IJ", 34, "Table 20, subregion 2a's backward T(p, h) (Eq. 22)"),
    ("R7-97:21", "t2bph", "IJ", 38, "Table 21, subregion 2b's backward T(p, h) (Eq. 23)"),
    ("R7-97:22", "t2cph", "IJ", 23, "Table 22, subregion 2c's backward T(p, h) (Eq. 24)"),
    ("R7-97:25", "t2aps", "IJ", 46, "Table 25, subregion 2a's backward T(p, s) (Eq. 25); its I are quarters"),
    ("R7-97:26", "t2bps", "IJ", 44, "Table 26, subregion 2b's backward T(p, s) (Eq. 26)"),
    ("R7-97:27", "t2cps", "IJ", 30, "Table 27, subregion 2c's backward T(p, s) (Eq. 27)"),
    ("R7-97:30", "r3", "IJ", 40, "Table 30, region 3's φ (Eq. 28); n₁ multiplies ln δ"),
    ("R7-97:34", "r4", "", 10, "Table 34, region 4's n₁ to n₁₀ (Eqs. 29 to 31)"),
    ("R7-97:37", "r5o", "J", 6, "Table 37, region 5's ideal-gas part γ° (Eq. 33)"),
    ("R7-97:38", "r5", "IJ", 6, "Table 38, region 5's residual part γʳ (Eq. 34)"),
]
SR3 = [
    ("SR3-03:2", "h3ab", "", 4, "Table 2, the boundary h3ab(p) (Eq. 1)"),
    ("SR3-03:3", "t3aph", "IJ", 31, "Table 3, T3a(p, h) (Eq. 2)"),
    ("SR3-03:4", "t3bph", "IJ", 33, "Table 4, T3b(p, h) (Eq. 3)"),
    ("SR3-03:6", "v3aph", "IJ", 32, "Table 6, v3a(p, h) (Eq. 4)"),
    ("SR3-03:7", "v3bph", "IJ", 30, "Table 7, v3b(p, h) (Eq. 5)"),
    ("SR3-03:10", "t3aps", "IJ", 33, "Table 10, T3a(p, s) (Eq. 6)"),
    ("SR3-03:11", "t3bps", "IJ", 28, "Table 11, T3b(p, s) (Eq. 7)"),
    ("SR3-03:13", "v3aps", "IJ", 28, "Table 13, v3a(p, s) (Eq. 8)"),
    ("SR3-03:14", "v3bps", "IJ", 31, "Table 14, v3b(p, s) (Eq. 9)"),
    ("SR3-03:17", "p3sath", "IJ", 14, "Table 17, the boundary p3sat(h) (Eq. 10)"),
    ("SR3-03:19", "p3sats", "IJ", 10, "Table 19, the boundary p3sat(s) (Eq. 11)"),
]
# SR5-05's Tables 4 and 12 give each subregion's N.
SUB_N = dict(zip("abcdefghijklmnopqrst", [30, 32, 35, 38, 29, 42, 38, 29, 42, 29, 34, 43, 40, 39, 24, 27, 24, 27, 29, 33]))
SUB_N.update(dict(zip("uvwxyz", [38, 39, 35, 36, 20, 23])))
SR5 = [(f"SR5-05:A1.{k + 1}", f"v3{c}", "IJ", SUB_N[c], f"Table A1.{k + 1}, v3{c}(p, T) (Eq. {'5' if c == 'n' else '4'})")
       for k, c in enumerate("abcdefghijklmnopqrst")]
SR5 += [(f"SR5-05:A2.{k + 1}", f"v3{c}", "IJ", SUB_N[c], f"Table A2.{k + 1}, v3{c}(p, T) (Eq. 4)") for k, c in enumerate("uvwxyz")]
# Tables 1 and 9: each boundary equation's n and I, by its name.
BOUNDS = {"T3ab": 5, "T3cd": 4, "T3gh": 5, "T3ij": 5, "T3jk": 5, "T3mn": 4, "T3op": 5, "T3qu": 4, "T3rx": 4,
          "T3uv": 4, "T3wx": 5}


WANTED = {k for k, *_ in R7 + SR3 + SR5} | {"SR5-05:1", "SR5-05:9"}


def v_list(xs, ty):
    if ty == "int":
        items = [str(int(x)) for x in xs]
    else:
        items = []
        for x in xs:
            s = str(x)
            if ty == "f64" and "." not in s and "e" not in s:
                s += ".0"
            items.append(s)
    lines, cur = [], ""
    for it in items:
        piece = it + ", "
        if len(cur) + len(piece) > 96 and cur:
            lines.append(cur.rstrip())
            cur = ""
        cur += piece
    lines.append(cur.rstrip().rstrip(","))
    return "[" + "\n\t".join(lines) + "]"


def emit_table(out, es, name, cols, terms, what):
    if len(es) != terms:
        raise SystemExit(f"{what}: {len(es)} terms, not {terms}")
    idx = [e["ints"][0] for e in es if e["ints"]]
    if idx and idx != list(range(1, len(idx) + 1)) and not (name == "r4" and len(idx) == 10):
        if not (name == "r3" and idx == list(range(1, 41))):
            raise SystemExit(f"{what}: rows {idx}")
    out.append(f"// {what}.")
    if name == "r3":
        n1 = [e for e in es if len(e["ints"]) == 1]
        es = [e for e in es if len(e["ints"]) == 3]
        out.append(f"const r3_n1 = {n1[0]['n']}")
    if cols == "IJ":
        bad = [e for e in es if len(e["ints"]) != 3]
        if bad:
            raise SystemExit(f"{what}: rows without i, I and J: {bad}")
        ity = "f64" if name == "t2aps" else "int"
        out.append(f"const {name}_i = {v_list([e['ints'][1] for e in es], ity)}")
        out.append(f"const {name}_j = {v_list([e['ints'][2] for e in es], 'int')}")
    elif cols == "J":
        out.append(f"const {name}_j = {v_list([e['ints'][1] for e in es], 'int')}")
    out.append(f"const {name}_n = {v_list([e['n'] for e in es], 'num')}")
    out.append("")


HEADER = {
    "coefficients.v": """// IAPWS-IF97's coefficients, transcribed from IAPWS R7-97(2012), "Revised
// Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic
// Properties of Water and Steam" (Wagner et al., J. Eng. Gas Turbines Power
// 122, 150, 2000), by tools/extract_if97.py from the release's PDF. IAPWS
// allows the release's publication in whole or in part with attribution. A
// test holds every number to the `iapws` package's (testdata/m5_if97.json),
// and the release's verification tables hold the equations.
""",
    "coefficients_sr3.v": """// The coefficients of IAPWS SR3-03(2014), "Revised Supplementary Release on
// Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s)
// for Region 3 of the IAPWS Industrial Formulation 1997", transcribed by
// tools/extract_if97.py from the release's PDF. A test holds every number
// to the `iapws` package's (testdata/m5_if97.json), and the release's
// verification tables hold the equations.
""",
    "coefficients_sr5.v": """// The coefficients of IAPWS SR5-05(2016), "Revised Supplementary Release on
// Backward Equations for Specific Volume as a Function of Pressure and
// Temperature v(p,T) for Region 3 of the IAPWS Industrial Formulation 1997",
// transcribed by tools/extract_if97.py from the release's PDF: the
// subregion boundaries (Tables 1 and 9) and the equations of the subregions
// 3a to 3t (Tables A1.1 to A1.20) and 3u to 3z (Tables A2.1 to A2.6). A test
// holds every number to the `iapws` package's (testdata/m5_if97.json), and
// the release's verification tables hold the equations.
""",
}


def main():
    pdf97, pdf3, pdf5 = sys.argv[1:4]
    found = {}
    for path, tag in ((pdf97, "R7-97"), (pdf3, "SR3-03"), (pdf5, "SR5-05")):
        found.update(tables(path, tag))
    files = {"coefficients.v": R7, "coefficients_sr3.v": SR3, "coefficients_sr5.v": SR5}
    for fname, specs in files.items():
        out = ["module if97", "", HEADER[fname]]
        if fname == "coefficients.v":
            out[-1] += """
// Table 10's footnote: n₁° and n₂° of γ° when it's part of the metastable
// vapour's equation, Eq. (18), in place of Table 10's.
const r2m_n01 = -0.96937268393049e1
const r2m_n02 = 0.10087275970006e2
"""
        if fname == "coefficients_sr5.v":
            out.append("// Tables 1 and 9: the subregion boundaries' nᵢ and Iᵢ, Eqs. (1) and (2).")
            es = found["SR5-05:1"] + found["SR5-05:9"]
            for name, terms in BOUNDS.items():
                rows_ = [e for e in es if e["label"] == name]
                if len(rows_) != terms or [e["ints"][0] for e in rows_] != list(range(1, terms + 1)):
                    raise SystemExit(f"{name}: {[(e['ints'], e['n']) for e in rows_]}")
                low = name.lower()
                out.append(f"const {low}_i = {v_list([e['ints'][1] for e in rows_], 'int')}")
                out.append(f"const {low}_n = {v_list([e['n'] for e in rows_], 'num')}")
            out.append("")
        for key, name, cols, terms, what in specs:
            emit_table(out, found[key], name, cols, terms, what)
        text = "\n".join(out).rstrip() + "\n"
        (OUT / fname).write_bytes(text.encode("utf-8"))
        print(f"wrote {OUT / fname}: {sum(t for *_, t, _ in specs)} terms")


if __name__ == "__main__":
    main()
