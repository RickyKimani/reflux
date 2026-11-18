"""M6 step 0: which ThermoML archive entries hold phase-equilibrium data for
criterion 8's binaries. It asks the archive's public search API (NIST TRC,
https://trc.nist.gov/ThermoML/) for every entry with both compounds, keeps
the binaries with compositions at phase equilibrium, and lists each with
its DOI, year, title and data points.

    PYTHONIOENCODING=utf-8 python spike/thermoml/search.py

Not part of the simulator. Step 10 fetches the chosen files themselves.
"""
import json
import urllib.parse
import urllib.request

API = "https://trc.nist.gov/ThermoML-API/objects"
PAIRS = [("methanol", "water"), ("ethanol", "water"), ("1-butanol", "water")]
CONTEXT = json.dumps({"Citation": True, "data_summary": True, "Compound": True})


def search(a, b):
    q = f"type:TRCTml4 AND (/Compound/_/sCommonName/_:{a} AND /Compound/_/sCommonName/_:{b})"
    out, page = [], 0
    while True:
        url = API + "?" + urllib.parse.urlencode({"query": q, "pageNum": page, "pageSize": 100,
                                                  "requestContext": CONTEXT})
        req = urllib.request.Request(url, headers={"User-Agent": "curl/8.0"})  # the API refuses urllib's own
        with urllib.request.urlopen(req, timeout=120) as r:
            d = json.load(r)
        out += d["results"]
        if len(out) >= d["size"] or not d["results"]:
            return out
        page += 1


def walk(node, path=()):
    """Every property in a data summary, with its data points."""
    for k, v in node.items():
        if isinstance(v, dict):
            if "data_points" in v and not any(isinstance(x, dict) for x in v.values()):
                yield path + (k,), v["data_points"]
            else:
                yield from walk(v, path + (k,))


for a, b in PAIRS:
    hits = search(a, b)
    rows = []
    for h in hits:
        c = h["content"]
        binary = c.get("data_summary", {}).get("binary")
        if not binary or len(c.get("Compound", [])) != 2:
            continue
        props = [(p, n) for p, n in walk(binary) if any(w in " ".join(p) for w in
                                                       ["Composition", "equilibrium", "Equilibrium", "Boiling"])]
        if not props:
            continue
        cit = c["Citation"]
        title = (cit.get("sTitle") or [""])[0] if isinstance(cit.get("sTitle"), list) else cit.get("sTitle", "")
        rows.append((cit.get("TRCRefID", {}).get("yrYrPub", "?"), h["id"].split("/", 1)[1], title, props))
    rows.sort()
    print(f"\n== {a} + {b}: {len(hits)} entries, {len(rows)} binaries with phase-equilibrium compositions")
    for yr, doi, title, props in rows:
        summary = "; ".join(f"{p[-1]} ({n})" for p, n in props[:4])
        print(f"  {yr} {doi}: {title[:90]}")
        print(f"       {summary[:200]}")
