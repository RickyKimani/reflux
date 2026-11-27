"""M5b step 0: which ThermoML archive entries hold mixtures' viscosity,
thermal conductivity and surface tension, for candidate binaries. It asks
the archive's public search API (NIST TRC, https://trc.nist.gov/ThermoML/)
for every entry with both compounds, as M6's spike/thermoml/search.py did,
and lists per pair and property the binary data sets: year, DOI, points.

    PYTHONIOENCODING=utf-8 python -I spike/transport/thermoml.py

Not part of the simulator. The data step fetches the chosen files.
"""
import json
import urllib.parse
import urllib.request

API = "https://trc.nist.gov/ThermoML-API/objects"
CONTEXT = json.dumps({"Citation": True, "data_summary": True, "Compound": True})
PAIRS = [("methanol", "water"), ("ethanol", "water"), ("water", "1-propanol"), ("acetone", "water"),
         ("benzene", "toluene"), ("hexane", "heptane"), ("benzene", "cyclohexane"), ("toluene", "heptane"),
         ("methanol", "ethanol"), ("acetone", "methanol"),
         ("methane", "ethane"), ("methane", "propane"), ("methane", "butane"), ("nitrogen", "methane"),
         ("carbon dioxide", "methane"), ("hydrogen", "methane"), ("hydrogen", "nitrogen"),
         ("hydrogen", "carbon dioxide"), ("hydrogen", "carbon monoxide"), ("nitrogen", "oxygen"),
         ("helium", "argon"), ("argon", "methane")]
KINDS = {"viscosity": "Viscosity", "conductivity": "Thermal conductivity", "surface tension": "Surface tension"}


def search(a, b):
    q = f"type:TRCTml4 AND (/Compound/_/sCommonName/_:\"{a}\" AND /Compound/_/sCommonName/_:\"{b}\")"
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


def main():
    for a, b in PAIRS:
        hits = search(a, b)
        found = {k: [] for k in KINDS}
        for h in hits:
            c = h["content"]
            binary = c.get("data_summary", {}).get("binary")
            if not binary or len(c.get("Compound", [])) != 2:
                continue
            for p, n in walk(binary):
                s = " / ".join(p).lower()
                for k, w in KINDS.items():
                    if w.lower() in s:
                        year = c["Citation"].get("TRCRefID", {}).get("yrYrPub", "?")
                        found[k].append((year, h["id"].split("/", 1)[1], n))
        print(f"== {a} + {b}: {len(hits)} entries")
        for k, rows in found.items():
            if rows:
                pts = sum(r[2] for r in rows if isinstance(r[2], int))
                print(f"   {k}: {len(rows)} sets, {pts} points: {sorted(rows)}")


if __name__ == "__main__":
    main()
