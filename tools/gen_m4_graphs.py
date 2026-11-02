# /// script
# requires-python = ">=3.11"
# dependencies = ["networkx==3.4.2", "scipy==1.14.1"]
# ///
"""Generate testdata/m4_graphs.json, the graphs of M4's exit criterion 1.

Run from the repository root:

    uv run tools/gen_m4_graphs.py

Three groups of directed graphs, each with its strongly connected components,
its elementary cycles, the exact minimum number of tears, and the fewest
second breaks a minimum tear set can have. flowsheet/graph_test.v checks
reflux against those.

- literature: test problems 1 to 9 of T. Gundersen and T. Hertzberg,
  "Partitioning and tearing of networks, applied to process flowsheeting",
  Modeling, Identification and Control 4(3), 139-165 (1983),
  doi:10.4173/mic.1983.3.2. They're "the classical graphs from the chemical
  engineering literature", which are drawn in its Appendix A. The edge lists
  here were read off those figures by hand, so each is held to the paper's
  own Tables 1 and 3: its vertices, edges, maximal cyclical nets, elementary
  cycles and minimum tears. The script stops if any differs. Problem 10, the
  heavy water plant, is too dense a figure to read reliably.
- worked: small graphs of our own, each with its answer worked by hand here.
- random: seeded graphs, with self-loops and parallel streams.

The answers are computed independently of reflux: components and cycles by
networkx, and the tears by an exact integer program (scipy's HiGHS), over the
cycles. Parallel streams between two units are separate cycles, so a cycle
here is a list of edges, not of nodes.
"""

import datetime
import itertools
import json
import pathlib
import random

import networkx as nx
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m4_graphs.json"
RANDOM_GRAPHS = 40
RANDOM_SEED = 20261102


def chain(*nodes):
    """Edges along a path of nodes."""
    return [(a, b) for a, b in zip(nodes, nodes[1:])]


def net(a, b, c, d, e):
    """Barkley and Motard's five-unit net, as Problem 3 repeats it three times."""
    return chain(a, b, c, d, e) + [(d, a), (d, b), (a, c), (c, e), (e, b), (e, a)]


LITERATURE = [
    {
        "name": "problem 1: complete graph",
        "published": {"vertices": 6, "edges": 36, "nets": 1, "cycles": 415, "min_tears": 21},
        "edges": [(a, b) for a in range(1, 7) for b in range(1, 7)],
    },
    {
        "name": "problem 2: Pho and Lapidus",
        "published": {"vertices": 12, "edges": 21, "nets": 1, "cycles": 22, "min_tears": 2},
        "edges": chain(1, 2, 3, 4, 5) + chain(6, 7, 8, 9) + [
            (5, 1), (1, 4), (1, 6), (7, 2), (3, 8), (9, 4), (5, 10), (12, 10), (10, 6),
            (11, 12), (7, 11), (9, 11), (6, 4), (12, 4),
        ],
    },
    {
        "name": "problem 3: Barkley and Motard",
        "published": {"vertices": 15, "edges": 35, "nets": 3, "cycles": 27, "min_tears": 6},
        "edges": net(1, 2, 3, 4, 5) + net(6, 7, 8, 9, 10) + net(11, 12, 13, 14, 15) + [
            (2, 6), (5, 7), (8, 13), (10, 11), (9, 12),
        ],
    },
    {
        "name": "problem 4: Sargent and Westerberg",
        "published": {"vertices": 19, "edges": 31, "nets": 1, "cycles": 20, "min_tears": 6},
        "edges": [
            (4, 1), (1, 5), (6, 2), (2, 5), (6, 3), (3, 7), (4, 5), (5, 6), (6, 7), (4, 8),
            (8, 14), (14, 15), (15, 16), (16, 10), (10, 9), (9, 4), (6, 10), (11, 6), (17, 11),
            (6, 12), (12, 7), (12, 17), (17, 16), (17, 18), (18, 17), (18, 7), (19, 18), (7, 19),
            (19, 13), (13, 7), (7, 13),
        ],
    },
    {
        "name": "problem 5: Christensen and Rudd, the first",
        "published": {"vertices": 25, "edges": 32, "nets": 1, "cycles": 10, "min_tears": 3},
        "edges": chain(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) + chain(11, 21, 22, 23, 24, 25) + [
            (3, 11), (6, 12), (12, 11), (25, 14), (14, 13), (13, 23), (13, 15), (16, 15), (17, 16),
            (16, 8), (10, 17), (10, 20), (15, 18), (18, 19), (19, 20), (20, 1), (25, 20), (21, 20),
        ],
    },
    {
        "name": "problem 6: Jain and Eakman, HF alkylation",
        "published": {"vertices": 29, "edges": 37, "nets": 1, "cycles": 11, "min_tears": 5},
        "edges": [
            (23, 3), (22, 1), (8, 1), (7, 1), (1, 2), (2, 7), (2, 3), (3, 10), (3, 4), (3, 5),
            (5, 3), (3, 6), (4, 25), (6, 24), (9, 10), (10, 11), (11, 12), (12, 26), (10, 15),
            (7, 9), (9, 14), (13, 9), (15, 13), (15, 17), (17, 15), (15, 16), (16, 26), (15, 18),
            (18, 19), (19, 20), (20, 18), (20, 21), (21, 29), (19, 26), (26, 8), (14, 28), (14, 27),
        ],
    },
    {
        "name": "problem 7: Christensen and Rudd, the second",
        "published": {"vertices": 30, "edges": 42, "nets": 1, "cycles": 31, "min_tears": 3},
        "edges": chain(1, 2, 3, 4, 5, 6) + chain(7, 8, 9, 10) + chain(14, 15, 16, 17, 18, 19)
        + chain(23, 24, 25, 26) + [
            (1, 5), (6, 8), (10, 11), (11, 1), (11, 12), (12, 13), (13, 2), (4, 14), (14, 19),
            (17, 19), (19, 7), (13, 15), (13, 23), (18, 23), (12, 19), (25, 19), (26, 22),
            (22, 21), (21, 12), (10, 20), (20, 29), (29, 30), (29, 28), (28, 27), (27, 29), (30, 22),
        ],
    },
    {
        "name": "problem 8: Shannon, sulfuric acid",
        "published": {"vertices": 41, "edges": 61, "nets": 1, "cycles": 103, "min_tears": 5},
        "edges": chain(1, 2, 3, 4, 5, 6, 25, 27, 7, 8, 9, 10, 11, 12, 13, 14, 29, 30, 31, 33, 34, 35, 36, 37, 40)
        + chain(26, 28, 21, 20, 19, 17, 16, 15, 14) + chain(40, 39, 38, 37) + chain(2, 24, 26) + [
            (2, 6), (8, 5), (5, 7), (4, 6), (8, 27), (27, 26), (7, 18), (14, 7), (24, 20), (24, 16),
            (34, 32), (32, 31), (38, 33), (40, 30), (29, 41), (41, 40), (21, 22), (22, 23), (23, 1),
            (17, 22), (18, 21), (18, 17), (1, 33), (36, 1),
        ],
    },
    {
        "name": "problem 9: Jain and Eakman, vegetable oil",
        "published": {"vertices": 50, "edges": 79, "nets": 1, "cycles": 22, "min_tears": 8},
        "edges": [
            (38, 30), (33, 30), (30, 28), (28, 39), (28, 26), (26, 28), (30, 29), (29, 30), (20, 50),
            (50, 26), (26, 50), (26, 29), (29, 26), (29, 31), (31, 29), (31, 32), (32, 31), (29, 27),
            (27, 29), (50, 27), (27, 32), (32, 35), (34, 38), (34, 35), (35, 34), (38, 36), (18, 20),
            (17, 18), (18, 19), (19, 49), (49, 47), (19, 47), (17, 47), (13, 14), (14, 15), (15, 17),
            (16, 17), (7, 16), (12, 16), (7, 13), (11, 13), (10, 11), (9, 10), (8, 9), (7, 8), (1, 9),
            (2, 8), (8, 12), (3, 12), (6, 17), (17, 21), (21, 22), (22, 23), (23, 24), (24, 47),
            (24, 25), (21, 25), (12, 25), (5, 17), (5, 36), (36, 25), (48, 13), (48, 42), (4, 13),
            (4, 42), (35, 25), (35, 40), (40, 41), (41, 47), (41, 42), (35, 37), (37, 42), (42, 44),
            (44, 43), (43, 47), (43, 46), (45, 46), (46, 47), (36, 35),
        ],
    },
]

WORKED = [
    # No cycle at all: each unit is its own block, in order.
    {"name": "a chain", "vertices": 4, "edges": chain(1, 2, 3, 4),
     "answer": {"nets": 0, "cycles": 0, "min_tears": 0}},
    # One loop of three: any one of its streams will do.
    {"name": "one loop", "vertices": 4, "edges": chain(1, 2, 3, 4) + [(3, 2)],
     "answer": {"nets": 1, "cycles": 1, "min_tears": 1}},
    # A unit's outlet back to its own inlet: a net of one.
    {"name": "a self-loop", "vertices": 2, "edges": [(1, 2), (2, 2)],
     "answer": {"nets": 1, "cycles": 1, "min_tears": 1}},
    # 2-3-4-5-2 with 3-4-3 inside it. Tearing 3-4 breaks both, once each.
    {"name": "nested loops", "vertices": 5, "edges": chain(1, 2, 3, 4, 5) + [(5, 2), (4, 3)],
     "answer": {"nets": 1, "cycles": 2, "min_tears": 1}},
    # 1-2-3-1 and 2-3-4-2 share stream 2-3. It alone breaks both.
    {"name": "interlocking loops", "vertices": 4, "edges": chain(1, 2, 3, 4) + [(3, 1), (4, 2)],
     "answer": {"nets": 1, "cycles": 2, "min_tears": 1}},
    # 1-2-1 and 3-4-3, joined one way only: two nets, a tear each.
    {"name": "two nets in series", "vertices": 4, "edges": [(1, 2), (2, 1), (2, 3), (3, 4), (4, 3)],
     "answer": {"nets": 2, "cycles": 2, "min_tears": 2}},
    # Vapour and liquid both go from 2 back to 1: two cycles sharing 1-2.
    {"name": "parallel streams", "vertices": 2, "edges": [(1, 2), (2, 1), (2, 1)],
     "answer": {"nets": 1, "cycles": 2, "min_tears": 1}},
    # Where the smallest set can't help tearing a cycle twice. Cycles: a = 1-2-1,
    # b = 2-3-2, c = 1-2-3-1. Any single stream misses one of them, so 2 tears.
    # 2-1 with 3-2 cover all three with no second break; 1-2 with 2-3 tear c
    # twice.
    {"name": "a redundant pair to avoid", "vertices": 3, "edges": [(1, 2), (2, 1), (2, 3), (3, 2), (3, 1)],
     "answer": {"nets": 1, "cycles": 3, "min_tears": 2}},
    # The M4 reference flowsheet: mixer 1, cooler 2, drum 3, splitter 4, pump
    # 5, valve 6, drum 7. The liquid recycle is the one loop.
    {"name": "the reference flowsheet", "vertices": 7, "edges": chain(1, 2, 3, 4, 5, 1) + chain(4, 6, 7),
     "answer": {"nets": 1, "cycles": 1, "min_tears": 1}},
    # Pyomo's simple_recycle_model: mixer 1, unit 2, splitter 3. Its test says
    # any one of the three streams is a tear.
    {"name": "Pyomo's simple recycle", "vertices": 3, "edges": chain(1, 2, 3, 1),
     "answer": {"nets": 1, "cycles": 1, "min_tears": 1}},
]


def cycles_of(n, edges):
    """Every elementary cycle, as a sorted list of edge indices. Parallel edges
    make separate cycles."""
    g = nx.DiGraph()
    g.add_nodes_from(range(1, n + 1))
    between = {}
    for i, (a, b) in enumerate(edges):
        g.add_edge(a, b)
        between.setdefault((a, b), []).append(i)
    out = []
    for nodes in nx.simple_cycles(g):
        hops = list(zip(nodes, nodes[1:] + nodes[:1]))
        for choice in itertools.product(*(between[h] for h in hops)):
            out.append(sorted(choice))
    return sorted(out)


def nets_of(n, edges):
    """The maximal cyclical nets, as sorted lists of nodes: strongly connected
    components with a cycle inside."""
    g = nx.DiGraph()
    g.add_nodes_from(range(1, n + 1))
    g.add_edges_from(edges)
    loops = {a for a, b in edges if a == b}
    nets = [sorted(c) for c in nx.strongly_connected_components(g) if len(c) > 1 or set(c) & loops]
    return sorted(nets)


def solve(cycles, n_edges, weights, total=None):
    """The minimum of weights·x for x in {0, 1} per edge, hitting every cycle,
    with Σx = total if given."""
    a = lil_matrix((len(cycles) + (total is not None), n_edges))
    for row, cycle in enumerate(cycles):
        for e in cycle:
            a[row, e] = 1
    lower = [1.0] * len(cycles)
    upper = [np.inf] * len(cycles)
    if total is not None:
        a[len(cycles), :] = 1
        lower.append(total)
        upper.append(total)
    r = milp(np.array(weights, dtype=float), integrality=np.ones(n_edges), bounds=Bounds(0, 1),
             constraints=LinearConstraint(a.tocsr(), lower, upper))
    if not r.success:
        raise SystemExit(f"the integer program failed: {r.message}")
    return int(round(r.fun))


def answers(n, edges):
    cycles = cycles_of(n, edges)
    out = {"nets": nets_of(n, edges), "cycles": len(cycles), "cycle_list": cycles, "min_tears": 0,
           "second_breaks": 0}
    if cycles:
        k = solve(cycles, len(edges), [1.0] * len(edges))
        hits = [sum(e in c for c in cycles) for e in range(len(edges))]
        out["min_tears"] = k
        # Σ over cycles of (breaks − 1), for x a minimum tear set: Σ hits·x − cycles.
        out["second_breaks"] = solve(cycles, len(edges), hits, k) - len(cycles)
    return out


def row(name, source, n, edges):
    return {"name": name, "source": source, "vertices": n, "edges": [list(e) for e in edges], **answers(n, edges)}


def literature():
    out = []
    for g in LITERATURE:
        p = g["published"]
        r = row(g["name"], "Gundersen and Hertzberg 1983, Appendix A, Tables 1 and 3", p["vertices"], g["edges"])
        got = {"vertices": r["vertices"], "edges": len(r["edges"]), "nets": len(r["nets"]),
               "cycles": r["cycles"], "min_tears": r["min_tears"]}
        if got != p:
            raise SystemExit(f"{g['name']}: {got} but the paper has {p}: the edge list is misread")
        out.append(r)
    return out


def worked():
    out = []
    for g in WORKED:
        r = row(g["name"], "worked by hand", g["vertices"], g["edges"])
        a = g["answer"]
        got = {"nets": len(r["nets"]), "cycles": r["cycles"], "min_tears": r["min_tears"]}
        if got != a:
            raise SystemExit(f"{g['name']}: {got} but worked by hand it's {a}")
        out.append(r)
    return out


def seeded():
    rng = random.Random(RANDOM_SEED)
    out = []
    for k in range(RANDOM_GRAPHS):
        n = rng.randint(3, 12)
        m = rng.randint(n, 3 * n)
        edges = [(rng.randint(1, n), rng.randint(1, n)) for _ in range(m)]
        r = row(f"random {k}", f"random.Random({RANDOM_SEED}), graph {k}", n, edges)
        out.append(r)
    return out


def main():
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "by": "tools/gen_m4_graphs.py, with networkx 3.4.2 and scipy 1.14.1's HiGHS",
            "note": "nets are the strongly connected components with a cycle inside. "
                    "second_breaks is the fewest extra breaks a minimum tear set can have: "
                    "the sum over cycles of the times it is torn, less 1.",
        },
        "graphs": literature() + worked() + seeded(),
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    for g in data["graphs"]:
        print(f"{g['name']:48} {g['vertices']:3} units {len(g['edges']):3} streams {len(g['nets']):2} nets "
              f"{g['cycles']:4} cycles {g['min_tears']:3} tears {g['second_breaks']:4} second breaks")


if __name__ == "__main__":
    main()
