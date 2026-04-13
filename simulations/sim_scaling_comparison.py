#!/usr/bin/env python3
"""
Lattice Scaling Simulation: Hex vs Square
==========================================

How does the pentachoric merkabit scale with lattice size on both
hexagonal and square grid topologies?

Hex scaling:   3-node triangle → 7-node Eisenstein → 19-node (r=2)
Square scaling: 4-node (2x2) → 9-node (3x3) → 16-node (4x4) → 25-node (5x5)

At each cell size, measures:
  - Fano factor vs error rate (epsilon sweep)
  - Detection rate vs tau
  - Interior vs boundary Fano separation
  - Paired vs control Fano ratio (the merkabit advantage)

The key prediction for Paper 26: does the square grid scale as well
as (or better than) hex, supporting the claim that the merkabit
is geometric rather than topological?
"""

import numpy as np
from collections import defaultdict

T_CYCLE = 12
STEP_PHASE = 2 * np.pi / T_CYCLE
GATES = ['S', 'R', 'T', 'F', 'P']
NUM_GATES = 5


def absent_gate(base, chirality, t):
    return (base + chirality * t) % NUM_GATES


def find_valid_assignment(edges, num_nodes, seed=42):
    rng = np.random.default_rng(seed)
    for _ in range(100000):
        a = [int(rng.integers(0, NUM_GATES)) for _ in range(num_nodes)]
        if all(a[i] != a[j] for i, j in edges):
            return a
    # Fallback: greedy
    a = [0] * num_nodes
    nbrs = defaultdict(set)
    for i, j in edges:
        nbrs[i].add(j)
        nbrs[j].add(i)
    for i in range(num_nodes):
        used = {a[j] for j in nbrs[i] if j < i}
        for g in range(NUM_GATES):
            if g not in used:
                a[i] = g
                break
    return a


# ═══════════════════════════════════════════════════════════
# CELL BUILDERS
# ═══════════════════════════════════════════════════════════

def build_hex_cell(radius):
    """Eisenstein cell of given radius. r=0: 1 node, r=1: 7, r=2: 19, r=3: 37."""
    unit_vecs = [(1,0),(-1,0),(0,1),(0,-1),(-1,-1),(1,1)]
    coords = [(a, b) for a in range(-radius-1, radius+2)
              for b in range(-radius-1, radius+2)
              if a*a - a*b + b*b <= radius*radius]
    num_nodes = len(coords)
    idx = {c: i for i, c in enumerate(coords)}
    cs = set(coords)
    edges = []
    for i, (a, b) in enumerate(coords):
        for da, db in unit_vecs:
            nb = (a+da, b+db)
            if nb in cs and idx[nb] > i:
                edges.append((i, idx[nb]))
    chirality = [{0:0, 1:+1, 2:-1}[(a+b) % 3] for a, b in coords]
    # Coordination
    nbr_count = defaultdict(int)
    for i, j in edges:
        nbr_count[i] += 1
        nbr_count[j] += 1
    coord = [nbr_count[i] for i in range(num_nodes)]
    max_coord = max(coord) if coord else 0
    interior = [i for i in range(num_nodes) if coord[i] == max_coord]
    boundary = [i for i in range(num_nodes) if coord[i] < max_coord]
    return {
        "name": f"hex-r{radius}-{num_nodes}n",
        "topo": "hex",
        "num_nodes": num_nodes,
        "num_edges": len(edges),
        "edges": edges,
        "chirality": chirality,
        "coordination": coord,
        "interior": interior,
        "boundary": boundary,
    }


def build_hex_triangle():
    """Special case: 3-node triangle (not a radius cell)."""
    return {
        "name": "hex-tri-3n",
        "topo": "hex",
        "num_nodes": 3,
        "num_edges": 3,
        "edges": [(0,1), (0,2), (1,2)],
        "chirality": [0, +1, -1],
        "coordination": [2, 2, 2],
        "interior": [],
        "boundary": [0, 1, 2],
    }


def build_square_cell(side):
    """Square grid of side x side nodes."""
    n = side * side
    edges = []
    for r in range(side):
        for c in range(side):
            i = r * side + c
            if c < side - 1:
                edges.append((i, i + 1))
            if r < side - 1:
                edges.append((i, i + side))
    chirality = [{0:0, 1:+1, 2:-1}[(r+c) % 3] for r in range(side) for c in range(side)]
    # Coordination
    nbr_count = defaultdict(int)
    for i, j in edges:
        nbr_count[i] += 1
        nbr_count[j] += 1
    coord = [nbr_count[i] for i in range(n)]
    max_coord = max(coord) if coord else 0
    interior = [i for i in range(n) if coord[i] == max_coord]
    boundary = [i for i in range(n) if coord[i] < max_coord]
    return {
        "name": f"sq-{side}x{side}-{n}n",
        "topo": "square",
        "num_nodes": n,
        "num_edges": len(edges),
        "edges": edges,
        "chirality": chirality,
        "coordination": coord,
        "interior": interior,
        "boundary": boundary,
    }


# ═══════════════════════════════════════════════════════════
# SIMULATION ENGINE
# ═══════════════════════════════════════════════════════════

def simulate(cell, tau, epsilon, shots=10000, seed=42, paired=True):
    """Monte Carlo simulation with noise injection. Returns Fano + detection."""
    assignment = find_valid_assignment(cell["edges"], cell["num_nodes"], seed=seed)
    rng = np.random.default_rng(seed + 2000)

    weights = np.zeros(shots)
    for shot in range(shots):
        for t in range(tau):
            absents = []
            for i in range(cell["num_nodes"]):
                ab = absent_gate(assignment[i], cell["chirality"][i], t)
                if not paired:
                    # Control: remove chirality effect (all use chi=0)
                    ab = absent_gate(assignment[i], 0, t)
                if rng.random() < epsilon:
                    ab = (ab + rng.integers(1, NUM_GATES)) % NUM_GATES
                absents.append(ab)
            for (i, j) in cell["edges"]:
                if absents[i] == absents[j]:
                    weights[shot] += 1

    mean_w = np.mean(weights)
    var_w = np.var(weights, ddof=1)
    fano = var_w / mean_w if mean_w > 1e-10 else float('nan')
    det_rate = np.mean(weights > 0)
    return {"fano": fano, "detection": det_rate, "mean_weight": mean_w}


# ═══════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("LATTICE SCALING SIMULATION: HEX vs SQUARE")
    print("=" * 70)

    # Build all cells
    hex_cells = [
        build_hex_triangle(),
        build_hex_cell(1),   # 7 nodes
        build_hex_cell(2),   # 19 nodes
    ]
    sq_cells = [
        build_square_cell(2),  # 4 nodes
        build_square_cell(3),  # 9 nodes
        build_square_cell(4),  # 16 nodes
        build_square_cell(5),  # 25 nodes
    ]
    all_cells = hex_cells + sq_cells

    print(f"\nCells:")
    for cell in all_cells:
        int_frac = len(cell["interior"]) / cell["num_nodes"] * 100 if cell["num_nodes"] > 0 else 0
        print(f"  {cell['name']:>16}: {cell['num_nodes']:3d} nodes, {cell['num_edges']:3d} edges, "
              f"interior: {len(cell['interior']):2d}/{cell['num_nodes']} ({int_frac:.0f}%), "
              f"coord range: {min(cell['coordination'])}-{max(cell['coordination'])}")

    # ── Test 1: Fano vs cell size at fixed epsilon and tau ──
    print(f"\n{'='*70}")
    print("TEST 1: Fano vs cell size (tau=5, epsilon=0.10)")
    print(f"{'='*70}")
    eps = 0.10
    tau = 5
    shots = 15000

    print(f"\n{'Cell':>16} | {'Topo':>6} | {'Nodes':>5} | {'Edges':>5} | "
          f"{'F(paired)':>10} | {'F(ctrl)':>8} | {'Ratio':>6} | {'Det':>5}")
    print("-" * 80)

    scaling_data = []
    for cell in all_cells:
        r_p = simulate(cell, tau, eps, shots, paired=True)
        r_c = simulate(cell, tau, eps, shots, paired=False)
        ratio = r_p["fano"] / r_c["fano"] if r_c["fano"] > 0.01 else float('nan')
        print(f"{cell['name']:>16} | {cell['topo']:>6} | {cell['num_nodes']:>5} | "
              f"{cell['num_edges']:>5} | {r_p['fano']:>10.4f} | {r_c['fano']:>8.4f} | "
              f"{ratio:>6.3f} | {r_p['detection']:>5.3f}")
        scaling_data.append({
            "cell": cell["name"], "topo": cell["topo"],
            "nodes": cell["num_nodes"], "edges": cell["num_edges"],
            "fano_paired": r_p["fano"], "fano_ctrl": r_c["fano"],
            "ratio": ratio, "detection": r_p["detection"],
            "interior_frac": len(cell["interior"]) / cell["num_nodes"],
        })

    # ── Test 2: Fano vs epsilon (scaling with noise) ──
    print(f"\n{'='*70}")
    print("TEST 2: Fano vs epsilon at tau=5 (noise scaling)")
    print(f"{'='*70}")
    epsilons = [0.01, 0.05, 0.10, 0.20, 0.30]
    tau = 5

    # Pick representative cells
    rep_cells = [hex_cells[1], sq_cells[1], hex_cells[2], sq_cells[2]]
    print(f"\n{'Cell':>16} | " + " | ".join(f"eps={e:.2f}" for e in epsilons))
    print("-" * (20 + 11 * len(epsilons)))

    for cell in rep_cells:
        fanos = []
        for eps in epsilons:
            r = simulate(cell, tau, eps, shots=10000, paired=True)
            fanos.append(r["fano"])
        print(f"{cell['name']:>16} | " +
              " | ".join(f"{f:>8.4f}" for f in fanos))

    # ── Test 3: Fano vs tau (depth scaling) ──
    print(f"\n{'='*70}")
    print("TEST 3: Fano vs tau at epsilon=0.10 (depth scaling)")
    print(f"{'='*70}")
    tau_values = [1, 3, 5, 8, 12]
    eps = 0.10

    print(f"\n{'Cell':>16} | " + " | ".join(f"tau={t:>2}" for t in tau_values))
    print("-" * (20 + 10 * len(tau_values)))

    for cell in rep_cells:
        fanos = []
        for tau in tau_values:
            r = simulate(cell, tau, eps, shots=10000, paired=True)
            fanos.append(r["fano"])
        print(f"{cell['name']:>16} | " +
              " | ".join(f"{f:>7.4f}" for f in fanos))

    # ── Test 4: Interior vs boundary Fano ──
    print(f"\n{'='*70}")
    print("TEST 4: Interior vs boundary detection rate (tau=5, eps=0.10)")
    print(f"{'='*70}")

    tau = 5; eps = 0.10
    for cell in all_cells:
        if len(cell["interior"]) == 0:
            continue
        assignment = find_valid_assignment(cell["edges"], cell["num_nodes"], seed=42)
        rng = np.random.default_rng(2042)

        # Track per-node fire counts
        node_fires = np.zeros(cell["num_nodes"])
        total_rounds = 0
        n_shots = 10000

        for shot in range(n_shots):
            for t in range(tau):
                absents = []
                for i in range(cell["num_nodes"]):
                    ab = absent_gate(assignment[i], cell["chirality"][i], t)
                    if rng.random() < eps:
                        ab = (ab + rng.integers(1, NUM_GATES)) % NUM_GATES
                    absents.append(ab)
                for (i, j) in cell["edges"]:
                    if absents[i] == absents[j]:
                        node_fires[i] += 1
                        node_fires[j] += 1
                total_rounds += 1

        int_rate = np.mean(node_fires[cell["interior"]]) / (total_rounds) if cell["interior"] else float('nan')
        bnd_rate = np.mean(node_fires[cell["boundary"]]) / (total_rounds) if cell["boundary"] else float('nan')
        ratio = int_rate / bnd_rate if bnd_rate > 0 and not np.isnan(bnd_rate) else float('nan')
        int_coord = cell['coordination'][cell['interior'][0]] if cell['interior'] else '-'
        bnd_coord = cell['coordination'][cell['boundary'][0]] if cell['boundary'] else '-'
        print(f"  {cell['name']:>16}: interior rate={int_rate:.4f}  boundary={bnd_rate:.4f}  "
              f"ratio={ratio:.3f}  coord: int={int_coord} bnd={bnd_coord}")

    # ── Summary ──
    print(f"\n{'='*70}")
    print("SCALING SUMMARY")
    print(f"{'='*70}")

    print(f"\nFano at tau=5, eps=0.10:")
    print(f"{'Cell':>16} | {'Nodes':>5} | {'Topo':>6} | {'F(paired)':>10} | "
          f"{'F(ctrl)':>8} | {'P/C ratio':>9} | {'Int frac':>8}")
    print("-" * 75)
    for d in scaling_data:
        print(f"{d['cell']:>16} | {d['nodes']:>5} | {d['topo']:>6} | "
              f"{d['fano_paired']:>10.4f} | {d['fano_ctrl']:>8.4f} | "
              f"{d['ratio']:>9.3f} | {d['interior_frac']:>8.1%}")

    # Key comparison: does square scale as well as hex?
    hex_fanos = [(d['nodes'], d['fano_paired']) for d in scaling_data if d['topo'] == 'hex']
    sq_fanos = [(d['nodes'], d['fano_paired']) for d in scaling_data if d['topo'] == 'square']

    print(f"\nScaling trend (Fano vs nodes, tau=5, eps=0.10):")
    print(f"  Hex:    {' -> '.join(f'{n}n: F={f:.3f}' for n, f in hex_fanos)}")
    print(f"  Square: {' -> '.join(f'{n}n: F={f:.3f}' for n, f in sq_fanos)}")

    hex_improving = all(hex_fanos[i][1] >= hex_fanos[i+1][1]
                        for i in range(len(hex_fanos)-1)) if len(hex_fanos) > 1 else True
    sq_improving = all(sq_fanos[i][1] >= sq_fanos[i+1][1]
                       for i in range(len(sq_fanos)-1)) if len(sq_fanos) > 1 else True

    print(f"\n  Hex Fano improving with size? {hex_improving}")
    print(f"  Square Fano improving with size? {sq_improving}")

    if sq_improving:
        print(f"\n  PREDICTION: The merkabit scales on square grids.")
        print(f"  The anti-bunching strengthens with lattice size on BOTH topologies.")
    else:
        print(f"\n  The square grid does NOT show monotonic improvement.")
        print(f"  Scaling may be topology-dependent at larger sizes.")


if __name__ == "__main__":
    main()
