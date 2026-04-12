#!/usr/bin/env python3
"""
Square vs Hex with Noise — The Real Test
==========================================

The noiseless simulation gives Fano = 0 (deterministic). The real
question: when random errors are injected at rate epsilon, does
the pentachoric detection on a SQUARE grid produce sub-Poissonian
syndrome statistics the same way it does on a HEXAGONAL grid?

This matches what hardware does: the circuit has gate errors,
and the syndrome weight distribution's Fano factor tells us
whether the detection mechanism anti-bunches those errors.
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
    for _ in range(50000):
        a = [int(rng.integers(0, NUM_GATES)) for _ in range(num_nodes)]
        if all(a[i] != a[j] for i, j in edges):
            return a
    return [i % NUM_GATES for i in range(num_nodes)]


class TriangleCell:
    name = "hex-triangle-3n"
    def __init__(self):
        self.num_nodes = 3
        self.chirality = [0, +1, -1]
        self.edges = [(0, 1), (0, 2), (1, 2)]
        self.num_edges = 3


class SquareCell4:
    name = "square-2x2-4n"
    def __init__(self):
        self.num_nodes = 4
        self.chirality = [+1, -1, -1, +1]
        self.edges = [(0, 1), (0, 2), (1, 3), (2, 3)]
        self.num_edges = 4


class HexCell7:
    name = "hex-eisenstein-7n"
    def __init__(self):
        unit_vecs = [(1,0),(-1,0),(0,1),(0,-1),(-1,-1),(1,1)]
        coords = [(a, b) for a in range(-1, 2) for b in range(-1, 2)
                  if a*a - a*b + b*b <= 1]
        self.num_nodes = len(coords)
        idx = {c: i for i, c in enumerate(coords)}
        cs = set(coords)
        self.edges = []
        for i, (a, b) in enumerate(coords):
            for da, db in unit_vecs:
                nb = (a+da, b+db)
                if nb in cs and idx[nb] > i:
                    self.edges.append((i, idx[nb]))
        self.num_edges = len(self.edges)
        self.chirality = [{0:0,1:+1,2:-1}[(a+b)%3] for a,b in coords]


class SquareCell9:
    name = "square-3x3-9n"
    def __init__(self):
        self.num_nodes = 9
        self.edges = [(0,1),(1,2),(3,4),(4,5),(6,7),(7,8),
                      (0,3),(1,4),(2,5),(3,6),(4,7),(5,8)]
        self.num_edges = 12
        self.chirality = [{0:0,1:+1,2:-1}[(r+c)%3]
                          for r in range(3) for c in range(3)]


def simulate_noisy(cell, tau, epsilon, shots=10000, seed=42):
    """
    Monte Carlo simulation with random error injection.

    At each step t, for each node i:
      1. Compute the absent gate: absent_i = (base_i + chi_i * t) % 5
      2. With probability epsilon, flip the absent gate to a random OTHER gate
         (simulating a gate error at node i)
      3. For each edge (i,j): ancilla fires if absent_i == absent_j

    The syndrome weight per shot = total ancilla fires across all rounds.
    Fano = Var(weight) / Mean(weight).
    """
    assignment = find_valid_assignment(cell.edges, cell.num_nodes, seed=seed)
    rng = np.random.default_rng(seed + 1000)

    weights = np.zeros(shots)
    per_round_weights = np.zeros((shots, tau))

    for shot in range(shots):
        for t in range(tau):
            # Compute absent gates with possible errors
            absents = []
            for i in range(cell.num_nodes):
                ab = absent_gate(assignment[i], cell.chirality[i], t)
                if rng.random() < epsilon:
                    # Error: flip to a random DIFFERENT gate
                    ab = (ab + rng.integers(1, NUM_GATES)) % NUM_GATES
                absents.append(ab)

            # Check each edge
            round_w = 0
            for (i, j) in cell.edges:
                if absents[i] == absents[j]:
                    round_w += 1

            per_round_weights[shot, t] = round_w
            weights[shot] += round_w

    mean_w = np.mean(weights)
    var_w = np.var(weights, ddof=1)
    fano = var_w / mean_w if mean_w > 1e-10 else float('nan')
    det_rate = np.mean(weights > 0)

    # Per-round Fano
    pr_fano = []
    for t in range(tau):
        rm = np.mean(per_round_weights[:, t])
        rv = np.var(per_round_weights[:, t], ddof=1)
        pr_fano.append(rv / rm if rm > 1e-10 else float('nan'))

    return {
        "fano": fano,
        "detection_rate": det_rate,
        "mean_weight": mean_w,
        "per_round_fano": pr_fano,
        "sub_poissonian": fano < 1.0 if not np.isnan(fano) else False,
    }


def main():
    print("=" * 70)
    print("SQUARE vs HEX with NOISE — Merkabit Topology Test")
    print("=" * 70)

    cells = [TriangleCell(), SquareCell4(), HexCell7(), SquareCell9()]
    epsilons = [0.01, 0.05, 0.10, 0.20, 0.30]
    tau_values = [1, 5]
    shots = 20000

    print(f"\nError rates: {epsilons}")
    print(f"Tau values: {tau_values}")
    print(f"Shots: {shots}")

    # Run all combinations
    print(f"\n{'Cell':>20} | {'Topo':>6} | {'tau':>3} | {'eps':>5} | "
          f"{'Det':>6} | {'Fano':>8} | {'Sub-P':>5} | {'Mean wt':>8}")
    print("-" * 82)

    summary = {}
    for cell in cells:
        topo = "hex" if "hex" in cell.name else "square"
        for tau in tau_values:
            for eps in epsilons:
                r = simulate_noisy(cell, tau, eps, shots=shots)
                key = (cell.name, tau, eps)
                summary[key] = r

                sub = "YES" if r["sub_poissonian"] else "no"
                fano_str = f"{r['fano']:.4f}" if not np.isnan(r['fano']) else "nan"
                print(f"{cell.name:>20} | {topo:>6} | {tau:>3} | {eps:>5.2f} | "
                      f"{r['detection_rate']:>6.3f} | {fano_str:>8} | {sub:>5} | "
                      f"{r['mean_weight']:>8.3f}")

    # Comparison: hex vs square at matched parameters
    print(f"\n{'='*70}")
    print("DIRECT COMPARISON: Hex vs Square")
    print(f"{'='*70}")

    comparisons = [
        ("hex-triangle-3n", "square-2x2-4n", "3-node hex vs 4-node square"),
        ("hex-eisenstein-7n", "square-3x3-9n", "7-node hex vs 9-node square"),
    ]

    for hex_name, sq_name, label in comparisons:
        print(f"\n  {label}:")
        print(f"  {'tau':>3} | {'eps':>5} | {'F(hex)':>8} | {'F(sq)':>8} | "
              f"{'hex sub-P':>9} | {'sq sub-P':>9} | {'Same?':>6}")
        print("  " + "-" * 65)
        for tau in tau_values:
            for eps in epsilons:
                rh = summary[(hex_name, tau, eps)]
                rs = summary[(sq_name, tau, eps)]
                fh = rh['fano']
                fs = rs['fano']
                sh = "YES" if rh['sub_poissonian'] else "no"
                ss = "YES" if rs['sub_poissonian'] else "no"
                same = "YES" if sh == ss else "DIFF"
                print(f"  {tau:>3} | {eps:>5.2f} | {fh:>8.4f} | {fs:>8.4f} | "
                      f"{sh:>9} | {ss:>9} | {same:>6}")

    # Final verdict
    print(f"\n{'='*70}")
    print("VERDICT")
    print(f"{'='*70}")

    hex_sub_p = sum(1 for k, v in summary.items()
                    if "hex" in k[0] and v["sub_poissonian"])
    sq_sub_p = sum(1 for k, v in summary.items()
                   if "square" in k[0] and v["sub_poissonian"])
    hex_total = sum(1 for k in summary if "hex" in k[0])
    sq_total = sum(1 for k in summary if "square" in k[0])

    print(f"\n  Hex:    {hex_sub_p}/{hex_total} runs sub-Poissonian")
    print(f"  Square: {sq_sub_p}/{sq_total} runs sub-Poissonian")

    if sq_sub_p == sq_total and hex_sub_p == hex_total:
        print(f"\n  PREDICTION: The merkabit WILL produce sub-Poissonian")
        print(f"  statistics on Google Willow (square grid).")
        print(f"  The effect is topology-independent in simulation.")
    elif sq_sub_p > 0:
        print(f"\n  PREDICTION: Partial — some conditions produce sub-P")
        print(f"  on square grid. Topology matters at specific (tau, eps).")
    else:
        print(f"\n  PREDICTION: The merkabit does NOT produce sub-P on")
        print(f"  square grid. The effect is topology-dependent.")


if __name__ == "__main__":
    main()
