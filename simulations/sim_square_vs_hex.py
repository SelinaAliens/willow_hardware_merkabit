#!/usr/bin/env python3
"""
Square Grid vs Hexagonal Grid — Merkabit Simulation
=====================================================

Simulates the pentachoric ouroboros cycle on both topologies
and compares Fano factor, detection rate, and DTC signatures.

The question: does the merkabit produce sub-Poissonian statistics
on a square grid (Google Willow topology) or only on hexagonal
(IBM heavy-hex)?

Three levels of comparison:
  1. 2-qubit (topology-independent): P1b Ramsey, P2 stroboscopic, P5 DTC
     Should be identical on both — confirms the angle table doesn't care.
  2. Multi-merkabit on hex (3-node triangle, 6 neighbors max)
  3. Multi-merkabit on square (4-node square, 4 neighbors max)

Selina Stenberg with Claude Anthropic, April 2026
"""

import numpy as np
from collections import defaultdict
from itertools import product as iterproduct

# ═══════════════════════════════════════════════════════════
# GATE ANGLES (identical on both topologies)
# ═══════════════════════════════════════════════════════════

T_CYCLE = 12
STEP_PHASE = 2 * np.pi / T_CYCLE
GATES = ['S', 'R', 'T', 'F', 'P']
NUM_GATES = 5


def absent_gate(base, chirality, t):
    return (base + chirality * t) % NUM_GATES


def get_gate_angles(k, absent_idx):
    gl = GATES[absent_idx]
    p = STEP_PHASE
    sym = STEP_PHASE / 3
    w = 2 * np.pi * k / T_CYCLE
    rx = sym * (1.0 + 0.5 * np.cos(w))
    rz = sym * (1.0 + 0.5 * np.cos(w + 2 * np.pi / 3))
    if gl == 'S': rz *= 0.4; rx *= 1.3
    elif gl == 'R': rx *= 0.4; rz *= 1.3
    elif gl == 'T': rx *= 0.7; rz *= 0.7
    elif gl == 'P': p *= 0.6; rx *= 1.8; rz *= 1.5
    return p, rz, rx


def find_valid_assignment(edges, num_nodes, seed=42):
    rng = np.random.default_rng(seed)
    for _ in range(50000):
        assignment = [int(rng.integers(0, NUM_GATES)) for _ in range(num_nodes)]
        if all(assignment[i] != assignment[j] for i, j in edges):
            return assignment
    return [i % NUM_GATES for i in range(num_nodes)]


# ═══════════════════════════════════════════════════════════
# CELL DEFINITIONS
# ═══════════════════════════════════════════════════════════

class TriangleCell:
    """3-node hexagonal triangle: all three chirality classes."""
    name = "hex-triangle"
    def __init__(self):
        self.num_nodes = 3
        self.chirality = [0, +1, -1]
        self.edges = [(0, 1), (0, 2), (1, 2)]
        self.num_edges = 3
        self.coordination = [2, 2, 2]


class SquareCell4:
    """4-node square: minimal square grid (2x2)."""
    name = "square-2x2"
    def __init__(self):
        self.num_nodes = 4
        # 2x2 grid:  0 - 1
        #            |   |
        #            2 - 3
        self.edges = [(0, 1), (0, 2), (1, 3), (2, 3)]
        self.num_edges = 4
        # Checkerboard chirality (Z2, not Z3)
        self.chirality = [+1, -1, -1, +1]
        self.coordination = [2, 2, 2, 2]


class SquareCell9:
    """9-node square: 3x3 grid (same node count as triangle needs 9 qubits)."""
    name = "square-3x3"
    def __init__(self):
        self.num_nodes = 9
        # 3x3 grid:  0 - 1 - 2
        #            |   |   |
        #            3 - 4 - 5
        #            |   |   |
        #            6 - 7 - 8
        self.edges = [
            (0,1), (1,2), (3,4), (4,5), (6,7), (7,8),  # horizontal
            (0,3), (1,4), (2,5), (3,6), (4,7), (5,8),  # vertical
        ]
        self.num_edges = 12
        # Z3-like chirality on 3x3: (row+col) mod 3
        self.chirality = [
            (r + c) % 3 if (r + c) % 3 != 2 else -1
            for r in range(3) for c in range(3)
        ]
        # Map: 0->0, 1->+1, 2->-1
        self.chirality = [
            {0: 0, 1: +1, 2: -1}[(r + c) % 3]
            for r in range(3) for c in range(3)
        ]
        self.coordination = [
            2, 3, 2,
            3, 4, 3,
            2, 3, 2,
        ]


class HexCell7:
    """7-node Eisenstein cell (hex, radius 1). Centre + 6 periphery."""
    name = "hex-eisenstein-7"
    def __init__(self):
        self.num_nodes = 7
        # Standard Eisenstein cell
        unit_vecs = [(1,0),(-1,0),(0,1),(0,-1),(-1,-1),(1,1)]
        coords = []
        for a in range(-1, 2):
            for b in range(-1, 2):
                if a*a - a*b + b*b <= 1:
                    coords.append((a, b))
        self.num_nodes = len(coords)
        coord_to_idx = {c: i for i, c in enumerate(coords)}
        coord_set = set(coords)

        self.edges = []
        for i, (a1, b1) in enumerate(coords):
            for da, db in unit_vecs:
                nb = (a1 + da, b1 + db)
                if nb in coord_set:
                    j = coord_to_idx[nb]
                    if j > i:
                        self.edges.append((i, j))
        self.num_edges = len(self.edges)

        self.chirality = [
            {0: 0, 1: +1, 2: -1}[(a + b) % 3]
            for a, b in coords
        ]

        neighbours = defaultdict(list)
        for i, j in self.edges:
            neighbours[i].append(j)
            neighbours[j].append(i)
        self.coordination = [len(neighbours[i]) for i in range(self.num_nodes)]


# ═══════════════════════════════════════════════════════════
# SIMULATION ENGINE
# ═══════════════════════════════════════════════════════════

def simulate_cell(cell, tau_values=[1, 5], shots=8192, seed=42):
    """
    Classical simulation of the pentachoric circuit on a given cell.
    Returns detection rate and Fano factor at each tau.
    """
    assignment = find_valid_assignment(cell.edges, cell.num_nodes, seed=seed)
    rng = np.random.default_rng(seed)

    results = {}
    for tau in tau_values:
        # Monte Carlo: simulate syndrome measurements
        syndrome_weights = []

        for shot in range(shots):
            total_weight = 0
            per_round_weights = []

            for t in range(tau):
                round_weight = 0
                for e_idx, (i, j) in enumerate(cell.edges):
                    chi_i = cell.chirality[i]
                    chi_j = cell.chirality[j]
                    base_i = assignment[i]
                    base_j = assignment[j]
                    abs_i = absent_gate(base_i, chi_i, t)
                    abs_j = absent_gate(base_j, chi_j, t)

                    # Detection: ancilla fires if absent gates match
                    # (pentachoric closure failure)
                    if abs_i == abs_j:
                        round_weight += 1

                per_round_weights.append(round_weight)
                total_weight += round_weight

            syndrome_weights.append(total_weight)

        sw = np.array(syndrome_weights, dtype=float)
        mean_w = np.mean(sw)
        var_w = np.var(sw, ddof=1)
        fano = var_w / mean_w if mean_w > 1e-10 else float('nan')
        det_rate = np.mean(sw > 0)

        # Per-round analysis (for last tau)
        per_round_fano = []
        for t in range(tau):
            round_weights = []
            for shot in range(shots):
                rw = 0
                for e_idx, (i, j) in enumerate(cell.edges):
                    chi_i = cell.chirality[i]
                    chi_j = cell.chirality[j]
                    abs_i = absent_gate(assignment[i], chi_i, t)
                    abs_j = absent_gate(assignment[j], chi_j, t)
                    if abs_i == abs_j:
                        rw += 1
                round_weights.append(rw)
            rw_arr = np.array(round_weights, dtype=float)
            rm = np.mean(rw_arr)
            rv = np.var(rw_arr, ddof=1)
            rf = rv / rm if rm > 1e-10 else float('nan')
            per_round_fano.append(rf)

        # Edge fire rates
        edge_rates = []
        for e_idx, (i, j) in enumerate(cell.edges):
            fires = 0
            for t in range(tau):
                abs_i = absent_gate(assignment[i], cell.chirality[i], t)
                abs_j = absent_gate(assignment[j], cell.chirality[j], t)
                if abs_i == abs_j:
                    fires += 1
            edge_rates.append(fires / tau)

        chi_diffs = [abs(cell.chirality[i] - cell.chirality[j])
                     for i, j in cell.edges]

        results[f"tau_{tau}"] = {
            "detection_rate": float(det_rate),
            "fano_factor": float(fano),
            "mean_weight": float(mean_w),
            "per_round_fano": per_round_fano,
            "edge_fire_rates": edge_rates,
            "edge_chi_diffs": chi_diffs,
            "assignment": [GATES[a] for a in assignment],
        }

    return results


# ═══════════════════════════════════════════════════════════
# MAIN COMPARISON
# ═══════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("SQUARE GRID vs HEXAGONAL GRID — Merkabit Simulation")
    print("=" * 70)

    cells = [
        TriangleCell(),
        SquareCell4(),
        HexCell7(),
        SquareCell9(),
    ]

    tau_values = [1, 5, 12]

    print(f"\nCells to compare:")
    for cell in cells:
        chi_counts = {-1: 0, 0: 0, 1: 0}
        for c in cell.chirality:
            chi_counts[c] += 1
        print(f"  {cell.name}: {cell.num_nodes} nodes, {cell.num_edges} edges, "
              f"coord={cell.coordination}, "
              f"chi={{+1:{chi_counts[1]}, 0:{chi_counts[0]}, -1:{chi_counts[-1]}}}")

    all_results = {}
    for cell in cells:
        print(f"\n{'='*60}")
        print(f"  Simulating: {cell.name}")
        print(f"{'='*60}")

        results = simulate_cell(cell, tau_values=tau_values)
        all_results[cell.name] = results

        for tau in tau_values:
            k = f"tau_{tau}"
            r = results[k]
            print(f"\n  tau={tau}:")
            print(f"    Assignment: {r['assignment']}")
            print(f"    Detection rate: {r['detection_rate']:.4f}")
            print(f"    Fano factor: {r['fano_factor']:.4f} "
                  f"{'(sub-Poissonian)' if r['fano_factor'] < 1 else '(super-Poissonian)'}")
            print(f"    Mean weight: {r['mean_weight']:.4f}")
            if tau > 1:
                print(f"    Per-round Fano: {[f'{f:.4f}' for f in r['per_round_fano']]}")

            # Edge analysis by chirality difference
            counter_rates = [rate for rate, cd in zip(r['edge_fire_rates'], r['edge_chi_diffs']) if cd == 2]
            mixed_rates = [rate for rate, cd in zip(r['edge_fire_rates'], r['edge_chi_diffs']) if cd == 1]
            same_rates = [rate for rate, cd in zip(r['edge_fire_rates'], r['edge_chi_diffs']) if cd == 0]

            if counter_rates:
                print(f"    Counter-rotating edges (chi_diff=2): mean rate = {np.mean(counter_rates):.4f}")
            if mixed_rates:
                print(f"    Mixed edges (chi_diff=1): mean rate = {np.mean(mixed_rates):.4f}")
            if same_rates:
                print(f"    Same-chirality edges (chi_diff=0): mean rate = {np.mean(same_rates):.4f}")

    # ═══════════════════════════════════════════════════════
    # COMPARISON TABLE
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("COMPARISON: HEXAGONAL vs SQUARE GRID")
    print(f"{'='*70}")

    print(f"\n{'Cell':>20} | {'Topology':>8} | {'Nodes':>5} | {'Edges':>5} | "
          f"{'tau':>3} | {'Det rate':>8} | {'Fano':>8} | {'Sub-P?':>6}")
    print("-" * 85)

    for cell in cells:
        topo = "hex" if "hex" in cell.name else "square"
        for tau in tau_values:
            k = f"tau_{tau}"
            r = all_results[cell.name][k]
            sub_p = "YES" if r['fano_factor'] < 1 else "no"
            print(f"{cell.name:>20} | {topo:>8} | {cell.num_nodes:>5} | "
                  f"{cell.num_edges:>5} | {tau:>3} | "
                  f"{r['detection_rate']:>8.4f} | {r['fano_factor']:>8.4f} | {sub_p:>6}")

    # Key question
    print(f"\n{'='*70}")
    print("KEY QUESTION: Does the Fano factor differ between topologies?")
    print(f"{'='*70}")

    for tau in tau_values:
        k = f"tau_{tau}"
        hex_fano_tri = all_results["hex-triangle"][k]["fano_factor"]
        sq_fano_4 = all_results["square-2x2"][k]["fano_factor"]
        hex_fano_7 = all_results["hex-eisenstein-7"][k]["fano_factor"]
        sq_fano_9 = all_results["square-3x3"][k]["fano_factor"]

        print(f"\n  tau={tau}:")
        print(f"    hex-triangle (3n):    F = {hex_fano_tri:.4f}")
        print(f"    square-2x2 (4n):      F = {sq_fano_4:.4f}")
        print(f"    hex-eisenstein (7n):   F = {hex_fano_7:.4f}")
        print(f"    square-3x3 (9n):      F = {sq_fano_9:.4f}")

        all_sub = all(f < 1 for f in [hex_fano_tri, sq_fano_4, hex_fano_7, sq_fano_9])
        print(f"    All sub-Poissonian? {'YES' if all_sub else 'NO'}")

    print(f"\n{'='*70}")
    print("CONCLUSION")
    print(f"{'='*70}")
    # The simulation will tell us whether topology matters


if __name__ == "__main__":
    main()
