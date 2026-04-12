#!/usr/bin/env python3
"""
P1b Ramsey Berry Phase — Cirq Implementation for Google Hardware
=================================================================

Identical protocol to the IBM Qiskit version (merkabit_hardware_test/
experiments/run_p1_ramsey.py), translated to Google's Cirq framework.

Circuit: Ry(pi/2) -> U0_n -> Ry(-pi/2) -> measure
Forward vs reversed chirality (P sign swap) at each n.
Berry phase = (delta_fwd - delta_rev) / 2.

Zero two-qubit gates. All single-qubit. Depth 6 after compilation
(Cirq's optimiser merges consecutive single-qubit gates).

Usage:
  # Simulation only (no hardware)
  python run_p1b_ramsey_cirq.py --sim-only

  # On Google hardware (requires google.cloud.quantum access)
  python run_p1b_ramsey_cirq.py --project YOUR_PROJECT --processor PROCESSOR_ID

Authors: Stenberg & Hetland with Claude Anthropic, April 2026
"""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import cirq

RESULTS_DIR = Path(__file__).parent.parent / "outputs" / "p1_ramsey"

T_CYCLE = 12
STEP_PHASE = 2 * np.pi / T_CYCLE
GATES = ['S', 'R', 'T', 'F', 'P']


def get_gate_angles(k: int) -> tuple[float, float, float]:
    absent = k % 5
    gate_label = GATES[absent]
    p = STEP_PHASE
    sym = STEP_PHASE / 3
    w = 2 * np.pi * k / T_CYCLE
    rx = sym * (1.0 + 0.5 * np.cos(w))
    rz = sym * (1.0 + 0.5 * np.cos(w + 2 * np.pi / 3))
    if gate_label == 'S': rz *= 0.4;  rx *= 1.3
    elif gate_label == 'R': rx *= 0.4; rz *= 1.3
    elif gate_label == 'T': rx *= 0.7; rz *= 0.7
    elif gate_label == 'P': p  *= 0.6; rx *= 1.8; rz *= 1.5
    return p, rz, rx


# --- Circuit builders --------------------------------------------------------

def build_ramsey_circuit(n_steps: int, forward: bool = True) -> cirq.Circuit:
    """
    Ramsey interferometry on the merkabit dual-spinor pair.

    Ry(pi/2) -> U0_n -> Ry(-pi/2) -> measure

    q_plus = forward spinor, q_minus = inverse spinor.
    forward=True: P gate signs Rz(-p) on q+, Rz(+p) on q-
    forward=False: reversed chirality (swap P signs)
    """
    q_plus = cirq.LineQubit(0)
    q_minus = cirq.LineQubit(1)

    ops = []

    # Prepare |+> state: Ry(pi/2) on both qubits
    ops.append(cirq.ry(np.pi / 2).on(q_plus))
    ops.append(cirq.ry(np.pi / 2).on(q_minus))

    # Apply n ouroboros steps
    for k in range(n_steps):
        p, rz, rx = get_gate_angles(k)
        if forward:
            ops.append(cirq.rz(rz - p).on(q_plus))
            ops.append(cirq.rz(rz + p).on(q_minus))
        else:
            # Reversed chirality: swap P sign
            ops.append(cirq.rz(rz + p).on(q_plus))
            ops.append(cirq.rz(rz - p).on(q_minus))
        ops.append(cirq.rx(rx).on(q_plus))
        ops.append(cirq.rx(rx).on(q_minus))

    # Close Ramsey: Ry(-pi/2) on both
    ops.append(cirq.ry(-np.pi / 2).on(q_plus))
    ops.append(cirq.ry(-np.pi / 2).on(q_minus))

    # Measure
    ops.append(cirq.measure(q_plus, q_minus, key='result'))

    circuit = cirq.Circuit(ops)

    # Optimise: merge consecutive single-qubit gates
    circuit = cirq.optimize_for_target_gateset(
        circuit, gateset=cirq.CZTargetGateset()
    )

    return circuit


# --- Ideal simulation ---------------------------------------------------------

def simulate_ramsey(n_steps: int, forward: bool = True) -> dict:
    """Ideal simulation via Cirq statevector simulator."""
    q_plus = cirq.LineQubit(0)
    q_minus = cirq.LineQubit(1)

    # Build circuit without measurement for statevector
    ops = []
    ops.append(cirq.ry(np.pi / 2).on(q_plus))
    ops.append(cirq.ry(np.pi / 2).on(q_minus))
    for k in range(n_steps):
        p, rz, rx = get_gate_angles(k)
        if forward:
            ops.append(cirq.rz(rz - p).on(q_plus))
            ops.append(cirq.rz(rz + p).on(q_minus))
        else:
            ops.append(cirq.rz(rz + p).on(q_plus))
            ops.append(cirq.rz(rz - p).on(q_minus))
        ops.append(cirq.rx(rx).on(q_plus))
        ops.append(cirq.rx(rx).on(q_minus))
    ops.append(cirq.ry(-np.pi / 2).on(q_plus))
    ops.append(cirq.ry(-np.pi / 2).on(q_minus))

    circuit = cirq.Circuit(ops)
    result = cirq.Simulator().simulate(circuit)
    state = result.final_state_vector

    # Extract probabilities
    probs = np.abs(state) ** 2
    p00, p01, p10, p11 = probs[0], probs[1], probs[2], probs[3]

    z_fwd = float((p00 + p01) - (p10 + p11))
    z_inv = float((p00 + p10) - (p01 + p11))
    zz = float(p00 - p01 - p10 + p11)

    return {"z_fwd": z_fwd, "z_inv": z_inv, "zz": zz,
            "probs": {"00": float(p00), "01": float(p01),
                      "10": float(p10), "11": float(p11)}}


# --- Hardware runner ----------------------------------------------------------

def run_on_hardware(processor_id, project_id, n_steps_list, shots):
    """Run on Google quantum hardware via Cloud."""
    try:
        import google.cloud.quantum_v1alpha1 as qclient
        from cirq_google import get_engine
        engine = get_engine(project_id)
        processor = engine.get_processor(processor_id)
    except ImportError:
        print("ERROR: google-cloud-quantum not installed.")
        print("Install: pip install cirq-google google-cloud-quantum")
        return None

    results = []
    for n in n_steps_list:
        entry = {"n_steps": n}
        for direction, forward in [("fwd", True), ("rev", False)]:
            circuit = build_ramsey_circuit(n, forward)
            job = processor.run(circuit, repetitions=shots)
            counts = job.histogram(key='result')
            total = sum(counts.values())
            probs = {}
            for bitval in range(4):
                key = f"{bitval:02b}"
                probs[key] = counts.get(bitval, 0) / total
            p00 = probs.get("00", 0)
            p01 = probs.get("01", 0)
            p10 = probs.get("10", 0)
            p11 = probs.get("11", 0)
            z_fwd = (p00 + p01) - (p10 + p11)
            z_inv = (p00 + p10) - (p01 + p11)
            zz = p00 - p01 - p10 + p11
            entry[direction] = {
                "z_fwd": z_fwd, "z_inv": z_inv, "zz": zz,
                "probs": probs, "counts": dict(counts),
            }
        entry["z_fwd_diff"] = entry["fwd"]["z_fwd"] - entry["rev"]["z_fwd"]
        entry["z_inv_diff"] = entry["fwd"]["z_inv"] - entry["rev"]["z_inv"]
        entry["zz_diff"] = entry["fwd"]["zz"] - entry["rev"]["zz"]
        results.append(entry)
        print(f"  n={n}: z_fwd_diff={entry['z_fwd_diff']:+.4f}  "
              f"z_inv_diff={entry['z_inv_diff']:+.4f}")
    return results


# --- Main ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="P1b Ramsey Berry phase (Cirq, for Google hardware)")
    parser.add_argument("--project", default=None, help="Google Cloud project ID")
    parser.add_argument("--processor", default=None, help="Google processor ID")
    parser.add_argument("--shots", type=int, default=4096)
    parser.add_argument("--steps", nargs='+', type=int,
                        default=[1, 2, 3, 4, 6, 8, 10, 12])
    parser.add_argument("--sim-only", action="store_true")
    args = parser.parse_args()

    print("\n== P1b Ramsey Berry Phase (Cirq) =========================")
    print("Protocol: forward vs reversed chirality Ramsey interferometry")
    print("Berry phase = (delta_fwd - delta_rev) / 2")
    print(f"Steps: {args.steps}")
    print(f"Zero two-qubit gates. All single-qubit.\n")

    # Ideal predictions
    print("Ideal Ramsey signals (Cirq statevector simulator):")
    print(f"{'n':>3} | {'<Z+>_fwd':>9} | {'<Z+>_rev':>9} | {'diff':>7} | {'<ZZ>_fwd':>9}")
    print("-" * 55)
    for n in args.steps:
        sim_f = simulate_ramsey(n, forward=True)
        sim_r = simulate_ramsey(n, forward=False)
        diff = sim_f['z_fwd'] - sim_r['z_fwd']
        print(f"  {n:2d} | {sim_f['z_fwd']:>+9.4f} | {sim_r['z_fwd']:>+9.4f} | "
              f"{diff:>+7.4f} | {sim_f['zz']:>+9.4f}")

    # Show a circuit for inspection
    print(f"\nExample circuit (n=4, forward):")
    circuit = build_ramsey_circuit(4, forward=True)
    print(circuit)
    print(f"Depth: {len(circuit)}")
    print(f"Two-qubit gates: {sum(1 for op in circuit.all_operations() if len(op.qubits) > 1)}")

    if args.sim_only:
        print("\n[sim-only mode — no hardware submission]")
        return

    # Hardware
    if not args.project or not args.processor:
        print("\nERROR: --project and --processor required for hardware runs.")
        print("Example: python run_p1b_ramsey_cirq.py --project my-project --processor rainbow")
        return

    print(f"\nProcessor: {args.processor} (project: {args.project})")
    sweep = run_on_hardware(args.processor, args.project, args.steps, args.shots)

    if sweep:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        output = {
            "experiment": "P1b_ramsey_cirq",
            "platform": "google",
            "processor": args.processor,
            "project": args.project,
            "sweep": sweep,
            "shots": args.shots,
            "timestamp": datetime.now().isoformat(),
        }
        path = RESULTS_DIR / f"p1b_ramsey_{args.processor}_{ts}.json"
        with open(path, "w") as f:
            json.dump(output, f, indent=2, default=str)
        print(f"\nResults -> {path}")


if __name__ == "__main__":
    main()
