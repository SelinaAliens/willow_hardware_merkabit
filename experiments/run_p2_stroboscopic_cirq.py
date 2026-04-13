#!/usr/bin/env python3
"""
P2 Stroboscopic Quasi-Period — Cirq Implementation for Google Hardware
=======================================================================

Measures P(|00>) return probability across n = 1..59 ouroboros steps.
Predicts quasi-period peak at 3.3T = 39.6 steps.

Zero two-qubit gates. Depth 2 after optimisation (one PhXZ per qubit + measurement).

Usage:
  python run_p2_stroboscopic_cirq.py --sim-only
  python run_p2_stroboscopic_cirq.py --project YOUR_PROJECT --processor PROCESSOR_ID
"""

import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import cirq

RESULTS_DIR = Path(__file__).parent.parent / "outputs" / "p2_stroboscopic"

T_CYCLE = 12
STEP_PHASE = 2 * np.pi / T_CYCLE
GATES = ['S', 'R', 'T', 'F', 'P']
QUASI_PERIOD = 3.3 * T_CYCLE


def get_gate_angles(k: int) -> tuple[float, float, float]:
    absent = k % 5
    gl = GATES[absent]
    p = STEP_PHASE; sym = STEP_PHASE / 3
    w = 2 * np.pi * k / T_CYCLE
    rx = sym * (1.0 + 0.5 * np.cos(w))
    rz = sym * (1.0 + 0.5 * np.cos(w + 2 * np.pi / 3))
    if gl == 'S': rz *= 0.4; rx *= 1.3
    elif gl == 'R': rx *= 0.4; rz *= 1.3
    elif gl == 'T': rx *= 0.7; rz *= 0.7
    elif gl == 'P': p *= 0.6; rx *= 1.8; rz *= 1.5
    return p, rz, rx


def build_stroboscopic_circuit(n_steps: int) -> cirq.Circuit:
    """n steps of U0, then measure both qubits."""
    q_plus = cirq.LineQubit(0)
    q_minus = cirq.LineQubit(1)

    ops = []
    for k in range(n_steps):
        p, rz, rx = get_gate_angles(k)
        ops.append(cirq.rz(rz - p).on(q_plus))
        ops.append(cirq.rz(rz + p).on(q_minus))
        ops.append(cirq.rx(rx).on(q_plus))
        ops.append(cirq.rx(rx).on(q_minus))

    ops.append(cirq.measure(q_plus, q_minus, key='result'))

    circuit = cirq.Circuit(ops)
    circuit = cirq.optimize_for_target_gateset(
        circuit, gateset=cirq.CZTargetGateset()
    )
    return circuit


def compute_ideal(n_steps: int) -> dict:
    """Compute ideal P(|00>) via statevector."""
    q_plus = cirq.LineQubit(0)
    q_minus = cirq.LineQubit(1)
    ops = []
    for k in range(n_steps):
        p, rz, rx = get_gate_angles(k)
        ops.append(cirq.rz(rz - p).on(q_plus))
        ops.append(cirq.rz(rz + p).on(q_minus))
        ops.append(cirq.rx(rx).on(q_plus))
        ops.append(cirq.rx(rx).on(q_minus))

    circuit = cirq.Circuit(ops)
    result = cirq.Simulator().simulate(circuit)
    probs = np.abs(result.final_state_vector) ** 2
    return {"p_return": float(probs[0]), "n": n_steps}


def main():
    parser = argparse.ArgumentParser(
        description="P2 Stroboscopic quasi-period (Cirq, for Google hardware)")
    parser.add_argument("--project", default=None)
    parser.add_argument("--processor", default=None)
    parser.add_argument("--shots", type=int, default=4096)
    parser.add_argument("--n-max", type=int, default=60)
    parser.add_argument("--stride", type=int, default=2)
    parser.add_argument("--sim-only", action="store_true")
    args = parser.parse_args()

    step_list = list(range(1, args.n_max + 1, args.stride))

    print("\n== P2 Stroboscopic Quasi-Period (Cirq) ==================")
    print(f"Prediction: quasi-period = 3.3T = {QUASI_PERIOD:.1f} steps")
    print(f"Sweep: n = 1..{args.n_max} (stride {args.stride}, {len(step_list)} circuits)")
    print(f"Zero two-qubit gates.\n")

    # Ideal
    print(f"{'n':>3} | {'n/T':>5} | {'P_ideal':>8} | Note")
    print("-" * 45)
    for n in step_list:
        ideal = compute_ideal(n)
        note = ""
        if abs(n - QUASI_PERIOD) < 2: note = "<-- 3.3T"
        if n % 12 == 0: note = f"<-- {n//12}T"
        print(f"{n:>3} | {n/12:>5.1f} | {ideal['p_return']:>8.4f} | {note}")

    # Show circuit
    print(f"\nExample circuit (n=4):")
    circuit = build_stroboscopic_circuit(4)
    print(circuit)
    print(f"Depth: {len(circuit)}")

    if args.sim_only:
        print("\n[sim-only mode]")
        return

    if not args.project or not args.processor:
        print("\nERROR: --project and --processor required for hardware.")
        return

    # Hardware (batch submission)
    try:
        from cirq_google import get_engine
        engine = get_engine(args.project)
        processor = engine.get_processor(args.processor)
    except ImportError:
        print("ERROR: cirq-google not installed.")
        return

    print(f"\nBuilding {len(step_list)} circuits...")
    circuits = [build_stroboscopic_circuit(n) for n in step_list]

    results = []
    for i, (n, circuit) in enumerate(zip(step_list, circuits)):
        job = processor.run(circuit, repetitions=args.shots)
        counts = job.histogram(key='result')
        total = sum(counts.values())
        p00 = counts.get(0, 0) / total
        ideal = compute_ideal(n)
        fid = p00 / ideal["p_return"] * 100 if ideal["p_return"] > 0 else 0
        results.append({
            "n_steps": n,
            "p_return_hw": p00,
            "p_return_ideal": ideal["p_return"],
            "fidelity_pct": fid,
            "counts": dict(counts),
            "shots": total,
        })
        print(f"  n={n:3d}: P_hw={p00:.4f}  P_ideal={ideal['p_return']:.4f}  fid={fid:.1f}%")

    # Save
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = {
        "experiment": "P2_stroboscopic_cirq",
        "platform": "google",
        "processor": args.processor,
        "sweep": results,
        "shots": args.shots,
        "timestamp": datetime.now().isoformat(),
    }
    path = RESULTS_DIR / f"p2_strobo_{args.processor}_{ts}.json"
    with open(path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults -> {path}")


if __name__ == "__main__":
    main()
