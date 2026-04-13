#!/usr/bin/env python3
"""
P5 DTC Survival — Cirq Implementation for Google Hardware
==========================================================

Measures <ZZ>(n) time series for n = 1..48, Fourier analysis for
subharmonic at 1/(2T). Tests DTC survival under gate noise.

Three modes: paired (merkabit), unpaired (control), perturbed.

Zero two-qubit gates. Depth 2 after optimisation (one PhXZ per qubit + measurement).

Usage:
  python run_p5_dtc_cirq.py --sim-only
  python run_p5_dtc_cirq.py --project YOUR_PROJECT --processor PROCESSOR_ID --epsilon 0.1
"""

import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import cirq

RESULTS_DIR = Path(__file__).parent.parent / "outputs" / "p5_dtc"

T_CYCLE = 12
STEP_PHASE = 2 * np.pi / T_CYCLE
GATES = ['S', 'R', 'T', 'F', 'P']


def get_gate_angles(k: int, epsilon: float = 0.0,
                    rng: np.random.Generator = None) -> tuple[float, float, float]:
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
    if epsilon > 0 and rng is not None:
        p *= (1 + epsilon * rng.standard_normal())
        rz *= (1 + epsilon * rng.standard_normal())
        rx *= (1 + epsilon * rng.standard_normal())
    return p, rz, rx


def build_dtc_circuit(n_steps: int, paired: bool = True,
                      epsilon: float = 0.0, seed: int = 42) -> cirq.Circuit:
    """DTC circuit: n ouroboros steps, optional noise, measure."""
    rng = np.random.default_rng(seed) if epsilon > 0 else None
    q_plus = cirq.LineQubit(0)
    q_minus = cirq.LineQubit(1)

    ops = []
    for k in range(n_steps):
        p, rz, rx = get_gate_angles(k, epsilon, rng)
        if paired:
            ops.append(cirq.rz(rz - p).on(q_plus))
            ops.append(cirq.rz(rz + p).on(q_minus))
        else:
            ops.append(cirq.rz(rz).on(q_plus))
            ops.append(cirq.rz(rz).on(q_minus))
        ops.append(cirq.rx(rx).on(q_plus))
        ops.append(cirq.rx(rx).on(q_minus))

    ops.append(cirq.measure(q_plus, q_minus, key='result'))
    circuit = cirq.Circuit(ops)
    circuit = cirq.optimize_for_target_gateset(
        circuit, gateset=cirq.CZTargetGateset()
    )
    return circuit


def analyze_dtc_signal(zz_series: list, T: int = 12) -> dict:
    """Fourier analysis for DTC subharmonic."""
    zz = np.array(zz_series)
    n = len(zz)
    fft = np.fft.rfft(zz)
    freqs = np.fft.rfftfreq(n, d=1.0)
    power = np.abs(fft) ** 2

    f_T = 1.0 / T
    f_2T = 1.0 / (2 * T)
    idx_T = np.argmin(np.abs(freqs - f_T))
    idx_2T = np.argmin(np.abs(freqs - f_2T))

    if n > 2 * T:
        acf_T = float(np.mean(zz[:n - T] * zz[T:n]))
        acf_2T = float(np.mean(zz[:n - 2*T] * zz[2*T:n]))
    else:
        acf_T = acf_2T = float('nan')

    return {
        "power_at_T": float(power[idx_T]),
        "power_at_2T": float(power[idx_2T]),
        "dtc_ratio": float(power[idx_2T] / (power[idx_T] + 1e-30)),
        "acf_lag_T": acf_T,
        "acf_lag_2T": acf_2T,
        "mean_zz": float(np.mean(zz)),
        "std_zz": float(np.std(zz)),
    }


def compute_ideal_zz_series(n_max, paired=True, epsilon=0.0, seed=42):
    """Compute ideal <ZZ>(n) for n=1..n_max."""
    rng = np.random.default_rng(seed) if epsilon > 0 else None
    zz_series = []
    U = np.eye(4, dtype=complex)
    state0 = np.array([1, 0, 0, 0], dtype=complex)
    ZZ = np.diag([1, -1, -1, 1])

    for k in range(n_max):
        p, rz, rx = get_gate_angles(k, epsilon, rng)
        if paired:
            Rz_f = np.diag([np.exp(-1j*(rz-p)/2), np.exp(1j*(rz-p)/2)])
            Rz_i = np.diag([np.exp(-1j*(rz+p)/2), np.exp(1j*(rz+p)/2)])
        else:
            Rz_f = np.diag([np.exp(-1j*rz/2), np.exp(1j*rz/2)])
            Rz_i = Rz_f.copy()
        c = np.cos(rx/2); s = -1j*np.sin(rx/2)
        Rx = np.array([[c, s], [s, c]])
        U = np.kron(Rx @ Rz_f, Rx @ Rz_i) @ U
        psi = U @ state0
        zz_series.append(float(np.real(psi.conj() @ ZZ @ psi)))

    return zz_series


def main():
    parser = argparse.ArgumentParser(
        description="P5 DTC survival (Cirq, for Google hardware)")
    parser.add_argument("--project", default=None)
    parser.add_argument("--processor", default=None)
    parser.add_argument("--shots", type=int, default=4096)
    parser.add_argument("--n-max", type=int, default=48)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--epsilon", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sim-only", action="store_true")
    args = parser.parse_args()

    step_list = list(range(1, args.n_max + 1, args.stride))

    print("\n== P5 DTC Survival (Cirq) ================================")
    print(f"Drive period T = {T_CYCLE} steps")
    print(f"DTC signature: subharmonic at 2T = {2*T_CYCLE} steps")
    print(f"Sweep: n = 1..{args.n_max} (stride {args.stride}, {len(step_list)} circuits/mode)")
    if args.epsilon > 0:
        print(f"Perturbation: epsilon = {args.epsilon}")
    print(f"Zero two-qubit gates.\n")

    # Ideal
    ideal_clean = compute_ideal_zz_series(args.n_max, paired=True, epsilon=0.0)
    dtc_clean = analyze_dtc_signal(ideal_clean, T_CYCLE)
    ideal_ctrl = compute_ideal_zz_series(args.n_max, paired=False, epsilon=0.0)
    dtc_ctrl = analyze_dtc_signal(ideal_ctrl, T_CYCLE)

    print(f"Ideal DTC ratios:")
    print(f"  Paired clean:    {dtc_clean['dtc_ratio']:.2f}")
    print(f"  Unpaired control: {dtc_ctrl['dtc_ratio']:.2f}")

    # Show circuit
    print(f"\nExample circuit (n=4, paired):")
    circuit = build_dtc_circuit(4, paired=True)
    print(circuit)
    print(f"Depth: {len(circuit)}")

    if args.sim_only:
        print("\n[sim-only mode]")
        return

    if not args.project or not args.processor:
        print("\nERROR: --project and --processor required for hardware.")
        return

    # Hardware
    try:
        from cirq_google import get_engine
        engine = get_engine(args.project)
        processor = engine.get_processor(args.processor)
    except ImportError:
        print("ERROR: cirq-google not installed.")
        return

    all_results = {}
    modes = [
        ("paired", True, args.epsilon),
        ("unpaired", False, 0.0),
    ]
    if args.epsilon > 0:
        modes.append(("perturbed", True, args.epsilon))

    for mode_name, paired, eps in modes:
        print(f"\n=== {mode_name} ===")
        zz_series = []
        raw_data = []
        for n in step_list:
            circuit = build_dtc_circuit(n, paired, eps, args.seed)
            job = processor.run(circuit, repetitions=args.shots)
            counts = job.histogram(key='result')
            total = sum(counts.values())
            p00 = counts.get(0, 0) / total
            p01 = counts.get(1, 0) / total
            p10 = counts.get(2, 0) / total
            p11 = counts.get(3, 0) / total
            zz = p00 - p01 - p10 + p11
            zz_series.append(zz)
            raw_data.append({
                "n": n, "zz": zz,
                "p00": p00, "p01": p01, "p10": p10, "p11": p11,
                "counts": dict(counts), "shots": total,
            })
            if n <= 3 or n % 12 == 0 or n == step_list[-1]:
                print(f"  n={n:3d}: <ZZ>={zz:+.4f}  P(00)={p00:.4f}")

        dtc = analyze_dtc_signal(zz_series, T_CYCLE)
        all_results[mode_name] = {
            "zz_series": zz_series, "raw": raw_data, "dtc": dtc
        }
        print(f"  DTC ratio: {dtc['dtc_ratio']:.2f}")

    # Summary
    print(f"\n== P5 DTC Summary (Cirq) ================================")
    for mode, data in all_results.items():
        d = data["dtc"]
        print(f"  {mode:>12}: DTC ratio = {d['dtc_ratio']:.2f}  "
              f"power(2T) = {d['power_at_2T']:.2f}")

    # Save
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = {
        "experiment": "P5_dtc_cirq",
        "platform": "google",
        "processor": args.processor,
        "results": {m: {"dtc": d["dtc"], "zz_series": d["zz_series"]}
                    for m, d in all_results.items()},
        "epsilon": args.epsilon,
        "shots": args.shots,
        "n_max": args.n_max,
        "timestamp": datetime.now().isoformat(),
    }
    path = RESULTS_DIR / f"p5_dtc_{args.processor}_{ts}.json"
    with open(path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults -> {path}")


if __name__ == "__main__":
    main()
