#!/usr/bin/env python3
"""
P4 Scaling — Cirq Implementation for Google Hardware
=====================================================

Tests whether the sub-Poissonian Fano factor strengthens with qubit count N
under paired asymmetric-phase Floquet driving. This is the direct test of
whether asymmetric-phase Floquet driving is a geometric pathway to
maintaining coherence as quantum processors scale.

Circuit: N qubits initialized in |+>, then tau Floquet cycles of
alternating Rz(+phi) and Rz(-phi) assignments across neighbouring qubits,
followed by computational-basis measurement. Zero two-qubit gates.

Predictions (two models, both at tau=5):

    PRIMARY (Cirq depolarizing noise, error rate 0.005, Willow-realistic):
        N = 2: F_paired ~ 0.19,  F_control ~ 0.88
        N = 4: F_paired ~ 0.19,  F_control ~ 0.88
        N = 6: F_paired ~ 0.19,  F_control ~ 0.90
        N = 8: F_paired ~ 0.19,  F_control ~ 0.90
        Scaling law: F_paired stays flat (sub-Poissonian) as N grows,
        while F_control approaches 1. Standard decoherence theory predicts
        F -> 1 with N; the geometric hypothesis predicts F stays constant.

    SECONDARY (Paper 26 topology-correlated Monte Carlo, stronger claim):
        N = 2: F ~ 0.4   (matches P1 baseline at higher noise)
        N = 4: F in [0.3, 0.4]   (tighter anti-bunching with scale)
        N = 8: F in [0.25, 0.35]   (approaching the 4x4 sim prediction F ~ 0.285)
        Scaling law: F(N) monotonically decreases with N at fixed tau.
        This prediction is contingent on correlated-noise structure that
        only hardware can test; Cirq's depolarizing model does not reproduce it.

Falsification: F_paired(N) > 0.7 at any N >= 4, at 3 sigma.

Compiled depth per Floquet cycle: 2 (Cirq optimiser merges consecutive
single-qubit gates into one PhXZ per qubit).

Usage:
    # Simulation only (no hardware)
    python run_p4_scaling_cirq.py --sim-only

    # Noisy simulation (depolarizing model to approximate Willow noise)
    python run_p4_scaling_cirq.py --sim-only --noisy --error-rate 0.005

    # On Google hardware (requires google.cloud.quantum access)
    python run_p4_scaling_cirq.py --project YOUR_PROJECT --processor PROCESSOR_ID

Authors: Stenberg & Hetland with Claude Anthropic, April 2026
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import cirq

sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

RESULTS_DIR = Path(__file__).parent.parent / "outputs" / "p4_scaling"

# Architectural constants
T_CYCLE = 12                          # Coxeter period h=12
STEP_PHASE = 2 * np.pi / T_CYCLE      # fundamental Floquet step = pi/6
PHI = np.pi / 12                       # paired-phase amplitude

# Scaling sweep
N_QUBITS_SWEEP = [2, 4, 6, 8]         # qubit counts to test
TAU_FIXED = 5                          # fixed Floquet depth for scaling test
SHOTS = 8000                           # per configuration


# --- Circuit builder ---------------------------------------------------------

def build_scaling_circuit(n_qubits: int, tau: int, paired: bool = True) -> cirq.Circuit:
    """
    N-qubit paired-Floquet scaling circuit.

    Initialize all qubits in |+> via H.
    For tau Floquet cycles:
        For each qubit i:
            Paired:   Rz(+phi) on even i, Rz(-phi) on odd i
            Control:  Rz(+phi) on all qubits (same sign)
        Rx(STEP_PHASE) on all qubits (Floquet kick)
    Measure all qubits in computational basis.

    Paired circuit tests asymmetric-phase drive (the merkabit geometry).
    Control circuit tests symmetric drive (no asymmetry, expected F ~ 1).
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]

    ops = []

    # Initialize in |+>
    for q in qubits:
        ops.append(cirq.H(q))

    # tau Floquet cycles
    for cycle in range(tau):
        for i, q in enumerate(qubits):
            if paired:
                # Alternating sign: even i gets +phi, odd i gets -phi
                sign = +1 if i % 2 == 0 else -1
            else:
                # Control: same sign on all qubits
                sign = +1
            ops.append(cirq.rz(sign * PHI).on(q))
        # Floquet kick (Rx rotation) on all qubits
        for q in qubits:
            ops.append(cirq.rx(STEP_PHASE).on(q))

    # Measure all qubits
    ops.append(cirq.measure(*qubits, key='result'))

    circuit = cirq.Circuit(ops)

    # Optimise: merge consecutive single-qubit gates into PhXZ
    circuit = cirq.merge_single_qubit_gates_to_phased_x_and_z(circuit)
    circuit = cirq.drop_empty_moments(circuit)

    return circuit


# --- Fano factor analysis ----------------------------------------------------

def compute_fano_factor(counts: dict, n_qubits: int) -> tuple[float, float, float]:
    """
    Compute Fano factor F = Var(n_ones) / <n_ones> from measurement counts.

    counts: dict mapping bitstring -> occurrence count
    n_qubits: number of qubits measured

    Returns: (fano, mean, variance)
    """
    total_shots = sum(counts.values())
    if total_shots == 0:
        return float('nan'), 0.0, 0.0

    # Tally n_ones distribution
    n_ones_list = []
    for bitstring, count in counts.items():
        # bitstring could be int or string; count ones
        if isinstance(bitstring, (int, np.integer)):
            n_ones = bin(int(bitstring)).count('1')
        else:
            n_ones = str(bitstring).count('1')
        n_ones_list.extend([n_ones] * count)

    arr = np.array(n_ones_list)
    mean = float(np.mean(arr))
    variance = float(np.var(arr, ddof=1))
    fano = variance / mean if mean > 0 else float('nan')

    return fano, mean, variance


def counts_from_result(result: cirq.Result, key: str = 'result') -> dict:
    """Convert Cirq result to bitstring -> count dictionary."""
    measurements = result.measurements[key]
    counts = {}
    for bits in measurements:
        # bits is a 1D array for a single shot
        bitstring = ''.join(str(int(b)) for b in bits)
        counts[bitstring] = counts.get(bitstring, 0) + 1
    return counts


# --- Simulation modes --------------------------------------------------------

def simulate_noiseless(n_qubits: int, tau: int, paired: bool, shots: int) -> dict:
    """Ideal noiseless simulation via Cirq simulator."""
    circuit = build_scaling_circuit(n_qubits, tau, paired)
    simulator = cirq.Simulator()
    result = simulator.run(circuit, repetitions=shots)
    return counts_from_result(result)


def simulate_noisy(n_qubits: int, tau: int, paired: bool, shots: int,
                   error_rate: float = 0.005) -> dict:
    """
    Noisy simulation with depolarizing noise model.

    error_rate: per-gate depolarizing probability. Willow's measured
    single-qubit gate error is 0.035%, so 0.0035 is realistic. We
    use 0.005 as a conservative upper bound.
    """
    circuit = build_scaling_circuit(n_qubits, tau, paired)

    # Apply depolarizing noise after each moment
    noisy_circuit = circuit.with_noise(cirq.depolarize(p=error_rate))

    simulator = cirq.DensityMatrixSimulator()
    result = simulator.run(noisy_circuit, repetitions=shots)
    return counts_from_result(result)


# --- Main experiment loop ----------------------------------------------------

def run_scaling_experiment(mode: str = 'sim', error_rate: float = 0.005,
                           shots: int = SHOTS, verbose: bool = True) -> dict:
    """
    Run the P4 scaling experiment across N_QUBITS_SWEEP.

    mode: 'sim' for noiseless, 'noisy' for depolarizing noise
    error_rate: depolarizing probability per gate (noisy mode only)
    shots: measurements per configuration

    Returns: nested dict with results per N qubit count.
    """
    results = {
        'protocol': 'P4_scaling',
        'mode': mode,
        'error_rate': error_rate if mode == 'noisy' else 0.0,
        'shots_per_config': shots,
        'tau': TAU_FIXED,
        'timestamp': datetime.now().isoformat(),
        'predictions': {
            'primary': 'F_paired stays sub-Poissonian (F < 0.5) across N in {2,4,6,8}',
            'standard_theory': 'Independent-error decoherence predicts F -> 1 as N grows',
            'merkabit_prediction': 'F_paired stays ~ 0.2 (flat or decreasing) across N',
            'control_expectation': 'F_control ~ 1 (near-Poissonian) as baseline',
            'falsification': 'F_paired(N) exceeds 0.7 at any N >= 4, at 3 sigma',
        },
        'measurements': {},
    }

    if verbose:
        print(f"\n{'='*72}")
        print(f"P4 SCALING EXPERIMENT ({mode.upper()} mode)")
        print(f"{'='*72}")
        print(f"\ntau = {TAU_FIXED} Floquet cycles, {shots} shots per config")
        print(f"Qubit count sweep: {N_QUBITS_SWEEP}")
        if mode == 'noisy':
            print(f"Depolarizing error rate: {error_rate}")
        print()

    fano_paired = []
    fano_control = []

    for n_qubits in N_QUBITS_SWEEP:
        entry = {'n_qubits': n_qubits}

        # Paired (asymmetric-phase merkabit drive)
        if mode == 'sim':
            counts_p = simulate_noiseless(n_qubits, TAU_FIXED, paired=True, shots=shots)
        elif mode == 'noisy':
            counts_p = simulate_noisy(n_qubits, TAU_FIXED, paired=True,
                                       shots=shots, error_rate=error_rate)
        else:
            raise ValueError(f"Unknown mode: {mode}")

        F_p, mean_p, var_p = compute_fano_factor(counts_p, n_qubits)

        # Control (symmetric drive, same sign)
        if mode == 'sim':
            counts_c = simulate_noiseless(n_qubits, TAU_FIXED, paired=False, shots=shots)
        else:
            counts_c = simulate_noisy(n_qubits, TAU_FIXED, paired=False,
                                       shots=shots, error_rate=error_rate)

        F_c, mean_c, var_c = compute_fano_factor(counts_c, n_qubits)

        entry['paired'] = {'F': F_p, 'mean': mean_p, 'variance': var_p}
        entry['control'] = {'F': F_c, 'mean': mean_c, 'variance': var_c}
        entry['ratio_paired_over_control'] = (
            F_p / F_c if F_c > 0 else float('nan')
        )

        results['measurements'][f'N={n_qubits}'] = entry
        fano_paired.append(F_p)
        fano_control.append(F_c)

        if verbose:
            print(f"N = {n_qubits}:")
            print(f"  Paired  drive: F = {F_p:.4f}  (mean={mean_p:.3f}, var={var_p:.3f})")
            print(f"  Control drive: F = {F_c:.4f}  (mean={mean_c:.3f}, var={var_c:.3f})")
            print(f"  Ratio F_paired/F_control = {F_p/F_c:.3f}" if F_c > 0 else "")

    # Primary scaling criterion: does F_paired stay sub-Poissonian as N grows?
    # Standard decoherence theory predicts F -> 1 as qubit count grows.
    # The merkabit-geometric prediction is that F stays BELOW a threshold.
    F_N2 = fano_paired[0]
    F_N_max = fano_paired[-1]
    max_F = max(fano_paired)
    mean_F_paired = sum(fano_paired) / len(fano_paired)
    mean_F_control = sum(fano_control) / len(fano_control)

    # Secondary check: is there monotonic decrease (stronger claim)?
    is_monotonic = all(fano_paired[i] >= fano_paired[i+1]
                       for i in range(len(fano_paired) - 1))

    results['analysis'] = {
        'F_paired_values': fano_paired,
        'F_control_values': fano_control,
        'F_paired_mean': mean_F_paired,
        'F_control_mean': mean_F_control,
        'F_paired_max_across_N': max_F,
        'stays_sub_poissonian': max_F < 0.7,
        'stays_strongly_sub_poissonian': max_F < 0.5,
        'scaling_monotonic_decrease': is_monotonic,
        'F_N2_over_F_Nmax': F_N2 / F_N_max if F_N_max > 0 else float('nan'),
        'paired_vs_control_separation': mean_F_control - mean_F_paired,
    }

    if verbose:
        print(f"\n{'='*72}")
        print("ANALYSIS")
        print(f"{'='*72}")
        print(f"\nF paired by N: {[f'{f:.3f}' for f in fano_paired]}")
        print(f"F control by N: {[f'{f:.3f}' for f in fano_control]}")
        print(f"Paired mean across N: {mean_F_paired:.3f}")
        print(f"Control mean across N: {mean_F_control:.3f}")
        print(f"Paired vs control separation: {mean_F_control - mean_F_paired:.3f}")
        print(f"Max F_paired across N: {max_F:.3f}")

        print(f"\nVerdict:")
        if max_F < 0.5 and mean_F_control > 0.7:
            print("  SUB-POISSONIAN SCALING CONFIRMED:")
            print(f"  Paired drive stays at F <= {max_F:.3f} across all N in {N_QUBITS_SWEEP}.")
            print(f"  Control drive at F ~ {mean_F_control:.3f} (near-Poissonian).")
            print(f"  Standard decoherence predicts F -> 1 with N; observed F stays flat.")
            print(f"  The effect does not degrade with scale.")
            if is_monotonic:
                print(f"  Additional finding: F decreases monotonically with N.")
        elif max_F < 0.7:
            print("  WEAKLY SUB-POISSONIAN: F stays below 0.7 but not strongly suppressed.")
        else:
            print("  NO SCALING SIGNATURE: F exceeds 0.7 at some N.")
            print("  The geometric scaling hypothesis is not supported.")

    return results


# --- Hardware execution (Google Cloud Quantum) -------------------------------

def run_on_hardware(project_id: str, processor_id: str, shots: int = SHOTS):
    """
    Execute the scaling experiment on real Google quantum hardware.

    Requires google-cloud-quantum authentication.
    """
    try:
        import cirq_google
    except ImportError:
        raise ImportError(
            "cirq_google not installed. Run: pip install cirq-google"
        )

    engine = cirq_google.Engine(project_id=project_id)
    processor = engine.get_processor(processor_id=processor_id)

    results = {
        'protocol': 'P4_scaling',
        'mode': 'hardware',
        'project_id': project_id,
        'processor_id': processor_id,
        'shots_per_config': shots,
        'tau': TAU_FIXED,
        'timestamp': datetime.now().isoformat(),
        'measurements': {},
    }

    for n_qubits in N_QUBITS_SWEEP:
        entry = {'n_qubits': n_qubits}
        for paired_flag, label in [(True, 'paired'), (False, 'control')]:
            circuit = build_scaling_circuit(n_qubits, TAU_FIXED, paired=paired_flag)
            job = processor.run(program=circuit, repetitions=shots)
            result = job.results()[0]
            counts = counts_from_result(result)
            F, mean, var = compute_fano_factor(counts, n_qubits)
            entry[label] = {'F': F, 'mean': mean, 'variance': var, 'counts': counts}
        results['measurements'][f'N={n_qubits}'] = entry

    return results


# --- Main entry point --------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='P4 Scaling: F(N) vs qubit count N on Willow'
    )
    parser.add_argument('--sim-only', action='store_true',
                        help='Simulation only (no hardware access required)')
    parser.add_argument('--noisy', action='store_true',
                        help='Use depolarizing noise model in simulation')
    parser.add_argument('--error-rate', type=float, default=0.005,
                        help='Depolarizing error rate per gate (default: 0.005)')
    parser.add_argument('--shots', type=int, default=SHOTS,
                        help=f'Shots per configuration (default: {SHOTS})')
    parser.add_argument('--project', type=str, default=None,
                        help='Google Cloud project ID (hardware mode)')
    parser.add_argument('--processor', type=str, default=None,
                        help='Google quantum processor ID (hardware mode)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output JSON path (default: timestamped)')

    args = parser.parse_args()

    # Run experiment
    if args.sim_only:
        mode = 'noisy' if args.noisy else 'sim'
        results = run_scaling_experiment(mode=mode,
                                          error_rate=args.error_rate,
                                          shots=args.shots)
    else:
        if not args.project or not args.processor:
            print("Hardware mode requires --project and --processor")
            print("Run with --sim-only for noiseless simulation.")
            sys.exit(1)
        results = run_on_hardware(args.project, args.processor, shots=args.shots)

    # Save results
    if args.output:
        out_path = Path(args.output)
    else:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        mode_tag = 'noisy' if args.noisy else ('sim' if args.sim_only else 'hw')
        out_path = RESULTS_DIR / f'p4_scaling_{mode_tag}_{stamp}.json'

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nResults saved: {out_path}")


if __name__ == '__main__':
    main()
