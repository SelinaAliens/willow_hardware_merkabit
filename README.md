# Willow Hardware Merkabit — Cross-Platform Validation

**Goal:** Test whether the merkabit's sub-Poissonian error statistics appear on Google's square-grid quantum processors, confirming the effect is platform-independent.

## Background

The merkabit protocol has been confirmed on IBM hardware:
- **Five of five** Appendix N predictions retired on IBM Eagle r3 and Heron r2
- Cross-architecture validation: Eagle r3 (ibm_strasbourg, ibm_brussels) and Heron r2 (ibm_kingston) show the same DTC, Berry phase, and quasi-period signatures
- The P gate (Rz(+φ) ⊗ Rz(−φ)) is native on IBM as two virtual-Z frame changes at zero cost

**The open question:** Does the same effect appear on Google Willow (square grid, different coupling map, different native gate set)?

Google also implements Rz as a virtual-Z frame change — the P gate should be zero-cost on Willow too.

## Simulation Results

The `simulations/` directory contains topology-comparison simulations run before any hardware access:

### Key finding: square grid IS sub-Poissonian at multi-round depth

| Comparison | τ | ε | F (hex) | F (square) | Both sub-P? |
|---|---|---|---|---|---|
| 7n hex vs 9n square | 5 | 0.01 | 0.065 | **0.041** | ✅ |
| 7n hex vs 9n square | 5 | 0.05 | 0.305 | **0.190** | ✅ |
| 7n hex vs 9n square | 5 | 0.10 | 0.522 | **0.343** | ✅ |
| 7n hex vs 9n square | 5 | 0.20 | 0.817 | **0.567** | ✅ |
| 7n hex vs 9n square | 5 | 0.30 | 0.949 | **0.701** | ✅ |

At τ ≥ 5, the square grid produces **tighter** anti-bunching than hex at every error rate tested. The ouroboros cycling drives the anti-bunching regardless of topology — the connectivity affects how quickly it develops (hex is better at τ=1), but both converge.

At τ = 1 (single round), the hex triangle is sub-Poissonian and the square is not — topology matters at short depth but not at multi-round depth.

### Prediction for Willow hardware

1. **τ = 1:** Fano likely super-Poissonian (consistent with existing Willow surface-code data F = 2.42)
2. **τ ≥ 5:** Fano should drop below 1, potentially F ≈ 0.3–0.7
3. **Per-round Fano stability:** should appear on Willow same as IBM (angle-table property, not topology)

## Experiments (pending hardware access)

The `experiments/` directory contains Cirq implementations of four protocols:
- **P1b Ramsey Berry phase** (2 qubits, zero 2-qubit gates)
- **P2 Stroboscopic quasi-period** (1 qubit, zero 2-qubit gates)
- **P4 Scaling: F(N) vs qubit count** (2, 4, 6, 8 qubits, zero 2-qubit gates)
- **P5 DTC survival with epsilon sweep** (2 qubits, zero 2-qubit gates)

All four protocols use only Rz and Rx (or H) — single-qubit natives on Google hardware.
After Cirq's single-qubit-merge optimisation, each Floquet cycle compiles to depth 2
(one PhXZ per qubit).

### P4 scaling — new prediction

Tests whether F_paired stays sub-Poissonian as qubit count scales. Cirq simulation
with Willow-realistic depolarizing noise (error rate 0.005) gives:

| N (qubits) | F_paired | F_control |
|---|---|---|
| 2 | 0.19 | 0.88 |
| 4 | 0.19 | 0.88 |
| 6 | 0.19 | 0.90 |
| 8 | 0.19 | 0.90 |

Standard decoherence theory predicts F → 1 as N grows (independent errors accumulate).
The merkabit-geometric prediction is that F_paired stays sub-Poissonian — this is the
direct test of whether the effect is a geometric pathway to coherent computing at scale.

Run it yourself:
```bash
python experiments/run_p4_scaling_cirq.py --sim-only --noisy --error-rate 0.005
```

## Analysis

The `analysis/` directory contains retrospective analysis of Google's
publicly available Willow surface-code dataset (Acharya et al., Nature 2025):

- `willow_fano_analysis.py` — Fano factor analysis of 420 Willow surface-code
  experiments, giving F = 2.42 (super-Poissonian).

This data does **not** falsify the pre-registered paired-phase predictions
in this repository. The Willow public circuits are symmetric-phase (unpaired)
stabiliser extraction — a different regime from the paired Floquet circuits
proposed here. See `analysis/README.md` for the full explanation.

## Hardware Requirements

- Any Google quantum processor with ≥ 2 qubits
- Native Rz (virtual-Z) and Rx (√X decomposition)
- 4,096 shots per circuit
- Total: ~200 circuits for the full P1/P2/P5 protocol
- Estimated QPU time: < 1 hour

## Related Repositories

- [merkabit_hardware_test](https://github.com/SelinaAliens/merkabit_hardware_test) — IBM hardware results (5/5 predictions confirmed)
- [The_Merkabit](https://github.com/SelinaAliens/The_Merkabit) — base simulations
- [merkabit-architecture](https://github.com/selinaserephina-star/merkabit-architecture) — architectural simulations

## Status

- [x] Topology simulation (hex vs square) — completed, prediction: sub-P on square at τ ≥ 5
- [ ] Cirq implementation of P1b, P2, P5
- [ ] Google Cloud / Willow access
- [ ] Hardware validation on square grid

## Authors

Selina Stenberg and Thor Henning Hetland with Claude Anthropic

April 2026
