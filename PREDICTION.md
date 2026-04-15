# Pre-registered Prediction: Merkabit on Square Grid

**Date:** 12 April 2026
**Authors:** Selina Stenberg, Thor Henning Hetland, Claude Anthropic

## Prediction

The merkabit protocol (P gate = Rz(+φ) ⊗ Rz(−φ), ouroboros 12-step Coxeter cycle)
will produce sub-Poissonian syndrome statistics (Fano factor F < 1) on Google
square-grid quantum hardware at syndrome extraction depth τ ≥ 5.

## Specific numerical predictions

| Observable | Predicted range | Basis |
|---|---|---|
| Fano factor at τ=5 (2-qubit, paired) | F ≈ 0.3–0.7 | Noisy topology simulation |
| Fano factor at τ=5 (2-qubit, control) | F ≈ 1 (near-Poissonian) | Standard decoherence theory |
| **P4: F(N) scaling, N ∈ {2, 4, 6, 8}** | **F stays sub-Poissonian (< 0.5) across all N** | **Cirq simulation; geometric hypothesis** |
| DTC ratio (paired/control) | > 3× | IBM cross-architecture result (3.2–3.9× on Eagle+Heron) |
| P1b Ramsey sign flip | n=6→8 | Angle table (topology-independent) |
| P2 quasi-period peak | n=39 ± 2 (= 3.25T) | Angle table (topology-independent) |
| P2 local min (n=13) | P ≈ 0.696 ± 0.05 | Theoretical value; IBM Eagle gave 0.688, Heron gave 0.709 |

## Falsification conditions

1. **F > 1 at all τ tested** (including τ ≥ 5) would indicate topology-dependence
2. **F_paired(N) > 0.7 at any N ≥ 4** would falsify the scaling prediction — standard decoherence (F → 1 with N) would then be the correct model
3. **DTC ratio < 1.5× (paired/control)** would indicate P gate has no effect on this topology
4. **No sign flip in P1b Ramsey** would indicate the Berry phase doesn't accumulate correctly
5. **P2 local min at n=13 outside [0.65, 0.75]** would indicate the angle table behaves differently on this platform

## Evidence base for this prediction

- 5/5 Appendix N predictions confirmed on IBM (6–12 April 2026)
- Cross-architecture validation: Eagle r3 and Heron r2 give matching results
- Topology simulation (sim_square_vs_hex_noisy.py): square grid produces F < 1 at τ ≥ 5
- Google implements Rz as virtual-Z frame change (same as IBM) — P gate is zero-cost

## This prediction is deposited before any Willow hardware run

The simulation results, prediction document, and experimental protocols are committed
to this repository with timestamps predating any Google hardware access. This constitutes
a pre-registered prediction testable by any group with access to a Google quantum processor.
