# Analysis: Willow Public Dataset vs Pre-Registered Predictions

This directory contains retrospective analysis of Google's publicly available
Willow surface-code dataset (Acharya et al., *Nature* 638, 920–926, 2025;
Zenodo DOI: 10.5281/zenodo.13273331). The analysis establishes that the
existing Willow data does not falsify the pre-registered predictions for the
paired-phase Floquet experiments described in this repository.

## Context

A natural question for any reviewer evaluating the proposed Willow experiments
is: *does Google's existing public Willow data already refute the merkabit
hypothesis?*

The answer is no, but the reasoning requires understanding which circuit
classes were run on Willow and which are being proposed here.

## Files

- `willow_fano_analysis.py` — Full Fano factor analysis of 420 surface-code
  experiments from the Willow public dataset, along with burst scaling and
  correlation analyses.

## Key result

The Willow public data (420 surface-code experiments, distance 3 through 7)
shows Fano factor **F = 2.42** — super-Poissonian. Analysed in detail in
Paper 24 §7.

## Why this does NOT falsify the pre-registered prediction

The Willow public circuits are **symmetric-phase** (unpaired) stabiliser
extraction circuits. They apply the same Z rotation to every qubit at each
syndrome round. In the merkabit framework, this is the unpaired regime — the
P gate Rz(+φ) ⊗ Rz(−φ) is absent.

The pre-registered prediction in `PREDICTION.md` is for **paired** Floquet
circuits explicitly using asymmetric-phase drives. These circuits do not
exist in the Willow public dataset because no one has yet run paired-phase
Floquet experiments on Willow.

### Summary of the two regimes

| Regime | Willow data? | Measured F | Merkabit prediction |
|---|---|---|---|
| Symmetric-phase (unpaired) | Yes (Nature 2024, 420 experiments) | F = 2.42 | F ≈ 1 (or larger, if errors burst) |
| Asymmetric-phase (paired) | **Not yet run** | — | **F ∈ [0.3, 0.7] at τ ≥ 5** |

The paired regime is what this repository proposes to test on Willow.

## Paper 24 §7 correction

Paper 24 §7 originally concluded from this analysis that "sub-Poissonian
signal is architecture-specific to IBM's heavy-hex topology." Paper 26
§3 corrects this interpretation: the comparison was not hex-vs-square,
it was paired-vs-unpaired. The Willow super-Poissonian result is what
symmetric-phase circuits produce on any topology, which is consistent with
standard decoherence theory.

The merkabit-geometric prediction is that **paired** circuits on Willow
will show F < 1, not that **unpaired** circuits will. The public Willow
data tests the unpaired prediction (and confirms F > 1 as expected). The
paired prediction is what the proposed experiment tests.

## Running the analysis

The script reads Google's Willow public dataset (downloaded separately from
the Zenodo link in the reference). It computes:

1. Overall Fano factor across all 420 experiments
2. Per-experiment Fano with error bars
3. Burst scaling (how Fano changes with observation window)
4. Temporal correlation analysis

```bash
python willow_fano_analysis.py --data <path-to-willow-public-data>
```

## References

- Acharya et al., *Quantum error correction below the surface code threshold.*
  Nature 638, 920–926 (2025). [doi:10.1038/s41586-024-08449-y](https://doi.org/10.1038/s41586-024-08449-y)
- Willow public data: [Zenodo 10.5281/zenodo.13273331](https://zenodo.org/records/13273331)
- Paper 24: *The P Gate Is Native.* [Zenodo 10.5281/zenodo.19484743](https://zenodo.org/records/19484743)
- Paper 26: *The Merkabit Is Geometric.* [Zenodo 10.5281/zenodo.19554030](https://zenodo.org/records/19554030)
