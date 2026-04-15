"""
Google Willow 105Q Surface Code - Fano Factor Analysis
Compares syndrome statistics with IBM Eagle r3 results from Paper 15.
"""
import zipfile, json, sys, numpy as np
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')

ZIP_PATH = 'C:/Users/selin/Downloads/google_105Q_surface_code_d3_d5_d7.zip'
z = zipfile.ZipFile(ZIP_PATH)


def read_detection_counts(z, path, n_detectors, n_shots):
    raw = z.read(path)
    bytes_per_shot = (n_detectors + 7) // 8
    data = np.frombuffer(raw, dtype=np.uint8).reshape(n_shots, bytes_per_shot)
    counts = np.zeros(n_shots, dtype=np.int32)
    for byte_idx in range(bytes_per_shot):
        byte_col = data[:, byte_idx].astype(np.int32)
        for bit in range(min(8, n_detectors - byte_idx * 8)):
            counts += (byte_col >> bit) & 1
    return counts


def count_detectors(z, stim_path):
    text = z.read(stim_path).decode()
    return text.count('DETECTOR')


# Gather all experiments
metas = sorted([n for n in z.namelist() if n.endswith('metadata.json')])
results = []
by_distance = defaultdict(list)

print("=" * 100)
print("GOOGLE WILLOW 105Q SURFACE CODE - FANO FACTOR ANALYSIS")
print("Comparison with IBM Eagle r3 (Paper 15 / Paper 3)")
print("=" * 100)
print()

for i, meta_path in enumerate(metas):
    md = json.loads(z.read(meta_path))
    d = md['distance']
    rounds = md['rounds']
    shots = md['shots']
    basis = md['basis']
    parts = meta_path.split('/')
    patch = parts[1]

    stim_path = '/'.join(parts[:4]) + '/circuit_ideal.stim'
    n_det = count_detectors(z, stim_path)

    det_path = '/'.join(parts[:4]) + '/detection_events.b8'
    counts = read_detection_counts(z, det_path, n_det, shots)

    mean_c = np.mean(counts)
    var_c = np.var(counts, ddof=1)
    fano = var_c / mean_c if mean_c > 0 else np.nan

    rec = {
        'distance': d, 'patch': patch, 'basis': basis, 'rounds': rounds,
        'shots': shots, 'n_det': n_det, 'mean': mean_c, 'var': var_c,
        'fano': fano
    }
    results.append(rec)
    by_distance[d].append(rec)

    if (i + 1) % 50 == 0 or i == 0:
        print(f"  Processing {i+1}/{len(metas)}...", flush=True)

print(f"\nProcessed {len(results)} experiments total.\n")

# ── PART 1: by distance ──
print("=" * 100)
print("PART 1: FANO FACTOR BY CODE DISTANCE")
print("=" * 100)
print()
hdr = f"{'Dist':>5} | {'N':>5} | {'Mean F':>9} | {'Std F':>9} | {'Min':>8} | {'Max':>8} | {'Mean det':>9}"
print(hdr)
print("-" * len(hdr))

all_fanos = []
for d in sorted(by_distance.keys()):
    expts = by_distance[d]
    fanos = np.array([e['fano'] for e in expts])
    means = np.array([e['mean'] for e in expts])
    all_fanos.extend(fanos)
    print(f"{d:>5} | {len(expts):>5} | {np.mean(fanos):>9.4f} | {np.std(fanos):>9.4f} | "
          f"{np.min(fanos):>8.4f} | {np.max(fanos):>8.4f} | {np.mean(means):>9.2f}")

all_fanos = np.array(all_fanos)
se = np.std(all_fanos) / np.sqrt(len(all_fanos))
t_vs_1 = (np.mean(all_fanos) - 1.0) / se
t_vs_ibm = (np.mean(all_fanos) - 0.856) / se

print(f"\n  OVERALL: F = {np.mean(all_fanos):.4f} +/- {np.std(all_fanos):.4f}  (N = {len(all_fanos)})")
print(f"  t-test vs Poisson (F=1):  t = {t_vs_1:+.1f}")
print(f"  t-test vs IBM (F=0.856):  t = {t_vs_ibm:+.1f}")

if np.mean(all_fanos) > 1.0:
    print(f"\n  >>> SUPER-POISSONIAN: detection events are OVER-dispersed (F > 1)")
    print(f"  >>> OPPOSITE sign from IBM Eagle r3 (F = 0.856, sub-Poissonian)")
elif np.mean(all_fanos) < 1.0:
    print(f"\n  >>> SUB-POISSONIAN: detection events are UNDER-dispersed (F < 1)")
    print(f"  >>> SAME sign as IBM Eagle r3 (F = 0.856)")


# ── PART 2: by rounds ──
print()
print("=" * 100)
print("PART 2: FANO FACTOR BY QEC ROUNDS (does temporal depth drive the effect?)")
print("=" * 100)
print()

by_dr = defaultdict(list)
for r in results:
    by_dr[(r['distance'], r['rounds'])].append(r['fano'])

print(f"{'Dist':>5} | {'Rounds':>7} | {'Mean F':>9} | {'Std':>8} | {'N':>4}")
print("-" * 50)
for (d, rr) in sorted(by_dr.keys()):
    fs = by_dr[(d, rr)]
    print(f"{d:>5} | {rr:>7} | {np.mean(fs):>9.4f} | {np.std(fs):>8.4f} | {len(fs):>4}")


# ── PART 3: by patch ──
print()
print("=" * 100)
print("PART 3: FANO FACTOR BY PATCH (spatial location on chip)")
print("=" * 100)
print()

by_patch = defaultdict(list)
for r in results:
    by_patch[r['patch']].append(r['fano'])

print(f"{'Patch':>20} | {'Mean F':>9} | {'Std':>8} | {'N':>4}")
print("-" * 50)
for patch in sorted(by_patch.keys()):
    fs = by_patch[patch]
    print(f"{patch:>20} | {np.mean(fs):>9.4f} | {np.std(fs):>8.4f} | {len(fs):>4}")


# ── PART 4: burst scaling ──
print()
print("=" * 100)
print("PART 4: BURST SCALING (mean detection events vs code distance)")
print("=" * 100)
print()

for rounds in [1, 10, 13, 30, 50, 110, 250]:
    d_vals, m_vals = [], []
    for d in [3, 5, 7]:
        expts = [r for r in results if r['distance'] == d and r['rounds'] == rounds]
        if expts:
            d_vals.append(d)
            m_vals.append(np.mean([e['mean'] for e in expts]))

    if len(d_vals) == 3:
        d_arr = np.array(d_vals, dtype=float)
        m_arr = np.array(m_vals)

        # Linear fit: m = a*d + b
        A_lin = np.vstack([d_arr, np.ones(3)]).T
        coef_lin = np.linalg.lstsq(A_lin, m_arr, rcond=None)[0]
        pred_lin = A_lin @ coef_lin
        ss_res_lin = np.sum((m_arr - pred_lin)**2)
        ss_tot = np.sum((m_arr - np.mean(m_arr))**2)
        r2_lin = 1 - ss_res_lin / ss_tot if ss_tot > 0 else 0

        # Quadratic fit: m = a*d^2 + b
        A_quad = np.vstack([d_arr**2, np.ones(3)]).T
        coef_quad = np.linalg.lstsq(A_quad, m_arr, rcond=None)[0]
        pred_quad = A_quad @ coef_quad
        ss_res_quad = np.sum((m_arr - pred_quad)**2)
        r2_quad = 1 - ss_res_quad / ss_tot if ss_tot > 0 else 0

        r53 = m_vals[1] / m_vals[0]
        r75 = m_vals[2] / m_vals[1]
        print(f"  r={rounds:>3}: d3={m_vals[0]:>7.1f}  d5={m_vals[1]:>7.1f}  d7={m_vals[2]:>7.1f}")
        print(f"         R2_lin={r2_lin:.6f}  R2_quad={r2_quad:.6f}")
        print(f"         d5/d3={r53:.3f} (lin:1.667, quad:2.778)  d7/d5={r75:.3f} (lin:1.400, quad:1.960)")
        print()


# ── PART 5: Fano normalized by rounds ──
print("=" * 100)
print("PART 5: FANO FACTOR AT FIXED SHORT ROUNDS (r=1, r=10) vs LONG (r=250)")
print("  If super-Poisson comes from temporal correlations,")
print("  Fano should grow with rounds.")
print("=" * 100)
print()

for d in [3, 5, 7]:
    for rr in [1, 10, 250]:
        expts = [r for r in results if r['distance'] == d and r['rounds'] == rr]
        if expts:
            fanos = [e['fano'] for e in expts]
            print(f"  d={d}, r={rr:>3}: F = {np.mean(fanos):.4f} +/- {np.std(fanos):.4f}  (N={len(expts)})")
    print()


# ── PART 6: per-round Fano ──
print("=" * 100)
print("PART 6: DETECTION RATE PER ROUND PER DETECTOR")
print("  (normalizing out the effect of round count)")
print("=" * 100)
print()

for d in [3, 5, 7]:
    for rr in [1, 13, 50, 250]:
        expts = [r for r in results if r['distance'] == d and r['rounds'] == rr]
        if expts:
            rates = [e['mean'] / e['n_det'] for e in expts]
            print(f"  d={d}, r={rr:>3}: det_per_detector = {np.mean(rates):.4f} "
                  f"(per_round = {np.mean(rates)/rr:.6f})")
    print()


# ── SUMMARY ──
print("=" * 100)
print("SUMMARY COMPARISON")
print("=" * 100)
print()
print("  IBM Eagle r3 (Paper 15, 756 QEC runs):")
print("    Fano = 0.856 +/- 0.030   (SUB-Poissonian)")
print("    d=3: 0.855   d=5: 0.857   d=7: 0.855")
print("    ANOVA across distances: p = 0.79 (no trend)")
print("    Burst scaling: LINEAR (R2 = 0.9999)")
print()
print("  Google Willow (this analysis, 420 experiments):")
print(f"    Fano = {np.mean(all_fanos):.4f} +/- {np.std(all_fanos):.4f}")
for d in sorted(by_distance.keys()):
    fanos = [e['fano'] for e in by_distance[d]]
    print(f"    d={d}: {np.mean(fanos):.4f}")

# ANOVA across distances
from scipy import stats as sp_stats
groups = [np.array([e['fano'] for e in by_distance[d]]) for d in sorted(by_distance.keys())]
if len(groups) == 3:
    F_stat, p_anova = sp_stats.f_oneway(*groups)
    print(f"    ANOVA across distances: F={F_stat:.2f}, p={p_anova:.4f}")

print()
