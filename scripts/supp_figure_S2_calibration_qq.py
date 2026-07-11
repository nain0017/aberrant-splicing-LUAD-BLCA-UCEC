"""
step3c_S2_calibration_qq_v1.py
Paper 5A - Supplementary Figure S2.

Beta-binomial LRT calibration QQ-plot via label permutation.

Procedure (per cancer):
  1. Load PSI matrix + integer counts for coverage-pass events.
  2. Sample N_EVENTS random events.
  3. For each, permute the tumor/normal labels once, refit the BB LRT
     under shared-mu null vs free-mu alt, get a p-value from chi-sq df=1.
  4. Plot expected vs observed -log10(p) as a QQ-plot with 95% envelope.
  5. Report lambda (genomic inflation factor) per cancer.

Three panels: LUAD, BLCA, UCEC.
Saves: Stage2C_corrected/figures/S2_calibration_qq_v1.pdf + .png

Reads:
  data/processed/psi/{LUAD,BLCA,UCEC}_psi.npz
  data/processed/psi/{LUNG,BLADDER,UTERUS}_psi.npz
  data/processed/psi/suppa_events_parsed.parquet
  data/processed/psi/Stage2C_corrected/{LUAD,BLCA,UCEC}_corrected_results.parquet

USAGE IN COLAB:
  After mounting Drive, paste this file into one cell and run.
  Runtime: ~8-12 minutes (3 cancers x ~2000 events x 1 perm).
  Tweak N_EVENTS_PER_CANCER if you want faster / slower.
"""

import os
import sys
import time
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.optimize import minimize
from scipy.special import betaln
from scipy import stats


# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
PROJECT_ROOT = "/content/drive/MyDrive/Paper5_SharedAntigens"
PSI_DIR      = os.path.join(PROJECT_ROOT, "data/processed/psi")
STAGE_C_DIR  = os.path.join(PSI_DIR, "Stage2C_corrected")
FIG_DIR      = os.path.join(STAGE_C_DIR, "figures")

CANCERS = ["LUAD", "BLCA", "UCEC"]
NORMALS = {"LUAD": "LUNG", "BLCA": "BLADDER", "UCEC": "UTERUS"}

N_EVENTS_PER_CANCER = 2000   # subsampling for calibration
COVERAGE_FLOOR      = 10     # per-sample reads to be 'well-covered'
MIN_WELL_COVERED    = 20     # per side
RNG_SEED            = 42


# ----------------------------------------------------------------------
# Style
# ----------------------------------------------------------------------
STYLE = {
    'figure.dpi':       150,
    'savefig.dpi':      300,
    'font.family':      'sans-serif',
    'font.sans-serif':  ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size':        8.0,
    'axes.titlesize':   9.0,
    'axes.titleweight': 'bold',
    'axes.labelsize':   8.0,
    'xtick.labelsize':  7.0,
    'ytick.labelsize':  7.0,
    'axes.linewidth':   0.7,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'pdf.fonttype':     42,
    'ps.fonttype':      42,
}
COHORT_COLOR = {"LUAD": "#2E86AB", "BLCA": "#E07A5F", "UCEC": "#76A646"}


def apply_style():
    for k, v in STYLE.items():
        mpl.rcParams[k] = v


# ----------------------------------------------------------------------
# Beta-binomial LRT (identical to cell 107 implementation)
# ----------------------------------------------------------------------
def _bb_neg_loglik(params, k_t, n_t, k_n, n_n, mode):
    if mode == 'alt':
        mu_t, mu_n, log_phi = params
    else:
        mu_s, log_phi = params
        mu_t = mu_n = mu_s
    if not (1e-6 < mu_t < 1 - 1e-6 and 1e-6 < mu_n < 1 - 1e-6):
        return 1e10
    phi = np.exp(log_phi)
    a_t, b_t = mu_t * phi, (1 - mu_t) * phi
    a_n, b_n = mu_n * phi, (1 - mu_n) * phi
    ll_t = (betaln(k_t + a_t, n_t - k_t + b_t) - betaln(a_t, b_t)).sum()
    ll_n = (betaln(k_n + a_n, n_n - k_n + b_n) - betaln(a_n, b_n)).sum()
    return -(ll_t + ll_n)


def betabin_lrt_one_event(k_t, n_t, k_n, n_n):
    """Return p-value from chi-sq df=1 LRT, or NaN if fit fails."""
    # Sensible initial guesses
    mu_t0 = float(np.clip(k_t.sum() / max(n_t.sum(), 1), 0.05, 0.95))
    mu_n0 = float(np.clip(k_n.sum() / max(n_n.sum(), 1), 0.05, 0.95))
    mu_s0 = float(np.clip((k_t.sum() + k_n.sum()) /
                          max(n_t.sum() + n_n.sum(), 1), 0.05, 0.95))
    log_phi0 = 3.0
    try:
        alt = minimize(_bb_neg_loglik,
                       x0=[mu_t0, mu_n0, log_phi0],
                       args=(k_t, n_t, k_n, n_n, 'alt'),
                       method='Nelder-Mead',
                       options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 200})
        nul = minimize(_bb_neg_loglik,
                       x0=[mu_s0, log_phi0],
                       args=(k_t, n_t, k_n, n_n, 'null'),
                       method='Nelder-Mead',
                       options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 200})
        if not (alt.success and nul.success):
            return np.nan, np.nan
        chi2 = 2 * (nul.fun - alt.fun)
        if chi2 < 0:
            chi2 = 0.0
        p = stats.chi2.sf(chi2, df=1)
        return p, chi2
    except Exception:
        return np.nan, np.nan


# ----------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------
def load_npz(label):
    path = os.path.join(PSI_DIR, f"{label}_psi.npz")
    if not os.path.exists(path):
        sys.exit(f"Missing {path}")
    data = np.load(path, allow_pickle=True)
    keys = list(data.keys())
    print(f"  [{label}] {os.path.basename(path)}  keys={keys}")
    return data


def load_coverage_pass_events(c):
    path = os.path.join(STAGE_C_DIR, f"{c}_corrected_results.parquet")
    df = pd.read_parquet(path)
    cov = df['coverage_pass'].fillna(False).astype(bool)
    return df[cov]['event_id'].astype(str).values


def reconstruct_counts(psi, coverage):
    """rint reconstruction matching cell 107."""
    total_int = np.rint(coverage).astype(np.int32)
    psi_safe  = np.where(np.isnan(psi), 0.0, psi).astype(np.float64)
    inc_int   = np.rint(psi_safe * total_int.astype(np.float64)).astype(np.int32)
    inc_int   = np.clip(inc_int, 0, total_int)
    return inc_int, total_int


# ----------------------------------------------------------------------
# Per-cancer calibration run
# ----------------------------------------------------------------------
def run_calibration(c, suppa_eid_to_idx, n_events, rng):
    print(f"\n--- Calibration for {c} ---")
    t_npz = load_npz(c)
    n_npz = load_npz(NORMALS[c])

    # Required arrays
    psi_t = t_npz['psi']         # (n_events, n_t_samples)
    cov_t = t_npz['coverage']    # same shape
    psi_n = n_npz['psi']
    cov_n = n_npz['coverage']

    print(f"  psi_t shape: {psi_t.shape}    psi_n shape: {psi_n.shape}")

    # Coverage-pass event ids and their row positions in the PSI matrix
    cov_pass_eids = load_coverage_pass_events(c)
    cov_pass_idx = [suppa_eid_to_idx.get(e) for e in cov_pass_eids]
    cov_pass_idx = [i for i in cov_pass_idx if i is not None]
    print(f"  coverage-pass events tracked: {len(cov_pass_idx):,}")

    # Subsample
    n_take = min(n_events, len(cov_pass_idx))
    sample_idx = rng.choice(cov_pass_idx, size=n_take, replace=False)
    print(f"  sampling {n_take:,} events for calibration")

    # Reconstruct integer counts on the subsample only
    k_t_full, n_t_full = reconstruct_counts(psi_t[sample_idx], cov_t[sample_idx])
    k_n_full, n_n_full = reconstruct_counts(psi_n[sample_idx], cov_n[sample_idx])

    n_t_samples = psi_t.shape[1]
    n_n_samples = psi_n.shape[1]
    total_samples = n_t_samples + n_n_samples

    pvals = np.full(n_take, np.nan)
    chi2s = np.full(n_take, np.nan)

    t0 = time.time()
    n_fit_success = 0
    for i in range(n_take):
        kt_row = k_t_full[i]
        nt_row = n_t_full[i]
        kn_row = k_n_full[i]
        nn_row = n_n_full[i]

        # Well-covered masks (>= COVERAGE_FLOOR)
        ok_t = nt_row >= COVERAGE_FLOOR
        ok_n = nn_row >= COVERAGE_FLOOR

        if ok_t.sum() < MIN_WELL_COVERED or ok_n.sum() < MIN_WELL_COVERED:
            continue

        kt_v = kt_row[ok_t]
        nt_v = nt_row[ok_t]
        kn_v = kn_row[ok_n]
        nn_v = nn_row[ok_n]

        # ----- LABEL PERMUTATION -----
        # Pool samples then re-split at random into "tumor-sized" and
        # "normal-sized" groups. This preserves per-sample (k, n) pairs
        # so dispersion structure stays realistic.
        k_pool = np.concatenate([kt_v, kn_v])
        n_pool = np.concatenate([nt_v, nn_v])
        n_pool_total = len(k_pool)
        if n_pool_total < 2 * MIN_WELL_COVERED:
            continue
        perm = rng.permutation(n_pool_total)
        k_pool = k_pool[perm]
        n_pool = n_pool[perm]

        # Split into pseudo-T and pseudo-N of original sizes
        nT = len(kt_v)
        kt_perm = k_pool[:nT]
        nt_perm = n_pool[:nT]
        kn_perm = k_pool[nT:]
        nn_perm = n_pool[nT:]

        p, chi2 = betabin_lrt_one_event(kt_perm, nt_perm, kn_perm, nn_perm)
        if not np.isnan(p):
            pvals[i] = p
            chi2s[i] = chi2
            n_fit_success += 1

        if (i + 1) % 250 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / max(elapsed, 0.01)
            eta = (n_take - i - 1) / max(rate, 0.01)
            print(f"    {i+1:,}/{n_take:,}  ({rate:.1f}/s, ETA {eta:.0f}s)  "
                  f"successes={n_fit_success}")

    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s. Fit successes: {n_fit_success}/{n_take}")

    valid = ~np.isnan(pvals)
    print(f"  valid p-values: {valid.sum()}/{n_take}")
    return pvals[valid], chi2s[valid]


# ----------------------------------------------------------------------
# QQ-plot + lambda
# ----------------------------------------------------------------------
def compute_lambda(chi2s):
    """Genomic inflation lambda = median(chi2) / median expected under chi2_1."""
    if len(chi2s) == 0:
        return np.nan
    obs_median = np.median(chi2s)
    expected_median = stats.chi2.ppf(0.5, df=1)
    return obs_median / expected_median


def plot_qq(ax, c, pvals, lam):
    color = COHORT_COLOR[c]
    n = len(pvals)
    if n == 0:
        ax.text(0.5, 0.5, "No valid p-values",
                transform=ax.transAxes, ha='center', va='center')
        ax.set_axis_off()
        return

    sorted_p = np.sort(pvals)
    # Expected uniform quantiles
    expected_unif = (np.arange(1, n + 1) - 0.5) / n
    obs_log = -np.log10(np.clip(sorted_p, 1e-30, 1.0))
    exp_log = -np.log10(expected_unif)

    # 95% pointwise envelope (Beta(i, n-i+1) for i-th order stat under unif)
    i = np.arange(1, n + 1)
    lo_unif = stats.beta.ppf(0.025, i, n - i + 1)
    hi_unif = stats.beta.ppf(0.975, i, n - i + 1)
    lo_log = -np.log10(np.clip(hi_unif, 1e-30, 1.0))   # NB: invert
    hi_log = -np.log10(np.clip(lo_unif, 1e-30, 1.0))

    ax.fill_between(exp_log, lo_log, hi_log,
                    color='#CCCCCC', alpha=0.5, linewidth=0,
                    label='95% envelope')

    ax.scatter(exp_log, obs_log, s=5, c=color, alpha=0.45,
               edgecolors='none', rasterized=True)

    # Diagonal
    mx = max(exp_log.max(), obs_log.max()) * 1.05
    ax.plot([0, mx], [0, mx], color='black', linewidth=0.7, linestyle='--')

    # Aesthetics
    ax.set_xlim(0, mx)
    ax.set_ylim(0, mx)
    ax.set_xlabel(r'Expected $-\log_{10}(p)$ under uniform null')
    ax.set_ylabel(r'Observed $-\log_{10}(p)$')
    ax.set_title(f"{c}  (n={n:,} permuted events)",
                 loc='center', color=color)

    # Lambda inset
    verdict = "clean"
    if lam > 1.15:
        verdict = "elevated"
    elif lam > 1.05:
        verdict = "borderline"
    ax.text(0.04, 0.96,
            f"$\\lambda$ = {lam:.3f}\n({verdict})",
            transform=ax.transAxes, ha='left', va='top',
            fontsize=7.5,
            bbox=dict(boxstyle='round,pad=0.3',
                      facecolor='white', edgecolor='black',
                      linewidth=0.4))


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    if not os.path.exists("/content/drive/MyDrive"):
        sys.exit("Drive not mounted. Run drive.mount('/content/drive') first.")
    if not os.path.exists(STAGE_C_DIR):
        sys.exit(f"Stage2C_corrected not found: {STAGE_C_DIR}")
    os.makedirs(FIG_DIR, exist_ok=True)

    apply_style()
    rng = np.random.default_rng(RNG_SEED)
    print(f"Output: {FIG_DIR}")
    print(f"Date:   {date.today().isoformat()}")
    print(f"N_EVENTS_PER_CANCER = {N_EVENTS_PER_CANCER}")
    print(f"COVERAGE_FLOOR per sample = {COVERAGE_FLOOR}")
    print(f"MIN_WELL_COVERED per side = {MIN_WELL_COVERED}")
    print()

    # SUPPA events -> event_id to matrix row index
    suppa_path = os.path.join(PSI_DIR, "suppa_events_parsed.parquet")
    if not os.path.exists(suppa_path):
        sys.exit(f"Missing {suppa_path}")
    suppa = pd.read_parquet(suppa_path)
    id_col = None
    for c in ['event_id', 'event', 'event_name', 'id']:
        if c in suppa.columns:
            id_col = c; break
    if id_col is None:
        sys.exit(f"No event id column in SUPPA parquet. Got: {list(suppa.columns)}")
    suppa_eid_to_idx = {str(e): i for i, e in enumerate(suppa[id_col])}
    print(f"SUPPA events: {len(suppa):,}; id column='{id_col}'")
    print()

    results = {}
    for c in CANCERS:
        pvals, chi2s = run_calibration(c, suppa_eid_to_idx,
                                       N_EVENTS_PER_CANCER, rng)
        lam = compute_lambda(chi2s)
        results[c] = (pvals, chi2s, lam)
        print(f"  [{c}] lambda = {lam:.3f}")

    # Save raw outputs for tables / supplementary
    out_csv = os.path.join(STAGE_C_DIR, "S2_calibration_lambdas.csv")
    pd.DataFrame({
        'cancer': list(results.keys()),
        'n_events_perm': [len(results[c][0]) for c in results],
        'lambda':        [results[c][2]      for c in results],
    }).to_csv(out_csv, index=False)
    print(f"\n  Lambda summary saved: {out_csv}")

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.4))
    for ax, c in zip(axes, CANCERS):
        pvals, chi2s, lam = results[c]
        plot_qq(ax, c, pvals, lam)

    plt.tight_layout(w_pad=2.0)

    pdf_path = os.path.join(FIG_DIR, "S2_calibration_qq_v1.pdf")
    png_path = os.path.join(FIG_DIR, "S2_calibration_qq_v1.png")
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    print(f"\n  Saved PDF: {pdf_path}")
    print(f"  Saved PNG: {png_path}")
    print()
    print("S2 v1 done. Open the PNG, review, then send the image back.")
    print("Lambdas to interpret in Methods:")
    for c in CANCERS:
        print(f"  {c}: lambda = {results[c][2]:.3f}")


if __name__ == "__main__":
    main()