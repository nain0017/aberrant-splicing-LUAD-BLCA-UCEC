"""
step3c_F1_pipeline_overview_v2.py
Paper 5A — Figure 1, corrected-values-only version.

Replaces the cell 124 F1 by dropping Panel A's old-vs-new comparison. Four
panels, all sourced directly from the corrected parquets, no inflated counts
anywhere in the figure.

  Panel A: Stacked bars of primary + secondary aberrant counts per cancer,
           with per-cancer M2' threshold lines (LUAD>=500, BLCA>=300, UCEC>=300).
  Panel B: BB primary vs Wilcoxon sensitivity test concordance (both /
           BB-only / MW-only / neither) on coverage-pass events.
  Panel C: Stage funnel per cancer: total -> coverage-pass -> q<0.05 ->
           |dPSI|>=0.15 -> primary (3C-validated).
  Panel D: Boxplot of well-covered tumor samples per event in the primary
           aberrant set per cancer.

Saves: Stage2C_corrected/figures/F1_pipeline_overview_v2.pdf + .png

Reads only:
  Stage2C_corrected/{LUAD,BLCA,UCEC}_corrected_results.parquet

USAGE IN COLAB:
  Paste this file's contents into one cell and run.
  Runtime: < 20 seconds.
"""

import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
PROJECT_ROOT = "/content/drive/MyDrive/Paper5_SharedAntigens"
STAGE_C_DIR  = os.path.join(PROJECT_ROOT, "data/processed/psi/Stage2C_corrected")
FIG_DIR      = os.path.join(STAGE_C_DIR, "figures")

CANCERS = ["LUAD", "BLCA", "UCEC"]

# Per-cancer M2' thresholds (from STAGE_C_FIX_PRECOMMIT.md Section 3)
M2_PRIME = {"LUAD": 500, "BLCA": 300, "UCEC": 300}


# ----------------------------------------------------------------------
# Style sheet (locked, identical to cell 124)
# ----------------------------------------------------------------------
STYLE = {
    'figure.figsize':        (7.2, 5.4),
    'figure.dpi':            150,
    'savefig.dpi':           300,
    'font.family':           'sans-serif',
    'font.sans-serif':       ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size':             8.0,
    'axes.titlesize':        9.0,
    'axes.titleweight':      'bold',
    'axes.labelsize':        8.0,
    'xtick.labelsize':       7.0,
    'ytick.labelsize':       7.0,
    'legend.fontsize':       7.0,
    'axes.linewidth':        0.7,
    'xtick.major.width':     0.6,
    'ytick.major.width':     0.6,
    'axes.spines.top':       False,
    'axes.spines.right':     False,
    'axes.grid':             False,
    'pdf.fonttype':          42,
    'ps.fonttype':           42,
}

COHORT_COLOR = {
    "LUAD":  "#2E86AB",
    "BLCA":  "#E07A5F",
    "UCEC":  "#76A646",
    "UNION": "#3D405B",
}

TIER_COLOR = {
    "primary":   "#2E5A88",
    "secondary": "#B9D1E3",
}


def apply_style():
    for k, v in STYLE.items():
        mpl.rcParams[k] = v


# ----------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------
def load_corrected_results():
    results = {}
    for c in CANCERS:
        path = os.path.join(STAGE_C_DIR, f"{c}_corrected_results.parquet")
        if not os.path.exists(path):
            sys.exit(f"Missing parquet: {path}")
        df = pd.read_parquet(path)
        results[c] = df
        n_pri = int((df['tier'] == 'primary').sum())
        n_sec = int((df['tier'] == 'secondary').sum())
        print(f"  [{c}] loaded {len(df):,} events; primary={n_pri}, secondary={n_sec}")
    return results


# ----------------------------------------------------------------------
# Panel A — primary + secondary stacked per cancer, with M2' threshold lines
# ----------------------------------------------------------------------
def panel_A_tier_counts(ax, results):
    cancers = CANCERS
    n_pri = [int((results[c]['tier'] == 'primary').sum())   for c in cancers]
    n_sec = [int((results[c]['tier'] == 'secondary').sum()) for c in cancers]

    x = np.arange(len(cancers))
    w = 0.55

    # Stacked: primary (saturated cohort color) + secondary (lighter)
    ax.bar(x, n_pri, w,
           color=[COHORT_COLOR[c] for c in cancers],
           edgecolor='black', linewidth=0.5,
           label='Primary (3C-validated)')
    ax.bar(x, n_sec, w, bottom=n_pri,
           color=[COHORT_COLOR[c] for c in cancers], alpha=0.35,
           edgecolor='black', linewidth=0.5,
           label='Secondary (insufficient SCN)')

    # M2' threshold ticks on each bar
    for xi, c in zip(x, cancers):
        thr = M2_PRIME[c]
        ax.hlines(thr, xi - w/2, xi + w/2,
                  colors='black', linestyles='--', linewidth=0.8)
        ax.text(xi + w/2 + 0.05, thr, f"M2' = {thr}",
                fontsize=6.5, va='center', ha='left', color='black')

    # Per-bar count labels (primary on its segment, secondary on top of stack)
    for xi, p, s in zip(x, n_pri, n_sec):
        ax.text(xi, p/2, f"{p:,}", ha='center', va='center',
                fontsize=7.5, color='white', fontweight='bold')
        if s > 0:
            ax.text(xi, p + s + max(n_pri+n_sec)*0.015,
                    f"+{s:,}", ha='center', va='bottom',
                    fontsize=6.5, color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(cancers)
    ax.set_ylabel('Aberrant events')
    ax.set_title('a   Aberrant event counts by tier', loc='left')
    ax.set_ylim(0, max(np.array(n_pri) + np.array(n_sec)) * 1.18)
    ax.legend(loc='upper left', frameon=False, fontsize=6.5,
              handlelength=1.2, handletextpad=0.5)


# ----------------------------------------------------------------------
# Panel B — BB vs MW concordance on coverage-pass events
# ----------------------------------------------------------------------
def panel_B_concordance(ax, results):
    cancers = CANCERS
    both, bb_only, mw_only, neither = [], [], [], []
    for c in cancers:
        df = results[c]
        cov = df['coverage_pass'].fillna(False).astype(bool)
        bb  = df['passes_q'].fillna(False).astype(bool)
        mw  = (df['mwu_pvalue'] < 0.05) & df['mwu_pvalue'].notna()
        both.append(   int((cov &  bb &  mw).sum()))
        bb_only.append(int((cov &  bb & ~mw).sum()))
        mw_only.append(int((cov & ~bb &  mw).sum()))
        neither.append(int((cov & ~bb & ~mw).sum()))

    both    = np.array(both)
    bb_only = np.array(bb_only)
    mw_only = np.array(mw_only)
    neither = np.array(neither)
    total   = both + bb_only + mw_only + neither

    pct_both    = 100 * both    / total
    pct_bb_only = 100 * bb_only / total
    pct_mw_only = 100 * mw_only / total
    pct_neither = 100 * neither / total

    x = np.arange(len(cancers))
    w = 0.55

    ax.bar(x, pct_both,    w, color='#3D5A80', label='Both significant',
           edgecolor='black', linewidth=0.4)
    ax.bar(x, pct_bb_only, w, bottom=pct_both,
           color='#98C1D9', label='BB only',
           edgecolor='black', linewidth=0.4)
    ax.bar(x, pct_mw_only, w, bottom=pct_both + pct_bb_only,
           color='#EE6C4D', label='MW only',
           edgecolor='black', linewidth=0.4)
    ax.bar(x, pct_neither, w, bottom=pct_both + pct_bb_only + pct_mw_only,
           color='#E0E0E0', label='Neither',
           edgecolor='black', linewidth=0.4)

    for xi, pb, n in zip(x, pct_both, total):
        ax.text(xi, pb/2, f"{pb:.0f}%", ha='center', va='center',
                fontsize=7.5, color='white', fontweight='bold')
        ax.text(xi, 102, f"n={n:,}", ha='center', va='bottom',
                fontsize=6.5, color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(cancers)
    ax.set_ylabel('% of coverage-pass events')
    ax.set_title('b   Beta-binomial vs Wilcoxon concordance', loc='left')
    ax.set_ylim(0, 112)
    ax.legend(loc='lower right', frameon=False, fontsize=6.5,
              handlelength=1.0, handletextpad=0.5, ncol=2,
              columnspacing=0.8)


# ----------------------------------------------------------------------
# Panel C — Stage funnel per cancer
# ----------------------------------------------------------------------
def panel_C_funnel(ax, results):
    cancers = CANCERS
    stages = ['Total\nevents', 'Coverage\npass', 'q < 0.05\n(BH)',
              '|dPSI|\n>= 0.15', 'Primary\n(3C-val.)']

    counts = {c: [] for c in cancers}
    for c in cancers:
        df = results[c]
        cov = df['coverage_pass'].fillna(False).astype(bool)
        q   = df['passes_q'].fillna(False).astype(bool)
        dp  = df['passes_dpsi'].fillna(False).astype(bool)
        counts[c] = [
            len(df),
            int(cov.sum()),
            int((cov & q).sum()),
            int((cov & q & dp).sum()),
            int((df['tier'] == 'primary').sum()),
        ]

    x = np.arange(len(stages))
    w = 0.27
    offsets = {"LUAD": -w, "BLCA": 0.0, "UCEC": +w}

    for c in cancers:
        ax.bar(x + offsets[c], counts[c], w,
               color=COHORT_COLOR[c], edgecolor='black', linewidth=0.4,
               label=c)
        for xi, v in zip(x, counts[c]):
            if v >= 1000:
                lbl = f"{v/1000:.1f}K" if v < 10000 else f"{v//1000}K"
            else:
                lbl = f"{v}"
            ax.text(xi + offsets[c], v + max(counts[c]) * 0.018,
                    lbl, ha='center', va='bottom',
                    fontsize=5.8, rotation=0, color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    ax.set_yscale('log')
    ax.set_ylim(100, 80000)
    ax.set_ylabel('Events (log scale)')
    ax.set_title('c   Filtering funnel by cancer', loc='left')
    ax.legend(loc='upper right', frameon=False, fontsize=6.5,
              handlelength=1.0, handletextpad=0.5, ncol=3,
              columnspacing=0.8)


# ----------------------------------------------------------------------
# Panel D — Coverage of primary set (tumor well-covered samples)
# ----------------------------------------------------------------------
def panel_D_coverage(ax, results):
    cancers = CANCERS
    data, labels, colors = [], [], []
    for c in cancers:
        df = results[c]
        pri = df[df['tier'] == 'primary']
        data.append(pri['n_tumor_well_covered'].values)
        labels.append(c)
        colors.append(COHORT_COLOR[c])

    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, widths=0.55,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=2,
                                    markerfacecolor='black',
                                    markeredgecolor='none', alpha=0.4),
                    medianprops=dict(color='black', linewidth=1.1),
                    whiskerprops=dict(color='black', linewidth=0.7),
                    capprops=dict(color='black', linewidth=0.7),
                    boxprops=dict(edgecolor='black', linewidth=0.7))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)

    # Coverage gate line (>= 20 well-covered tumor samples)
    ax.axhline(20, color='red', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.text(len(cancers) + 0.45, 20, 'Gate >= 20',
            fontsize=6.5, va='center', ha='left', color='red')

    ax.set_ylabel('Well-covered tumor samples\nper primary event')
    ax.set_title('d   Coverage of primary aberrant set', loc='left')
    ax.set_ylim(0, max(max(d) for d in data) * 1.05)


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
    print(f"Output: {FIG_DIR}")
    print(f"Date:   {date.today().isoformat()}")
    print()
    print("Loading corrected results parquets ...")
    results = load_corrected_results()

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    panel_A_tier_counts(axes[0, 0], results)
    panel_B_concordance(axes[0, 1], results)
    panel_C_funnel(axes[1, 0], results)
    panel_D_coverage(axes[1, 1], results)

    plt.tight_layout(w_pad=2.2, h_pad=2.4)

    pdf_path = os.path.join(FIG_DIR, "F1_pipeline_overview_v2.pdf")
    png_path = os.path.join(FIG_DIR, "F1_pipeline_overview_v2.png")
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    print(f"\n  Saved PDF: {pdf_path}")
    print(f"  Saved PNG: {png_path}")
    print()
    print("F1 v2 done. Open the PNG, review, then send the image back.")


if __name__ == "__main__":
    main()