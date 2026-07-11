"""
step3c_F3_volcanoes_v1.py
Paper 5A — Figure 3, corrected-values-only.

Three panels: LUAD, BLCA, UCEC volcano plots from the corrected beta-binomial
LRT. Signed dPSI (mean-based) on x; -log10(q_BB) on y. Coverage-pass events
that fail gates are grey low-alpha; events passing q<0.05 AND |dPSI|>=0.15
are cancer-colored. Primary-tier (3C-validated) events get a black edge ring.
Hallmark genes that surface in primary tier are labeled.

Saves: Stage2C_corrected/figures/F3_volcanoes_v1.pdf + .png

Reads only:
  Stage2C_corrected/{LUAD,BLCA,UCEC}_corrected_results.parquet

USAGE IN COLAB:
  Paste this file's contents into one cell and run. Run the Drive-mount
  block first if Drive is not yet mounted.
  Runtime: ~30-60s depending on adjustText availability.
"""

import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    from adjustText import adjust_text
    HAS_ADJUST = True
except ImportError:
    HAS_ADJUST = False


# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
PROJECT_ROOT = "/content/drive/MyDrive/Paper5_SharedAntigens"
STAGE_C_DIR  = os.path.join(PROJECT_ROOT, "data/processed/psi/Stage2C_corrected")
FIG_DIR      = os.path.join(STAGE_C_DIR, "figures")

CANCERS = ["LUAD", "BLCA", "UCEC"]

# Per-cancer sample sizes (T vs N) for subtitles
SAMPLE_SIZES = {
    "LUAD": (542, 655),
    "BLCA": (414, 21),
    "UCEC": (554, 159),
}

# Hallmark genes locked in cell 107 (28 genes)
HALLMARK_GENES = {
    "FAS", "CD44", "FGFR1", "FGFR2", "FGFR3", "BCL2L1", "MDM4", "MDM2", "VEGFA",
    "KLF6", "PKM", "MET", "EGFR", "AR", "TP53", "RAC1", "RON", "SYK", "TIA1",
    "NUMB", "CTNNB1", "BIN1", "CASP9", "CASP8", "TRAF3", "ESR1", "BRCA1", "PIK3CA",
}

# Gate thresholds (locked)
Q_THR     = 0.05
DPSI_THR  = 0.15


# ----------------------------------------------------------------------
# Style
# ----------------------------------------------------------------------
STYLE = {
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
    'legend.fontsize':       6.8,
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
}


def apply_style():
    for k, v in STYLE.items():
        mpl.rcParams[k] = v


# ----------------------------------------------------------------------
# Data loading + prep
# ----------------------------------------------------------------------
def load_one(c):
    path = os.path.join(STAGE_C_DIR, f"{c}_corrected_results.parquet")
    if not os.path.exists(path):
        sys.exit(f"Missing parquet: {path}")
    df = pd.read_parquet(path)
    # Keep only coverage-pass events for the volcano (others have no test
    # result we can plot). This matches what the betabin actually tested.
    cov = df['coverage_pass'].fillna(False).astype(bool)
    df = df[cov].copy()

    # Required columns: BB q-value (Stage 2C corrected pipeline writes 'q_value'
    # from the BB LRT + BH-FDR), dpsi (mean-based, signed), tier, gene_symbol.
    # Fall back across naming variants defensively.
    qcol = None
    for cand in ['q_value', 'bb_qvalue', 'qvalue_bb', 'q_bb', 'qvalue']:
        if cand in df.columns:
            qcol = cand; break
    if qcol is None:
        sys.exit(f"[{c}] cannot find BB q-value column. Got: {list(df.columns)}")

    dpcol = None
    for cand in ['delta_psi', 'dpsi', 'mean_dpsi']:
        if cand in df.columns:
            dpcol = cand; break
    if dpcol is None:
        sys.exit(f"[{c}] cannot find signed dPSI column. Got: {list(df.columns)}")

    gsym = None
    for cand in ['gene_symbol', 'gene_name', 'symbol', 'hgnc_symbol']:
        if cand in df.columns:
            gsym = cand; break

    df = df.rename(columns={qcol: 'q', dpcol: 'dpsi'})
    if gsym and gsym != 'gene_symbol':
        df = df.rename(columns={gsym: 'gene_symbol'})
    elif gsym is None:
        df['gene_symbol'] = ''

    # Clean
    df['q'] = df['q'].astype(float).clip(lower=1e-300)  # avoid log(0)
    df['neglog10q'] = -np.log10(df['q'])
    df['dpsi'] = df['dpsi'].astype(float)
    df['gene_symbol'] = df['gene_symbol'].fillna('').astype(str)

    return df


def label_hallmark_in_primary(df_primary):
    """Pick one row per hallmark gene to label: the row with max |dpsi|."""
    if df_primary.empty:
        return df_primary.iloc[0:0]
    hits = df_primary[df_primary['gene_symbol'].isin(HALLMARK_GENES)].copy()
    if hits.empty:
        return hits
    hits['abs_dpsi'] = hits['dpsi'].abs()
    # one row per gene = the strongest by |dpsi|
    hits = (hits.sort_values('abs_dpsi', ascending=False)
                .drop_duplicates('gene_symbol')
                .reset_index(drop=True))
    return hits


# ----------------------------------------------------------------------
# Per-panel plotter
# ----------------------------------------------------------------------
def plot_volcano(ax, c, df, y_cap):
    color = COHORT_COLOR[c]

    # Slice events into visual categories
    cov_total   = len(df)
    passes_q    = df['q'] < Q_THR
    passes_dp   = df['dpsi'].abs() >= DPSI_THR
    is_primary  = df['tier'] == 'primary'
    is_signif   = passes_q & passes_dp

    # Capped y for display only
    y = df['neglog10q'].clip(upper=y_cap)
    n_at_cap = int((df['neglog10q'] > y_cap).sum())

    # Background grey (everything not passing both gates)
    grey_idx = ~is_signif
    ax.scatter(df.loc[grey_idx, 'dpsi'], y[grey_idx],
               s=4, c='#CCCCCC', alpha=0.45,
               edgecolors='none', rasterized=True)

    # Significant non-primary (secondary or scn-discordant)
    sig_nonpri = is_signif & ~is_primary
    ax.scatter(df.loc[sig_nonpri, 'dpsi'], y[sig_nonpri],
               s=6, c=color, alpha=0.55,
               edgecolors='none', rasterized=True)

    # Primary tier on top with black ring
    pri = is_primary
    ax.scatter(df.loc[pri, 'dpsi'], y[pri],
               s=10, c=color, alpha=0.95,
               edgecolors='black', linewidths=0.25,
               rasterized=True)

    # Gate lines
    ax.axhline(-np.log10(Q_THR), color='black', linestyle=':',
               linewidth=0.5, alpha=0.5)
    ax.axvline( DPSI_THR, color='black', linestyle=':',
               linewidth=0.5, alpha=0.5)
    ax.axvline(-DPSI_THR, color='black', linestyle=':',
               linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='black', linewidth=0.4, alpha=0.4)

    # Hallmark labels (drawn from primary tier only) — vertical orientation
    # to avoid overlap when many hallmarks pile up at the y-cap line.
    df_pri = df[pri].copy()
    df_pri['neglog10q_capped'] = df_pri['neglog10q'].clip(upper=y_cap)
    hits = label_hallmark_in_primary(df_pri)

    # Sort by x so labels read in dPSI order along the cap; helps visual scan
    hits = hits.sort_values('dpsi').reset_index(drop=True)

    # Vertical offset above each marker so the text starts above the point
    # and reads upward. Tweak label_offset_frac if labels still touch markers.
    label_offset_frac = 0.025  # fraction of y_cap above each marker
    label_offset = label_offset_frac * y_cap

    for _, row in hits.iterrows():
        # Purple ring marker
        ax.scatter([row['dpsi']], [row['neglog10q_capped']],
                   s=22, facecolor='#C5A3D8', edgecolor='black',
                   linewidth=0.5, zorder=5)
        # Vertical label, anchored bottom so it grows upward from the marker
        ax.text(row['dpsi'],
                row['neglog10q_capped'] + label_offset,
                row['gene_symbol'],
                fontsize=6.5, fontstyle='italic', zorder=6,
                color='black',
                rotation=90,
                rotation_mode='anchor',
                ha='center', va='bottom')

    # Title with sample sizes
    nT, nN = SAMPLE_SIZES[c]
    ax.set_title(f"{c}  (n={nT} vs {nN})",
                 loc='center', color=color, fontweight='bold')

    # Axis labels
    ax.set_xlabel(r'$\Delta$PSI (tumor - normal)')
    ax.set_ylabel(r'$-\log_{10}(q)$')
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(0, y_cap * 1.04)

    # Counts inset
    n_primary   = int(pri.sum())
    n_secondary = int(((df['tier'] == 'secondary')).sum())
    inset = f"primary: {n_primary:,}\nsecondary: {n_secondary:,}"
    if n_at_cap > 0:
        inset += f"\n({n_at_cap:,} at y-cap)"
    ax.text(0.95, 0.90, inset, transform=ax.transAxes,
            ha='right', va='top', fontsize=6.5,
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
    print(f"Output: {FIG_DIR}")
    print(f"Date:   {date.today().isoformat()}")
    print(f"adjustText available: {HAS_ADJUST}")
    print()
    print("Loading corrected results ...")
    data = {c: load_one(c) for c in CANCERS}
    for c, df in data.items():
        n_pri = int((df['tier'] == 'primary').sum())
        n_sec = int((df['tier'] == 'secondary').sum())
        max_y = df['neglog10q'].replace(np.inf, np.nan).max()
        p99   = df['neglog10q'].replace(np.inf, np.nan).quantile(0.99)
        print(f"  [{c}] coverage-pass={len(df):,}  primary={n_pri}  "
              f"secondary={n_sec}  max(-log10q)={max_y:.1f}  "
              f"p99(-log10q)={p99:.1f}")

    # Shared y-cap across panels: max of per-cancer p99, capped at 50
    p99s = [data[c]['neglog10q'].replace(np.inf, np.nan).quantile(0.99)
            for c in CANCERS]
    y_cap = float(min(50.0, max(p99s) * 1.05))
    y_cap = max(y_cap, 10.0)  # never less than 10 for readability
    print(f"\n  Shared y-cap (display only): {y_cap:.1f}")

    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.6), sharey=True)
    for ax, c in zip(axes, CANCERS):
        plot_volcano(ax, c, data[c], y_cap)

    plt.tight_layout(w_pad=1.6)

    pdf_path = os.path.join(FIG_DIR, "F3_volcanoes_v1.pdf")
    png_path = os.path.join(FIG_DIR, "F3_volcanoes_v1.png")
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    print(f"\n  Saved PDF: {pdf_path}")
    print(f"  Saved PNG: {png_path}")
    print()
    if not HAS_ADJUST:
        print("  NOTE: adjustText not installed; labels may overlap. "
              "Install with: !pip install adjustText -q")
    print("F3 v1 done. Open the PNG, review, then send the image back.")


if __name__ == "__main__":
    main()