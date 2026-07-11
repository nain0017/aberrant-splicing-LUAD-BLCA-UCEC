"""
step3c_F5_hallmark_case_studies_v1.py
Paper 5A - Figure 5, corrected-values-only.

Six-panel grid of hallmark cancer-splicing genes showing per-sample PSI
distributions in tumor vs normal, for each cancer where the gene surfaces
in primary tier.

Panels are picked from the corrected primary-tier sets across the three
cancers. For each candidate gene, the script picks the event_id with the
strongest |dPSI| per cancer and pulls the per-sample PSI vector from the
cell-60 npz files. Mean PSI lines are drawn per group (matching the
corrected pipeline convention of mean, not median).

Saves: Stage2C_corrected/figures/F5_hallmark_case_studies_v1.pdf + .png

Reads:
  Stage2C_corrected/{LUAD,BLCA,UCEC}_corrected_results.parquet
  data/processed/suppa/suppa_events_parsed.parquet      (event ordering)
  data/processed/psi/{LUAD,BLCA,UCEC}_tumor_psi.npz     (per-sample PSI)
  data/processed/psi/{LUNG,BLADDER,UTERUS}_normal_psi.npz

USAGE IN COLAB:
  Paste this file into one cell and run after mounting Drive.
  Runtime: ~30-40s.

NOTE on input paths:
  npz filenames vary across pipeline iterations. The script tries several
  plausible filename patterns and reports which it loaded. If none match,
  it prints the candidate filenames and exits cleanly so we can patch.
"""

import os
import sys
import glob
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
PSI_DIR      = os.path.join(PROJECT_ROOT, "data/processed/psi")
SUPPA_DIR    = PSI_DIR  # suppa_events_parsed.parquet lives in data/processed/psi/
FIG_DIR      = os.path.join(STAGE_C_DIR, "figures")

CANCERS = ["LUAD", "BLCA", "UCEC"]
NORMALS = {"LUAD": "LUNG", "BLCA": "BLADDER", "UCEC": "UTERUS"}

# Candidate hallmarks to consider for the 6 panels. The actual 6 picked
# is determined from the data: prefer genes primary in >=2 cancers, then by
# total |dPSI|. This keeps the selection data-driven and reproducible.
HALLMARK_CANDIDATES = [
    "CD44", "NUMB", "FAS", "FGFR2", "FGFR3", "VEGFA", "MDM4", "SRSF2",
    "FGFR1", "SYK", "BIN1", "RON", "CASP8", "BCL2L1",
]
N_PANELS = 6


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
    "LUAD": "#2E86AB",
    "BLCA": "#E07A5F",
    "UCEC": "#76A646",
}


def apply_style():
    for k, v in STYLE.items():
        mpl.rcParams[k] = v


# ----------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------
def find_npz(label, kind):
    """kind in {'tumor', 'normal'}. Tries common filename patterns."""
    patterns = [
        f"{label}_{kind}_psi.npz",
        f"{label}_psi.npz",
        f"{label.lower()}_{kind}_psi.npz",
        f"{label.lower()}_psi.npz",
        f"{label}_{kind}.npz",
    ]
    for p in patterns:
        cand = os.path.join(PSI_DIR, p)
        if os.path.exists(cand):
            return cand
    # Fallback: glob for any file starting with this label
    matches = glob.glob(os.path.join(PSI_DIR, f"{label}*psi*.npz"))
    if matches:
        return matches[0]
    return None


def load_psi_matrix(label, kind):
    """Returns a dict with at least 'psi' (n_events x n_samples)."""
    path = find_npz(label, kind)
    if path is None:
        print(f"  [{label}/{kind}] npz NOT FOUND in {PSI_DIR}")
        return None
    data = np.load(path, allow_pickle=True)
    keys = list(data.keys())
    print(f"  [{label}/{kind}] loaded {os.path.basename(path)}  keys={keys}")
    if 'psi' not in keys:
        print(f"    No 'psi' key. Available: {keys}")
        return None
    psi = data['psi']  # shape (n_events, n_samples)
    return psi


def load_suppa_events():
    path = os.path.join(SUPPA_DIR, "suppa_events_parsed.parquet")
    if not os.path.exists(path):
        sys.exit(f"Missing SUPPA events parquet: {path}")
    df = pd.read_parquet(path)
    print(f"  SUPPA events: {len(df):,} rows; columns={list(df.columns)[:8]} ...")
    return df


def load_corrected_primary(c):
    path = os.path.join(STAGE_C_DIR, f"{c}_corrected_results.parquet")
    df = pd.read_parquet(path)
    if 'gene_name' in df.columns and 'gene_symbol' not in df.columns:
        df = df.rename(columns={'gene_name': 'gene_symbol'})
    return df[df['tier'] == 'primary'].copy()


# ----------------------------------------------------------------------
# Gene + event selection
# ----------------------------------------------------------------------
def pick_genes_and_events(primary_sets):
    """For each candidate hallmark gene, find which cancers it's primary in
    and the event_id with strongest |dPSI| in each. Then pick top N_PANELS
    by total |dPSI| across cancers, breaking ties by number of cancers hit."""
    info = {}  # gene -> {cancer: {event_id, delta_psi, abs_dpsi}}
    for g in HALLMARK_CANDIDATES:
        for c in CANCERS:
            df = primary_sets[c]
            sub = df[df['gene_symbol'] == g].copy()
            if sub.empty:
                continue
            sub['abs_dpsi'] = sub['delta_psi'].abs()
            best = sub.sort_values('abs_dpsi', ascending=False).iloc[0]
            info.setdefault(g, {})[c] = {
                'event_id':  str(best['event_id']),
                'delta_psi': float(best['delta_psi']),
                'abs_dpsi':  float(best['abs_dpsi']),
            }

    if not info:
        sys.exit("No hallmark candidates found in any primary set.")

    # Rank genes: more cancers hit > total |dPSI|
    def score(g):
        n_cancers = len(info[g])
        total_abs = sum(d['abs_dpsi'] for d in info[g].values())
        return (n_cancers, total_abs)
    genes_ranked = sorted(info.keys(), key=score, reverse=True)
    picked = genes_ranked[:N_PANELS]

    print(f"\n  Hallmark gene selection (top {N_PANELS}):")
    for g in picked:
        per_c = info[g]
        cc = ', '.join(f"{c}:{d['delta_psi']:+.2f}"
                       for c, d in per_c.items())
        print(f"    {g}  [{len(per_c)} cancer(s)]  {cc}")
    print()

    return {g: info[g] for g in picked}


# ----------------------------------------------------------------------
# PSI vector retrieval
# ----------------------------------------------------------------------
def build_event_index(suppa_df):
    """event_id -> row position in SUPPA events. Try multiple id columns."""
    for col in ['event_id', 'event', 'event_name', 'id']:
        if col in suppa_df.columns:
            id_col = col; break
    else:
        sys.exit(f"No event id column in SUPPA parquet. Got: {list(suppa_df.columns)}")
    print(f"  Event id column in SUPPA: {id_col}")
    return {str(eid): i for i, eid in enumerate(suppa_df[id_col])}, id_col


def get_psi_vector(psi_matrix, event_idx):
    if psi_matrix is None or event_idx is None:
        return None
    if event_idx >= psi_matrix.shape[0]:
        return None
    v = psi_matrix[event_idx, :].astype(float)
    # Drop NaN samples for strip plotting
    return v[~np.isnan(v)]


# ----------------------------------------------------------------------
# Per-panel plotter
# ----------------------------------------------------------------------
def plot_gene_panel(ax, gene, gene_info, event_id_to_idx, psi_tumor, psi_normal):
    """Strip plot of per-sample PSI for tumor and normal in each cancer where
    the gene is primary. mean lines drawn per group, light line connecting
    T-mean to N-mean within each cancer."""
    cancers_present = [c for c in CANCERS if c in gene_info]

    # Columns: [T_LUAD, N_LUAD, T_BLCA, N_BLCA, ...] in cancers_present order
    x_positions = []
    x_labels    = []
    pair_centers = []
    deltas = {}

    col_x = 0
    for c in cancers_present:
        eid = gene_info[c]['event_id']
        idx = event_id_to_idx.get(eid)
        if idx is None:
            print(f"    [{gene}/{c}] event_id {eid} not found in SUPPA index")
            continue

        t_vec = get_psi_vector(psi_tumor[c], idx)
        n_vec = get_psi_vector(psi_normal[c], idx)
        if t_vec is None or n_vec is None or len(t_vec) == 0 or len(n_vec) == 0:
            print(f"    [{gene}/{c}] empty PSI vectors at event_id {eid}")
            continue

        color = COHORT_COLOR[c]

        # Tumor strip
        xt = col_x + np.random.uniform(-0.18, 0.18, size=len(t_vec))
        ax.scatter(xt, t_vec, s=4, c=color, alpha=0.30,
                   edgecolors='none', rasterized=True)
        # Normal strip
        xn = (col_x + 1) + np.random.uniform(-0.18, 0.18, size=len(n_vec))
        ax.scatter(xn, n_vec, s=4, c=color, alpha=0.30,
                   edgecolors='none', rasterized=True)

        # Mean lines
        mt = float(np.mean(t_vec))
        mn = float(np.mean(n_vec))
        ax.hlines(mt, col_x - 0.30, col_x + 0.30,
                  colors='black', linewidth=1.6, zorder=5)
        ax.hlines(mn, col_x + 1 - 0.30, col_x + 1 + 0.30,
                  colors='black', linewidth=1.6, zorder=5)

        # Connector
        ax.plot([col_x, col_x + 1], [mt, mn],
                color='grey', linewidth=0.6, alpha=0.7, zorder=4)

        x_positions.extend([col_x, col_x + 1])
        x_labels.extend(['T', 'N'])
        pair_centers.append((col_x + 0.5, c, mt - mn))
        deltas[c] = mt - mn

        col_x += 2.6  # gap between cancer groups

    # Strip-level styling
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=7.0)
    ax.set_ylim(-0.04, 1.04)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel('PSI', fontsize=7.5)

    # Cancer labels under each pair
    for center, c, _ in pair_centers:
        ax.text(center, -0.18, c,
                ha='center', va='top', fontsize=8.0,
                color=COHORT_COLOR[c], fontweight='bold',
                transform=ax.get_xaxis_transform())

    # Per-cancer dPSI text above each pair
    for center, c, dpsi in pair_centers:
        sign = '+' if dpsi >= 0 else ''
        ax.text(center, 1.06, f"Δ={sign}{dpsi:.2f}",
                ha='center', va='bottom', fontsize=6.8,
                color='black')

    # Title: gene name only, italic, centered
    ax.set_title(gene, loc='center', fontstyle='italic', fontweight='bold',
                 pad=18)

    # Tidy: shrink x range
    if x_positions:
        ax.set_xlim(min(x_positions) - 0.6, max(x_positions) + 0.6)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    if not os.path.exists("/content/drive/MyDrive"):
        sys.exit("Drive not mounted. Run drive.mount('/content/drive') first.")
    if not os.path.exists(STAGE_C_DIR):
        sys.exit(f"Stage2C_corrected not found: {STAGE_C_DIR}")
    os.makedirs(FIG_DIR, exist_ok=True)

    np.random.seed(42)  # jitter reproducibility
    apply_style()
    print(f"Output: {FIG_DIR}")
    print(f"Date:   {date.today().isoformat()}")
    print()

    # Primary tier per cancer
    print("Loading corrected primary-tier sets ...")
    primary_sets = {c: load_corrected_primary(c) for c in CANCERS}
    for c, df in primary_sets.items():
        print(f"  [{c}] primary={len(df):,}")
    print()

    # Pick six hallmark genes from the data
    print("Selecting hallmark genes for case study panels ...")
    picked = pick_genes_and_events(primary_sets)

    # SUPPA events for event_id -> matrix-row mapping
    print("Loading SUPPA events parquet ...")
    suppa = load_suppa_events()
    event_id_to_idx, _ = build_event_index(suppa)
    print()

    # PSI matrices (lazy: only load cancers we actually need)
    needed_cancers = set()
    for g, info in picked.items():
        needed_cancers.update(info.keys())
    print(f"Loading PSI matrices for cancers: {sorted(needed_cancers)} ...")
    psi_tumor  = {c: load_psi_matrix(c, 'tumor')          for c in needed_cancers}
    psi_normal = {c: load_psi_matrix(NORMALS[c], 'normal') for c in needed_cancers}
    if any(v is None for v in psi_tumor.values()) or any(v is None for v in psi_normal.values()):
        sys.exit("\nOne or more PSI npz files could not be loaded. See messages above.")
    print()

    # Figure
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 6.2))
    axes_flat = axes.flatten()
    for i, (gene, info) in enumerate(picked.items()):
        plot_gene_panel(axes_flat[i], gene, info,
                        event_id_to_idx, psi_tumor, psi_normal)

    # Panel letters
    for ax, letter in zip(axes_flat, 'abcdef'):
        ax.text(-0.10, 1.06, letter, transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom', ha='left')

    plt.tight_layout(w_pad=2.6, h_pad=3.4)

    pdf_path = os.path.join(FIG_DIR, "F5_hallmark_case_studies_v1.pdf")
    png_path = os.path.join(FIG_DIR, "F5_hallmark_case_studies_v1.png")
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    print(f"  Saved PDF: {pdf_path}")
    print(f"  Saved PNG: {png_path}")
    print()
    print("F5 v1 done. Open the PNG, review, then send the image back.")


if __name__ == "__main__":
    main()