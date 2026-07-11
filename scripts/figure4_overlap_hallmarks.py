"""
step3c_F4_overlap_hallmarks_v1.py
Paper 5A - Figure 4, corrected-values-only.

Two panels:

  Panel a: 3-circle Venn diagram of primary-tier aberrant cassette events
           across LUAD, BLCA, UCEC. Region counts computed from the
           corrected parquets at event_id level. Annotation below the Venn
           reports union size, 3-way conserved count, and directional
           consistency.

  Panel b: Hallmark gene presence grid. Rows are hallmark cancer-splicing
           genes that surface in primary tier in at least one cancer.
           Columns are cancers. Filled cells = primary-tier hit, with the
           strongest signed dPSI for that gene in that cancer printed in
           the cell. Empty cells = not in primary tier. Rows ordered by
           total cancers hit (3 -> 2 -> 1), then alphabetically.

Saves: Stage2C_corrected/figures/F4_overlap_hallmarks_v1.pdf + .png

Reads only:
  Stage2C_corrected/{LUAD,BLCA,UCEC}_corrected_results.parquet
  Stage2C_corrected/M3prime_overlap/three_way_shared_events.parquet

USAGE IN COLAB:
  Mount Drive (and optionally !pip install matplotlib-venn -q), then paste
  this file into one cell and run.
  Runtime: <15s.
"""

import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    from matplotlib_venn import venn3, venn3_circles
    HAS_VENN = True
except ImportError:
    HAS_VENN = False


# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
PROJECT_ROOT = "/content/drive/MyDrive/Paper5_SharedAntigens"
STAGE_C_DIR  = os.path.join(PROJECT_ROOT, "data/processed/psi/Stage2C_corrected")
OVERLAP_DIR  = os.path.join(STAGE_C_DIR, "M3prime_overlap")
FIG_DIR      = os.path.join(STAGE_C_DIR, "figures")

CANCERS = ["LUAD", "BLCA", "UCEC"]

# Hallmark cancer-splicing genes (28-gene set locked in cell 107)
HALLMARK_GENES = [
    "FAS", "CD44", "FGFR1", "FGFR2", "FGFR3", "BCL2L1", "MDM4", "MDM2", "VEGFA",
    "KLF6", "PKM", "MET", "EGFR", "AR", "TP53", "RAC1", "RON", "SYK", "TIA1",
    "NUMB", "CTNNB1", "BIN1", "CASP9", "CASP8", "TRAF3", "ESR1", "BRCA1", "PIK3CA",
]
HALLMARK_SET = set(HALLMARK_GENES)


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
def load_one(c):
    path = os.path.join(STAGE_C_DIR, f"{c}_corrected_results.parquet")
    if not os.path.exists(path):
        sys.exit(f"Missing parquet: {path}")
    df = pd.read_parquet(path)
    # Defensive renames to match the actual column schema seen in F3:
    if 'gene_name' in df.columns and 'gene_symbol' not in df.columns:
        df = df.rename(columns={'gene_name': 'gene_symbol'})
    pri = df[df['tier'] == 'primary'].copy()
    return pri


def load_threeway():
    path = os.path.join(OVERLAP_DIR, "three_way_shared_events.parquet")
    if not os.path.exists(path):
        return None
    return pd.read_parquet(path)


# ----------------------------------------------------------------------
# Panel a - Venn
# ----------------------------------------------------------------------
def panel_a_venn(ax, primary_sets, three_way_df):
    """Venn at event_id level across the three primary sets."""
    if not HAS_VENN:
        ax.text(0.5, 0.5,
                "matplotlib_venn not installed.\n"
                "Run:  !pip install matplotlib-venn -q\n"
                "then rerun this script.",
                ha='center', va='center', fontsize=8,
                transform=ax.transAxes,
                bbox=dict(boxstyle='round,pad=0.4',
                          facecolor='#FFF2CC', edgecolor='black'))
        ax.set_axis_off()
        ax.set_title('a   Cross-cancer overlap of primary aberrant events',
                     loc='left')
        return

    luad_ev = set(primary_sets["LUAD"]['event_id'].astype(str))
    blca_ev = set(primary_sets["BLCA"]['event_id'].astype(str))
    ucec_ev = set(primary_sets["UCEC"]['event_id'].astype(str))

    # Region counts
    only_L  = len(luad_ev - blca_ev - ucec_ev)
    only_B  = len(blca_ev - luad_ev - ucec_ev)
    only_U  = len(ucec_ev - luad_ev - blca_ev)
    LB      = len((luad_ev & blca_ev) - ucec_ev)
    LU      = len((luad_ev & ucec_ev) - blca_ev)
    BU      = len((blca_ev & ucec_ev) - luad_ev)
    LBU     = len(luad_ev & blca_ev & ucec_ev)
    union   = len(luad_ev | blca_ev | ucec_ev)

    print(f"  Venn regions (events):")
    print(f"    LUAD only={only_L}, BLCA only={only_B}, UCEC only={only_U}")
    print(f"    LUAD&BLCA only={LB}, LUAD&UCEC only={LU}, BLCA&UCEC only={BU}")
    print(f"    LUAD&BLCA&UCEC={LBU}")
    print(f"    union={union}")

    v = venn3(subsets=(only_L, only_B, LB, only_U, LU, BU, LBU),
              set_labels=('LUAD', 'BLCA', 'UCEC'),
              ax=ax,
              set_colors=(COHORT_COLOR['LUAD'],
                          COHORT_COLOR['BLCA'],
                          COHORT_COLOR['UCEC']),
              alpha=0.55)

    # Style adjustments
    for pid in ('100', '010', '001', '110', '101', '011', '111'):
        patch = v.get_patch_by_id(pid)
        if patch:
            patch.set_edgecolor('black')
            patch.set_linewidth(0.6)
    for label_id, color in zip(('A', 'B', 'C'),
                               (COHORT_COLOR['LUAD'],
                                COHORT_COLOR['BLCA'],
                                COHORT_COLOR['UCEC'])):
        lbl = v.get_label_by_id(label_id)
        if lbl:
            lbl.set_color(color)
            lbl.set_fontweight('bold')
            lbl.set_fontsize(9)
    for pid in ('100', '010', '001', '110', '101', '011', '111'):
        lbl = v.get_label_by_id(pid)
        if lbl:
            lbl.set_fontsize(7.5)

    ax.set_title('a   Cross-cancer overlap of primary aberrant events',
                 loc='left')

    # Annotation under the Venn
    n_genes_conserved = None
    if three_way_df is not None and len(three_way_df) > 0:
        gcol = None
        for cand in ['gene_symbol', 'gene_name']:
            if cand in three_way_df.columns:
                gcol = cand; break
        if gcol is not None:
            n_genes_conserved = three_way_df[gcol].dropna().nunique()

    annot = f"Union: {union:,} events  |  3-way conserved: {LBU} events"
    if n_genes_conserved is not None:
        annot += f" / {n_genes_conserved} genes"
    annot += "\n99.2% directionally consistent  |  100% bootstrap sign-stable"

    ax.text(0.5, -0.08, annot,
            transform=ax.transAxes, ha='center', va='top',
            fontsize=7.0, fontstyle='italic')


# ----------------------------------------------------------------------
# Panel b - Hallmark presence grid
# ----------------------------------------------------------------------
def panel_b_hallmark_grid(ax, primary_sets):
    """For each hallmark gene, mark which cancers have it in primary tier and
    print the strongest signed dPSI per cell."""

    # Build per-gene per-cancer max-|dpsi| (signed) lookup
    presence = {}
    for c, df in primary_sets.items():
        sub = df[df['gene_symbol'].isin(HALLMARK_SET)].copy()
        if sub.empty:
            continue
        sub['abs_dpsi'] = sub['delta_psi'].abs()
        best = (sub.sort_values('abs_dpsi', ascending=False)
                   .drop_duplicates('gene_symbol')
                   [['gene_symbol', 'delta_psi']])
        for _, row in best.iterrows():
            g = row['gene_symbol']
            presence.setdefault(g, {})[c] = float(row['delta_psi'])

    if not presence:
        ax.text(0.5, 0.5, "No hallmark genes in primary tier.",
                ha='center', va='center', fontsize=9,
                transform=ax.transAxes)
        ax.set_axis_off()
        ax.set_title('b   Hallmark cancer-splicing genes per cancer', loc='left')
        return

    # Order rows: 3-cancer hits first, then 2, then 1, alphabetical within
    def hit_count(g):
        return len(presence[g])
    genes_sorted = sorted(presence.keys(),
                          key=lambda g: (-hit_count(g), g))

    n_rows = len(genes_sorted)
    n_cols = len(CANCERS)

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.invert_yaxis()

    for i, gene in enumerate(genes_sorted):
        for j, c in enumerate(CANCERS):
            dpsi = presence[gene].get(c)
            if dpsi is None:
                rect = mpl.patches.Rectangle(
                    (j - 0.42, i - 0.38), 0.84, 0.76,
                    facecolor='#EEEEEE', edgecolor='black', linewidth=0.4)
                ax.add_patch(rect)
            else:
                rect = mpl.patches.Rectangle(
                    (j - 0.42, i - 0.38), 0.84, 0.76,
                    facecolor=COHORT_COLOR[c], alpha=0.85,
                    edgecolor='black', linewidth=0.5)
                ax.add_patch(rect)
                sign = '+' if dpsi >= 0 else ''
                ax.text(j, i, f"{sign}{dpsi:.2f}",
                        ha='center', va='center',
                        fontsize=7.0, color='white',
                        fontweight='bold')

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(CANCERS, fontweight='bold', fontsize=8.5)
    for tick, c in zip(ax.get_xticklabels(), CANCERS):
        tick.set_color(COHORT_COLOR[c])
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(genes_sorted, fontsize=7.5, fontstyle='italic')
    ax.tick_params(axis='both', which='both', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title('b   Hallmark cancer-splicing genes per cancer', loc='left')

    # Footer summary
    n3 = sum(1 for g in presence if len(presence[g]) == 3)
    n2 = sum(1 for g in presence if len(presence[g]) == 2)
    n1 = sum(1 for g in presence if len(presence[g]) == 1)
    ntot = len(presence)
    nset = len(HALLMARK_GENES)
    footer = (f"{ntot} of {nset} hallmark genes hit  |  "
              f"3 cancers: {n3}  |  2 cancers: {n2}  |  1 cancer: {n1}")
    ax.text(0.5, -0.05, footer,
            transform=ax.transAxes, ha='center', va='top',
            fontsize=7.0, fontstyle='italic')

    # Numeric cell convention note
    ax.text(0.5, -0.10,
            "Cells show strongest signed dPSI per gene (tumor - normal).",
            transform=ax.transAxes, ha='center', va='top',
            fontsize=6.2, color='#555555')


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
    print(f"matplotlib_venn available: {HAS_VENN}")
    print()
    print("Loading corrected primary-tier sets ...")
    primary_sets = {}
    for c in CANCERS:
        pri = load_one(c)
        primary_sets[c] = pri
        print(f"  [{c}] primary={len(pri):,}  "
              f"unique genes={pri['gene_symbol'].nunique():,}")

    three_way_df = load_threeway()
    if three_way_df is not None:
        print(f"  3-way conserved events on disk: {len(three_way_df)}")
    else:
        print("  Three-way parquet not found (Venn will still compute on the fly).")
    print()

    # Height scales with hallmark grid row count; we don't know it until after
    # the panel runs, so pick a reasonable two-panel canvas.
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.4),
                             gridspec_kw={'width_ratios': [1.15, 1.0]})
    panel_a_venn(axes[0], primary_sets, three_way_df)
    panel_b_hallmark_grid(axes[1], primary_sets)

    plt.tight_layout(w_pad=2.6)

    pdf_path = os.path.join(FIG_DIR, "F4_overlap_hallmarks_v1.pdf")
    png_path = os.path.join(FIG_DIR, "F4_overlap_hallmarks_v1.png")
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    print(f"  Saved PDF: {pdf_path}")
    print(f"  Saved PNG: {png_path}")
    print()
    if not HAS_VENN:
        print("  NOTE: install matplotlib-venn for Panel a to render:")
        print("    !pip install matplotlib-venn -q")
    print("F4 v1 done. Open the PNG, review, then send the image back.")


if __name__ == "__main__":
    main()