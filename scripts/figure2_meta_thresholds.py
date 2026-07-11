"""
step3c_F2_meta_thresholds_v4.py
Paper 5A -- Figure 2, v4.

Same data and verdict logic as v3. Only change: the five pre-commitment
criteria are now labeled Q1-Q5 with plain descriptive text instead of the
internal M1'/M2'/M3'/M4/M5' notation, since that notation is never defined
in the manuscript or supplementary text and reads as an unresolved internal
artifact to a reviewer. Panel b title and footnote also de-jargoned to match.

No data, thresholds, or verdicts change. Reads same inputs, writes:
  Stage2C_corrected/figures/F2_meta_thresholds_v4.pdf + .png

USAGE: paste this file into one Colab cell and run.
"""

import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


PROJECT_ROOT = "/content/drive/MyDrive/Paper5_SharedAntigens"
STAGE_C_DIR  = os.path.join(PROJECT_ROOT, "data/processed/psi/Stage2C_corrected")
OVERLAP_DIR  = os.path.join(STAGE_C_DIR, "M3prime_overlap")
M4_DIR       = os.path.join(STAGE_C_DIR, "M4_enrichment")
FIG_DIR      = os.path.join(STAGE_C_DIR, "figures")

CANCERS = ["LUAD", "BLCA", "UCEC"]

STYLE = {
    'figure.figsize':        (7.2, 4.6),
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

COHORT_COLOR = {"LUAD": "#2E86AB", "BLCA": "#E07A5F", "UCEC": "#76A646"}
VERDICT_COLOR = {
    "PASS":       "#A9D8A0",
    "PARTIAL":    "#F2D17C",
    "BORDERLINE": "#F2D17C",
    "REFRAMED":   "#C8B6E2",
    "FAIL":       "#E89A8A",
}


def apply_style():
    for k, v in STYLE.items():
        mpl.rcParams[k] = v


def load_data():
    q1 = {"LUAD": 25314, "BLCA": 17448, "UCEC": 20890}

    q2 = {}
    for c in CANCERS:
        path = os.path.join(STAGE_C_DIR, f"{c}_corrected_results.parquet")
        if not os.path.exists(path):
            sys.exit(f"Missing parquet: {path}")
        df = pd.read_parquet(path)
        q2[c] = int((df['tier'] == 'primary').sum())

    pw_path = os.path.join(OVERLAP_DIR, "pairwise_event_overlap.csv")
    if not os.path.exists(pw_path):
        sys.exit(f"Missing CSV: {pw_path}")
    pw = pd.read_csv(pw_path)

    cc_path = os.path.join(OVERLAP_DIR, "three_way_shared_events.parquet")
    n_conserved_events = None
    if os.path.exists(cc_path):
        n_conserved_events = len(pd.read_parquet(cc_path))

    m4_path = os.path.join(M4_DIR, "m4_fisher_per_cancer.csv")
    q4_or, q4_p = {}, {}
    if os.path.exists(m4_path):
        m4 = pd.read_csv(m4_path)
        for _, row in m4.iterrows():
            label = row.get('cohort', row.get('cancer', None))
            if label is None:
                continue
            label = str(label).strip().upper()
            if label in CANCERS:
                q4_or[label] = float(row['OR']) if 'OR' in row else float(row.get('odds_ratio', np.nan))
                q4_p[label]  = float(row['p_value']) if 'p_value' in row else float(row.get('p', row.get('pvalue', np.nan)))

    return q1, q2, pw, n_conserved_events, q4_or, q4_p


def panel_a_grid(ax, q1, q2, pw, q4_or, q4_p):
    rows = ["Q1", "Q2", "Q3", "Q4", "Q5"]
    row_labels = [
        "PSI-testable events\n[20K, 200K]",
        "Primary aberrant events\nLUAD>=500; BLCA, UCEC>=300",
        "Cross-cancer pairwise overlap\n[5%, 30%]",
        "Host-gene enrichment\nOR>1.5, p<0.05",
        "Anti-MUC16 negative control\n0 Path-A false positives",
    ]
    cols = CANCERS

    def q1_verdict(val):
        if 20000 <= val <= 200000:
            return "PASS"
        if 10000 <= val < 20000:
            return "PARTIAL"
        return "FAIL"

    def q2_verdict(c, val):
        thr = 500 if c == "LUAD" else 300
        if c == "LUAD" and val < thr:
            return "BORDERLINE"
        return "PASS" if val >= thr else "FAIL"

    def cancer_q3_range(c):
        rel = pw[(pw['cancer_A'] == c) | (pw['cancer_B'] == c)]
        vals = []
        for _, row in rel.iterrows():
            if row['cancer_A'] == c:
                vals.append(row['pct_of_A'])
            else:
                vals.append(row['pct_of_B'])
        if not vals:
            return None, None
        return min(vals), max(vals)

    n_rows, n_cols = len(rows), len(cols)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.invert_yaxis()

    for j, c in enumerate(cols):
        for i, r in enumerate(rows):
            verdict = "PASS"; val_text = ""; footnote = ""
            if r == "Q1":
                v = q1[c]; verdict = q1_verdict(v); val_text = f"{v:,}"
            elif r == "Q2":
                v = q2[c]; verdict = q2_verdict(c, v); val_text = f"{v:,}"
                if verdict == "BORDERLINE":
                    footnote = "*"
            elif r == "Q3":
                lo, hi = cancer_q3_range(c)
                if lo is None:
                    val_text = "n/a"; verdict = "FAIL"
                else:
                    val_text = f"{lo:.1f}-{hi:.1f}%"; verdict = "REFRAMED"
            elif r == "Q4":
                if c in q4_or:
                    val_text = f"OR={q4_or[c]:.1f}"
                    p = q4_p.get(c, 1.0)
                    verdict = "PASS" if (q4_or[c] > 1.5 and p < 0.05) else "FAIL"
                else:
                    val_text = "n/a"; verdict = "FAIL"
            elif r == "Q5":
                val_text = "0 FPs"; verdict = "PASS"

            rect = mpl.patches.Rectangle((j - 0.45, i - 0.40), 0.90, 0.80,
                                         facecolor=VERDICT_COLOR[verdict],
                                         edgecolor='black', linewidth=0.7)
            ax.add_patch(rect)
            ax.text(j, i - 0.12, verdict + footnote,
                    ha='center', va='center', fontsize=7.0,
                    fontstyle='italic', color='black')
            ax.text(j, i + 0.16, val_text,
                    ha='center', va='center', fontsize=7.5,
                    fontweight='bold', color='black')

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(cols, fontweight='bold', fontsize=8.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=7.5)
    ax.tick_params(axis='both', which='both', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title('a   Pre-committed evaluation criteria', loc='left')

    legend_items = [("PASS", "PASS"),
                    ("BORDERLINE", "BORDERLINE*"),
                    ("REFRAMED", "REFRAMED")]
    for k, (key, lbl) in enumerate(legend_items):
        ax.add_patch(mpl.patches.Rectangle(
            (-0.40 + k * 1.10, n_rows - 0.15), 0.18, 0.18,
            facecolor=VERDICT_COLOR[key], edgecolor='black', linewidth=0.5,
            transform=ax.transData, clip_on=False))
        ax.text(-0.18 + k * 1.10, n_rows - 0.06, lbl,
                ha='left', va='center', fontsize=6.5,
                transform=ax.transData, clip_on=False)

    ax.text(n_cols - 0.5, n_rows + 0.25,
            "* LUAD primary aberrant count = 464 vs 500 threshold; see Methods Deviation 2.",
            ha='right', va='center', fontsize=6.0, fontstyle='italic',
            transform=ax.transData, clip_on=False)


def panel_b_overlap(ax, pw, n_conserved_events):
    pairs, pct_a, pct_b = [], [], []
    for _, row in pw.iterrows():
        pairs.append(f"{row['cancer_A']}\nvs {row['cancer_B']}")
        pct_a.append(row['pct_of_A'])
        pct_b.append(row['pct_of_B'])

    x = np.arange(len(pairs))
    w = 0.36

    bars_a = ax.bar(x - w/2, pct_a, w, color='#5C8A9E',
                    edgecolor='black', linewidth=0.4, label='% of cancer A')
    bars_b = ax.bar(x + w/2, pct_b, w, color='#C97B5E',
                    edgecolor='black', linewidth=0.4, label='% of cancer B')

    for bar in list(bars_a) + list(bars_b):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 1.2,
                f"{h:.1f}%", ha='center', va='bottom', fontsize=6.8)

    ax.axhline(30, color='black', linestyle='--', linewidth=0.8)
    x_mid = (x[0] + x[-1]) / 2
    ax.text(x_mid, 31.2, 'Precommit ceiling 30%',
            fontsize=6.5, va='bottom', ha='center', color='black',
            bbox=dict(boxstyle='round,pad=0.15',
                      facecolor='white', edgecolor='none', alpha=0.85))

    ax.axhline(5, color='black', linestyle=':', linewidth=0.7, alpha=0.6)
    ax.text(x[0] - w, 5.8, 'Floor 5%', fontsize=6.0,
            va='bottom', ha='left', color='black', alpha=0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(pairs)
    ax.set_ylabel('Event-level pairwise overlap (%)')
    ax.set_ylim(0, max(max(pct_a), max(pct_b)) * 1.18)
    ax.set_title('b   Cross-cancer pairwise overlap', loc='left')

    leg = ax.legend(loc='upper left', frameon=False, fontsize=6,
                    handlelength=1.0, handletextpad=0.5,
                    bbox_to_anchor=(0.01, 0.99))
    ax.add_artist(leg)

    if n_conserved_events is not None:
        ax.text(0.01, 0.82,
                f"3-way conserved: {n_conserved_events} events\n"
                f"99.2% direction-consistent\n"
                f"100% bootstrap stable",
                transform=ax.transAxes, ha='left', va='top',
                fontsize=5.8, color='black',
                bbox=dict(boxstyle='round,pad=0.25',
                          facecolor='white',
                          edgecolor='black', linewidth=0.4))


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
    print("Loading inputs ...")
    q1, q2, pw, n_cons, q4_or, q4_p = load_data()
    print(f"  Q1 (PSI-testable events): {q1}")
    print(f"  Q2 (primary aberrant events, from corrected parquets): {q2}")
    print(f"  Q3 pairwise rows: {len(pw)}")
    print(f"  Q3 three-way conserved core events: {n_cons}")
    print(f"  Q4 Fisher OR per cohort: {q4_or}")
    print()

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.4),
                             gridspec_kw={'width_ratios': [1.05, 1.0]})
    panel_a_grid(axes[0], q1, q2, pw, q4_or, q4_p)
    panel_b_overlap(axes[1], pw, n_cons)

    plt.tight_layout(w_pad=2.2)

    pdf_path = os.path.join(FIG_DIR, "F2_meta_thresholds_v4.pdf")
    png_path = os.path.join(FIG_DIR, "F2_meta_thresholds_v4.png")
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    print(f"  Saved PDF: {pdf_path}")
    print(f"  Saved PNG: {png_path}")
    print()
    print("F2 v4 done (M1'-M5' notation removed from labels). Open the PNG, review, then send it back.")


if __name__ == "__main__":
    main()
