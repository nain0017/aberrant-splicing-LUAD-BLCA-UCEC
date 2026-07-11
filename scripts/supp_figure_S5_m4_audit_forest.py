import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = "/content/drive/MyDrive/Paper5_SharedAntigens"
STAGE_C = os.path.join(ROOT, "data/processed/psi/Stage2C_corrected")
FIG = os.path.join(STAGE_C, "figures")
os.makedirs(FIG, exist_ok=True)

# Locked from cell 122 (M4_AUDIT_REPORT.txt)
BACKGROUND = 19782
ABERRANT = 940

TESTS = [
    ("Test 1: 71-gene curated",            "baseline",     71,  29, 14.25, 2.83e-20),
    ("Test 2: MSigDB EMT (n=200)",         "independent", 200,  22,  2.51, 2.22e-04),
    ("Test 3: MSigDB Myogenesis (n=200)",  "independent", 200,  27,  3.19, 9.86e-07),
    ("Test 6: de-circularized (n=54)",     "defensive",    54,  12,  5.79, 6.68e-06),
    ("Test 4a: KRAS Signaling Up",         "neg_control", 200,  10,  1.06, 4.81e-01),
    ("Test 4b: DNA Repair",                "neg_control", 150,   8,  1.13, 4.21e-01),
    ("Test 4c: Allograft Rejection",       "neg_control", 200,  15,  1.64, 5.47e-02),
    ("Test 4d: Pancreas Beta Cells",       "neg_control",  40,   1,  0.51, 8.58e-01),
    ("Test 4e: Heme Metabolism",           "neg_control", 200,  13,  1.40, 1.58e-01),
]

GROUP_COLOR = {
    "baseline":     "#2E5A88",
    "independent":  "#2E86AB",
    "defensive":    "#5F8A6C",
    "neg_control":  "#888888",
}
GROUP_LABEL = {
    "baseline":     "Baseline (curated reference)",
    "independent":  "Independent curation (MSigDB)",
    "defensive":    "De-circularized subset",
    "neg_control":  "Negative-control hallmarks",
}

STYLE = {
    'figure.dpi': 150, 'savefig.dpi': 300,
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 8.0,
    'axes.titlesize': 9.0, 'axes.titleweight': 'bold',
    'axes.labelsize': 8.0,
    'xtick.labelsize': 7.0, 'ytick.labelsize': 7.5,
    'axes.linewidth': 0.7,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
}
for k, v in STYLE.items():
    mpl.rcParams[k] = v

def or_ci(ref, ov, a_total=ABERRANT, bg=BACKGROUND):
    a = ov; b = a_total - ov; c = ref - ov; d = bg - a_total - c
    if any(x == 0 for x in [a, b, c, d]):
        a += 0.5; b += 0.5; c += 0.5; d += 0.5
    OR = (a * d) / (b * c)
    se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    lo = float(np.exp(np.log(OR) - 1.96 * se))
    hi = float(np.exp(np.log(OR) + 1.96 * se))
    return float(OR), lo, hi

n_rows = len(TESTS)
y_positions = np.arange(n_rows)[::-1]

# Compute all OR/CI values up front
results = []
for label, group, ref, ov, or_rep, p_rep in TESTS:
    OR, lo, hi = or_ci(ref, ov)
    results.append((label, group, OR, lo, hi, p_rep))

# Two-panel figure: forest on left, table on right
fig, (axL, axR) = plt.subplots(
    1, 2, figsize=(11.5, 5.4),
    gridspec_kw={'width_ratios': [2.0, 1.0], 'wspace': 0.05}
)

# === LEFT PANEL: forest plot ===
for y, (label, group, OR, lo, hi, p_rep) in zip(y_positions, results):
    color = GROUP_COLOR[group]
    axL.hlines(y, lo, hi, color=color, linewidth=1.8, alpha=0.85)
    axL.scatter([OR], [y], marker='D', s=52, color=color,
                edgecolor='black', linewidth=0.5, zorder=4)

axL.axvline(1.0, color='black', linestyle='--', linewidth=0.7, alpha=0.7)

# Group dividers
prev_group = None
for i, (_, group, *_) in enumerate(TESTS):
    if prev_group is not None and group != prev_group:
        y_div = (y_positions[i-1] + y_positions[i]) / 2.0
        axL.axhline(y_div, color='lightgrey', linewidth=0.5, alpha=0.7, zorder=0)
    prev_group = group

axL.set_xscale('log')
axL.set_xlim(0.3, 30)
axL.set_xticks([0.5, 1, 2, 5, 10, 20])
axL.set_xticklabels(['0.5', '1', '2', '5', '10', '20'])
axL.set_xlabel('Odds ratio (log scale; vertical dashed = no enrichment)')

axL.set_yticks(y_positions)
axL.set_yticklabels([t[0] for t in TESTS])
axL.set_ylim(-0.7, n_rows - 0.3)
axL.spines['top'].set_visible(False)
axL.spines['right'].set_visible(False)
axL.set_title("M4 host-gene enrichment 6-test audit", loc='left', pad=10)

# === RIGHT PANEL: text table aligned to forest rows ===
axR.set_xlim(0, 1)
axR.set_ylim(-0.7, n_rows - 0.3)  # matches left panel y-range exactly
axR.axis('off')

# Column headers (top of panel, in data y-coordinates)
header_y = n_rows - 0.45
axR.text(0.04, header_y, "OR  (95% CI)",
         ha='left', va='center',
         fontsize=8.0, fontweight='bold')
axR.text(0.62, header_y, "p-value",
         ha='left', va='center',
         fontsize=8.0, fontweight='bold')

# Header underline
axR.plot([0.02, 0.95], [header_y - 0.35, header_y - 0.35],
         color='black', linewidth=0.6, transform=axR.transData)

# Row text
for y, (label, group, OR, lo, hi, p_rep) in zip(y_positions, results):
    or_str = f"{OR:>5.2f}  ({lo:.2f}–{hi:.2f})"
    if p_rep < 1e-3:
        p_str = f"{p_rep:.1e}"
    elif p_rep < 1e-2:
        p_str = f"{p_rep:.3f}"
    else:
        p_str = f"{p_rep:.2f}"
    axR.text(0.04, y, or_str, ha='left', va='center',
             fontsize=7.5, family='DejaVu Sans Mono')
    axR.text(0.62, y, p_str, ha='left', va='center',
             fontsize=7.5, family='DejaVu Sans Mono')

# Mirror the group divider rules on the right panel
prev_group = None
for i, (_, group, *_) in enumerate(TESTS):
    if prev_group is not None and group != prev_group:
        y_div = (y_positions[i-1] + y_positions[i]) / 2.0
        axR.axhline(y_div, color='lightgrey', linewidth=0.5, alpha=0.7,
                    xmin=0.02, xmax=0.95, zorder=0)
    prev_group = group

# === LEGEND below the forest panel (in figure coords) ===
from matplotlib.patches import Patch
handles = [Patch(facecolor=GROUP_COLOR[g], edgecolor='black',
                 linewidth=0.4, label=GROUP_LABEL[g])
           for g in ['baseline', 'independent', 'defensive', 'neg_control']]
fig.legend(handles=handles,
           loc='lower left', bbox_to_anchor=(0.08, -0.02),
           ncol=4, frameon=True, fontsize=7.5,
           framealpha=0.95, edgecolor='black')

# === FOOTER: Test 5 randomization ===
fig.text(0.08, -0.08,
         "Test 5 (randomization): 0 of 1,000 random 71-gene reference sets "
         "reached the observed OR = 14.25; "
         "null distribution max = 3.31; empirical permutation p = 0.001.",
         fontsize=7.5, fontstyle='italic')

plt.tight_layout()

pdf_path = os.path.join(FIG, "S8_m4_audit_forest_v2.pdf")
png_path = os.path.join(FIG, "S8_m4_audit_forest_v2.png")
fig.savefig(pdf_path, bbox_inches='tight')
fig.savefig(png_path, bbox_inches='tight', dpi=300)
plt.show()

print(f"Saved PDF: {pdf_path}")
print(f"Saved PNG: {png_path}")
print("S8 v2 done.")