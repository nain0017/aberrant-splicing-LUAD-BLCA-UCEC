import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm

ROOT = "/content/drive/MyDrive/Paper5_SharedAntigens"
STAGE_C = os.path.join(ROOT, "data/processed/psi/Stage2C_corrected")
FIG = os.path.join(STAGE_C, "figures")
os.makedirs(FIG, exist_ok=True)

CANCERS = ["LUAD", "BLCA", "UCEC"]
COHORT_COLOR = {"LUAD": "#2E86AB", "BLCA": "#E07A5F", "UCEC": "#76A646"}

STYLE = {
    'figure.dpi': 150, 'savefig.dpi': 300,
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 8.0, 'axes.titlesize': 9.0, 'axes.titleweight': 'bold',
    'axes.labelsize': 8.0, 'xtick.labelsize': 7.0, 'ytick.labelsize': 7.0,
    'axes.linewidth': 0.7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
}
for k, v in STYLE.items():
    mpl.rcParams[k] = v

print(f"Output dir: {FIG}")
print("Loading coverage-pass events per cancer ...")

data = {}
for c in CANCERS:
    df = pd.read_parquet(os.path.join(STAGE_C, f"{c}_corrected_results.parquet"))
    cov = df['coverage_pass'].fillna(False).astype(bool)
    df = df[cov].copy()
    df['n_min_well_covered'] = np.minimum(df['n_tumor_well_covered'],
                                          df['n_normal_well_covered'])
    data[c] = df
    diff = (df['delta_psi'] - df['delta_psi_median']).abs()
    p95 = float(diff.quantile(0.95))
    n_div = int((diff > 0.10).sum())
    max_med = float(df['delta_psi_median'].abs().max())
    print(f"  [{c}] coverage-pass={len(df):,}  "
          f"|mean-median| p95={p95:.3f}  "
          f"diverge>0.10={n_div:,}  "
          f"max|median|={max_med:.3f}")

fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.4), sharey=True)
scs = []

for ax, c in zip(axes, CANCERS):
    df = data[c]
    color = COHORT_COLOR[c]
    n = len(df)

    cvals = df['n_min_well_covered'].values.astype(float)
    cvals = np.clip(cvals, 1, cvals.max() if cvals.max() > 1 else 2)

    sc = ax.scatter(df['delta_psi'], df['delta_psi_median'],
                    s=4, c=cvals, cmap='viridis_r',
                    norm=LogNorm(vmin=max(cvals.min(), 1.0), vmax=cvals.max()),
                    alpha=0.45, edgecolors='none', rasterized=True)
    scs.append(sc)

    ax.plot([-1, 1], [-1, 1], color='black', linewidth=0.6,
            linestyle='--', alpha=0.6)
    ax.axhline(1.0, color='red', linewidth=0.6, linestyle=':', alpha=0.7)
    ax.axhline(-1.0, color='red', linewidth=0.6, linestyle=':', alpha=0.7)
    ax.axhline(0.0, color='black', linewidth=0.4, alpha=0.3)
    ax.axvline(0.0, color='black', linewidth=0.4, alpha=0.3)

    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.02, 1.02)
    ax.set_xlabel(r'mean $\Delta$PSI (corrected pipeline)')
    if ax.get_subplotspec().is_first_col():
        ax.set_ylabel(r'median $\Delta$PSI (old convention)')
    ax.set_title(f"{c}  (n={n:,})", loc='center', color=color)

    diff = (df['delta_psi'] - df['delta_psi_median']).abs()
    n_diverge = int((diff > 0.10).sum())
    pct_diverge = 100.0 * n_diverge / n if n > 0 else 0.0
    p95 = float(diff.quantile(0.95))
    max_med = float(df['delta_psi_median'].abs().max())

    inset = (f"|mean-median| p95: {p95:.3f}\n"
             f"diverge >0.10: {n_diverge:,} ({pct_diverge:.1f}%)\n"
             f"max |median|: {max_med:.3f}\n"
             f"(no snapping at +-1.0)")
    ax.text(0.04, 0.96, inset, transform=ax.transAxes,
            ha='left', va='top', fontsize=6.5,
            bbox=dict(boxstyle='round,pad=0.3',
                      facecolor='white', edgecolor='black', linewidth=0.4))

cbar = fig.colorbar(scs[-1], ax=axes, location='right',
                    shrink=0.85, aspect=18, pad=0.025)
cbar.set_label('min(well-covered T, well-covered N)\nper event (log scale)',
               fontsize=7.0)
cbar.ax.tick_params(labelsize=6.5)

pdf_path = os.path.join(FIG, "S4_snapping_diagnostic_v2.pdf")
png_path = os.path.join(FIG, "S4_snapping_diagnostic_v2.png")
fig.savefig(pdf_path, bbox_inches='tight')
fig.savefig(png_path, bbox_inches='tight', dpi=300)
plt.show()

print(f"\nSaved PDF: {pdf_path}")
print(f"Saved PNG: {png_path}")
print("S4 v2 done.")