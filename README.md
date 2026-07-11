# aberrant-splicing-LUAD-BLCA-UCEC

Code and analysis notebook accompanying:

**"A dispersion-aware, expression-gated pipeline reveals a conserved
cassette-exon program across three solid tumors"**
Md. Zulkarnain Sajid
Submitted to NAR Genomics and Bioinformatics.

ORCID: 0009-0007-9421-3016

## What this is

A four-stage pipeline that combines TPM-based gene expression gating,
SUPPA2 cassette-exon cataloging, per-sample PSI computation with strict
coverage requirements, and a dispersion-aware beta-binomial likelihood-ratio
test to identify tumor-aberrant cassette-exon splicing events in LUAD,
BLCA, and UCEC (TCGA), using GTEx and TCGA solid-tissue-normal cohorts as
controls. All RNA-seq inputs are recount3-derived, uniformly re-aligned by
the Monorail pipeline; no re-alignment was performed by this project.

## Repository contents

```
paper_5__4_.ipynb              Locked analysis notebook (152 cells).
                                All headline numbers in the manuscript
                                trace to outputs printed in this notebook.
                                This notebook is the full development
                                history, including earlier superseded
                                analysis attempts; the scripts/ folder
                                below contains only the final versions
                                that produced the published figures.

STAGE_C_DEVIATIONS.md          Deviation register (D1-D7). Every
                                methodological deviation from the original
                                pre-commitment, logged before the resulting
                                numbers were used in the manuscript.

M4_AUDIT_REPORT.txt            Six-test audit of the host-gene enrichment
                                (M4) result: circularity check, independent
                                curation cross-check, five negative
                                controls, randomization control, and a
                                de-circularized reference subset.

SupplementaryTables.xlsx       Supplementary Tables T1-T10 (event
                                catalogs, conserved-core annotations,
                                enrichment audit, pathway enrichment,
                                hallmark gene coverage, sample manifest,
                                software versions).

scripts/                       Standalone figure- and table-generation
                                scripts, one per published figure/table.
                                Each can be pasted into a single Colab
                                cell and re-run against the
                                Stage2C_corrected/ outputs on Drive to
                                regenerate the corresponding figure or
                                table exactly. Only the final version
                                used in the published manuscript is
                                included; earlier iterations are part of
                                the notebook's history above, not
                                duplicated here.

  figure1_pipeline_overview.py       -> Figure 1
  figure2_meta_thresholds.py         -> Figure 2
  figure3_volcanoes.py               -> Figure 3
  figure4_overlap_hallmarks.py       -> Figure 4
  figure5_hallmark_case_studies.py   -> Figure 5
  supp_figure_S2_calibration_qq.py   -> Supplementary Figure S2
  supp_figure_S4_snapping_diagnostic.py -> Supplementary Figure S4
  supp_figure_S8_m4_audit_forest.py  -> Supplementary Figure S8
  generate_supplementary_tables.py   -> Supplementary Tables T1-T10

requirements.txt                Pinned Python package versions.

LICENSE                         MIT.
```

## Data availability

All RNA-seq data used are publicly available through the recount3 resource
(http://rna.recount.bio), sourced from TCGA and GTEx. This repository does
not re-host raw sequencing data. Intermediate processed outputs (corrected
per-cancer result parquets, PSI matrices, enrichment audit CSVs) referenced
by the scripts above live in the author's Google Drive project folder
(`Paper5_SharedAntigens/data/processed/psi/Stage2C_corrected/`) and are
available from the corresponding author on reasonable request; a full data
snapshot may be added to the Zenodo archive of this repository if journal
or reviewer requirements call for it.

## Reproducibility statement

Every headline number in the manuscript (primary/secondary aberrant event
counts, conserved-core statistics, host-gene enrichment odds ratios,
calibration lambdas, hallmark gene coverage) traces to a printed output
cell in `paper_5__4_.ipynb`. No figure in this project is generated from a
hardcoded results dictionary or a synthetic distribution parameterized by
already-known summary statistics; all figures read from saved parquet/CSV
outputs of the upstream analysis cells. The deviation register
(`STAGE_C_DEVIATIONS.md`) documents every methodological change made after
the original pre-commitment, each dated and logged before use.

## Citation

A citable Zenodo DOI for this repository will be added upon archival.
