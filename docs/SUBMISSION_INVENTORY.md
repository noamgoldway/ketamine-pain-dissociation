# NPP revision submission inventory

**Gate:** Only assets listed here are tracked in [ketamine-pain-dissociation](https://github.com/noamgoldway/ketamine-pain-dissociation). Everything else stays on Box under `revision/`.

**Source documents** (canonical copies in `text/npp_revision_2026/`):

| File | Role |
|------|------|
| `Main_revision_clean.docx` | Main manuscript (clean export) |
| `Main_revision.docx` | Main manuscript (full) |
| `Supplementary_Information_revision_clean.docx` | Supplementary Information |
| `rebuttal_letter.docx` | Response to reviewers |

Inventory generated from embedded media and in-text references in `revision/Final_files/` (May 2026).

## Main manuscript figures

| Asset | In-text refs | Embedded media (main docx) | Committed output |
|-------|--------------|--------------------------|------------------|
| Figure 1 | Study design / CADSS / task | `image1.png` | `output/revision/figures/Main_Figure_1.pdf` |
| Figure 2 | Behavioral results | `image2.png` | `output/revision/figures/Main_Figure_2.{pdf,png}` |
| Figure 3 | Univariate / brain–pain | `image3.png` | `output/revision/figures/Main_Figure_3.pdf` |
| Figure 4 | NPS / connectivity | `image4.png` | `output/revision/figures/Main_Figure_4.pdf` |

Main figures are **exported snapshots** aligned with the submitted PDF; regeneration uses `code/01_pain_ketamine_analysis_temp_covariate.r` → `output/figures/` then manual or scripted copy to `output/revision/figures/Main_Figure_*`.

## Supplementary figures

| Asset | Description | Reproducible | Committed output |
|-------|-------------|--------------|------------------|
| **Supplementary Figure S1** | Body weight vs **ketamine post-bolus** CADSS total (dose equivalence) | Yes | `output/revision/figures/Figure_S_dose_equivalence_weight_cadss.{png,pdf}`; legacy names `Figure_S_calibtemp_pain_slopes.*` |
| **Supplementary Figure S2** | ROI activation vs calibrated temperature (placebo vs ketamine) | **Frozen** | `output/revision/figures/Supplementary_Figure_S2_roi_calib_temperature.png` |
| **Supplementary Figure S3** | Within-network connectivity (ketamine vs placebo) | Yes | `output/figures/connectivity_all_networks.*` → `output/revision/figures/Supplementary_Figure_S3_connectivity.*` |

Supplementary docx embeds three images (`image1–3.png`); captions reference S1–S3.

## Supplementary tables

Tables **S1–S14** (and related) are embedded in `Supplementary_Information_revision_clean.docx` (17+ in-text “Table S” references). Statistics reproduced in git:

| Analysis | Script | Output |
|----------|--------|--------|
| Temperature-adjusted mixed models | `code/01_pain_ketamine_analysis_temp_covariate.r` | `output/tables/supp_table_tempadj_*`, `model_comparison_*calib*` |
| Steiger tests (pain–ROI vs dissociation) | `code/02_supplementary_revision.Rmd` | `output/revision/tables/steiger_*.csv` |
| ROI × calib_temp interaction / posthocs | `02_supplementary_revision.Rmd` | `output/revision/tables/roi_m*.csv`, `posthoc_m2_calibtemp_*.csv` |

## Main-text revision statistics (no separate figure)

| Claim | Script / data |
|-------|----------------|
| Demographics (weight 61.6 ± 9.4 kg) | `data/participants.csv` |
| NPS–pain Spearman ρ (.64 / .60) | `output/tables/nps_pain_corr_results.csv` |
| Pain `calib_temp` p = .070 | `output/tables/supp_table_tempadj_pain_model.csv` |
| ROI LRT χ²(1) = 9.11, F(1,681) = 7.16 | `output/tables/model_comparison_roi_calib_temp.csv`, `roi_m2_nice_anova.csv` |
| Steiger all ROI p > .11 | `output/revision/tables/steiger_roi_session_corr_diff.csv` |

## Rebuttal-only content

| Item | In rebuttal docx | Repo |
|------|------------------|------|
| Post-bolus CADSS descriptives (mean 24.94, SD 10.77, range 5–48) | Yes | `data/CADSS.csv` + `verify_manuscript_numbers.R` |
| Prior ketamine use counts | If cited | `data/prior_ketamine_use.csv` |

## Data files in git (Tier 1)

**Core (main paper):** `demog.csv`, `pain_ratings.csv`, `CADSS.csv`, `roi_beta_values_by_condition.csv`, `NPS.csv`, `pain_calibration.xlsx`, `within_connectivity.xlsx`

**Revision supplement:** `participants.csv`, `prior_ketamine_use.csv`, `CADSS_Weight_DoseEquivalence_Data.csv`

## Explicitly excluded from git

- `revision/manuscript/`, `revision/reviewer_materials/`, bulk `revision/tables/`, bulk `revision/figures/`
- `submission/NPP/`, `submission/Nature Communications/`, `submission/brain/`
- CADSS infusion-timing / arrow-task exploratory outputs (not in Final_files S1–S3)
- Raw DICOM, Qualtrics, fMRIPrep, first-level fMRI
- Internal markdown (`*_ADDITIONS.md`, `COMPLETION_STATUS.md`)
