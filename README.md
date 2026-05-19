# The Analgesic and Dissociative Properties of Ketamine are Separate and Correspond to Distinct Neural Mechanisms

_Goldway et al. — analysis code and derived data for the NPP revision resubmission._

Repository: [github.com/noamgoldway/ketamine-pain-dissociation](https://github.com/noamgoldway/ketamine-pain-dissociation)

Submitted manuscripts (snapshots): `text/npp_revision_2026/`

---

This repository reproduces **Tier 1** analyses from derived CSVs (behavioral ratings, ROI betas, CADSS, NPS, connectivity summaries). fMRI preprocessing is upstream of the files in `data/`.

Layout follows [tobywise/interactive-avoidance-mental-health_public](https://github.com/tobywise/interactive-avoidance-mental-health_public): clear structure, install instructions, and **numbered reproduction steps**.

## Repository structure

```
├── README.md
├── Makefile
├── requirements.txt          # Python (revision figures only)
├── docs/
│   ├── SUBMISSION_INVENTORY.md
│   └── REPRODUCTION_MANIFEST.md
├── data/                     # Analysis inputs (derived summaries)
├── code/
│   ├── 00_pain_ketamine_analysis_legacy.r   # Original OSF pipeline
│   ├── 01_pain_ketamine_analysis_temp_covariate.r
│   ├── 02_supplementary_revision.Rmd
│   └── revision/
│       └── plot_dose_equivalence.R
├── output/
│   ├── tables/               # Main + supplementary statistics
│   ├── figures/              # Main paper figures
│   └── revision/             # NPP resubmission figures & Steiger tables
└── text/
    ├── Main.pdf
    └── npp_revision_2026/    # Submitted Word files
```

## Installing dependencies

### R

Use R ≥ 4.2. Install packages used by the scripts (or use `renv` if `renv.lock` is present):

```r
install.packages(c(
  "tidyverse", "readxl", "broom", "broom.mixed", "lme4", "emmeans",
  "afex", "Hmisc", "patchwork", "janitor", "corrplot", "psych", "rmarkdown"
))
```

### Python (optional; Supplementary Figure S1 only)

```bash
pip install -r requirements.txt
```

## Reproducing the analyses

**Run all commands from this directory** (repository root).

### 1. Main behavioral and fMRI summary analyses

```bash
export ROOT_DIR="$(pwd)"
Rscript code/01_pain_ketamine_analysis_temp_covariate.r
```

Writes to `output/tables/` and `output/figures/`.

Legacy pipeline (original submission, without calibration-temperature extensions):

```bash
Rscript code/00_pain_ketamine_analysis_legacy.r
```

### 2. Supplementary revision analyses (Steiger tests, ROI × calibration models)

```bash
Rscript -e 'rmarkdown::render("code/02_supplementary_revision.Rmd", quiet = TRUE)'
```

Writes to `output/revision/tables/` (and related `output/tables/` from the Rmd).

### 3. Supplementary Figure S1 (weight vs post-bolus CADSS)

```bash
Rscript code/revision/plot_dose_equivalence.R
```

Writes `output/revision/figures/Figure_S_dose_equivalence_weight_cadss.png` (and `.pdf`; legacy copies as `Figure_S_calibtemp_pain_slopes.*`).

### All steps

```bash
make all
make verify   # spot-check key stats vs submitted manuscript
make check
```

See `docs/VERIFICATION.md` for what is validated and one caption–data note for Supplementary Figure S1.

## Data

The following derived files are included in `data/`:

- `demog.csv`, `pain_ratings.csv`, `CADSS.csv`, `roi_beta_values_by_condition.csv`, `NPS.csv`
- `pain_calibration.xlsx`, `within_connectivity.xlsx`
- `participants.csv`, `prior_ketamine_use.csv`, `CADSS_Weight_DoseEquivalence_Data.csv`

## What is not rebuilt here

See `docs/SUBMISSION_INVENTORY.md`. Supplementary Figure S2 is committed as a **frozen** export matching the submitted supplementary document. CADSS infusion-timing exploratory analyses are not part of the Final_files figure list.

## Citation

If you use this code, please cite the paper (preprint/manuscript as available).
