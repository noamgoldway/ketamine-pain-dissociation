# The Analgesic and Dissociative Properties of Ketamine are Separate and Correspond to Distinct Neural Mechanisms

_Goldway et al._

This repository contains the analysis code and data needed to reproduce the statistical results in our ketamine–pain dissociation study ([Goldway et al., bioRxiv](https://www.biorxiv.org/content/10.1101/2025.07.25.666594v1)). It is meant for readers of the paper who want to rerun the models, inspect the tables, or build on the analysis pipeline.

## What is included

- **Behavioral and summary neuroimaging data** in `data/` (pain ratings, CADSS, NPS, ROI betas, connectivity matrices, demographics).
- **R scripts** that fit the mixed models, post-hoc contrasts, Steiger tests, and supplementary analyses reported in the main text and supplement.
- **Committed outputs** in `output/` (tables and figures corresponding to manuscript results).
- **Manuscript files** in `text/manuscript/` (main text and supplementary information Word documents, for reference).

fMRI preprocessing and first-level modeling are **not** part of this repository. ROI activations and connectivity were computed upstream; we provide the derived summaries used in the published analyses.

## Quick start

Clone the repository and run all analyses from the repository root:

```bash
git clone https://github.com/noamgoldway/ketamine-pain-dissociation.git
cd ketamine-pain-dissociation

R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
R -e "renv::restore()"

export ROOT_DIR="$(pwd)"
make all
```

`make all` runs the main analysis, supplementary analyses, and supplementary figures S1–S3. Results are written to `output/tables/`, `output/figures/`, and `output/revision/`.

To spot-check that key statistics match the manuscript tables:

```bash
make verify
```

## Repository layout

```
├── README.md
├── Makefile
├── renv.lock
├── data/                 # Analysis inputs
├── code/
│   ├── 01_pain_ketamine_analysis_temp_covariate.r   # Main models & figures
│   ├── 02_supplementary_revision.Rmd                # Steiger tests, ROI × temperature
│   └── revision/                                    # Supplementary figure scripts
├── output/
│   ├── tables/           # Main-text statistics
│   ├── figures/          # Main-text figures
│   └── revision/         # Supplementary tables & figures
├── text/manuscript/      # Main text & supplementary Word files
└── docs/                 # Additional documentation for full reproduction checks
```

## Reproducing the analyses

All commands below are run from the repository root.

### 1. Main analyses (pain, CADSS, NPS, ROI activation, connectivity)

```bash
make main
```

Equivalent: `Rscript code/01_pain_ketamine_analysis_temp_covariate.r`

### 2. Supplementary analyses (Steiger tests, calibration-temperature models)

```bash
make supplementary
```

### 3. Supplementary Figure S1 (body weight vs. post-infusion CADSS)

```bash
make dose-equiv
```

### 4. Run everything

```bash
make all
```

For a detailed mapping between scripts, outputs, and manuscript tables, see [`docs/REPRODUCTION_MANIFEST.md`](docs/REPRODUCTION_MANIFEST.md).

## Data files

| File | Description |
|------|-------------|
| `pain_ratings.csv` | Trial-level pain ratings |
| `CADSS.csv` | Dissociation (CADSS) scores by session and time point |
| `NPS.csv` | Neurophysiological pain signature scores |
| `roi_beta_values_by_condition.csv` | ROI activation betas by condition |
| `within_connectivity.xlsx` | Resting-state connectivity summaries |
| `pain_calibration.xlsx` | Calibration temperatures and pain slopes |
| `demog.csv` | Demographics |
| `participants.csv`, `prior_ketamine_use.csv` | Participant metadata |
| `CADSS_Weight_DoseEquivalence_Data.csv` | Data for Supplementary Figure S1 |

## Installing dependencies

### R (recommended: renv)

```bash
R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
R -e "renv::restore()"
```

If `renv` is unavailable, install packages manually (R ≥ 4.2):

```r
install.packages(c(
  "tidyverse", "readxl", "broom", "broom.mixed", "lme4", "emmeans",
  "afex", "Hmisc", "patchwork", "janitor", "corrplot", "psych", "knitr", "rmarkdown"
))
```

### Python

Python is **optional** and only needed if you want to run the full manuscript cross-check (`make verify-all`). The statistical analyses are entirely in R.

## Citation

Goldway, Noam, Talma Hendler, Itamar Jalon, Yotam Pasternak, Roy Sar-El, Dan Mirelman, Noam Sarna, Nili Green, Yara Agbaria, and Haggai Sharon. "The Analgesic and Dissociative Properties of Ketamine are Separate and Correspond to Distinct Neural Mechanisms." *bioRxiv* (2025): 2025-07.
