# The Analgesic and Dissociative Properties of Ketamine are Separate and Correspond to Distinct Neural Mechanisms

_Goldway et al._

Analysis code and data to reproduce the statistical results in our ketamine–pain dissociation study ([Goldway et al., bioRxiv](https://www.biorxiv.org/content/10.1101/2025.07.25.666594v1)).

## What is included

- **Data** in `data/` — pain ratings, CADSS, NPS, ROI betas, connectivity summaries, demographics
- **R scripts** in `code/` — mixed models, post-hoc contrasts, Steiger tests, supplementary analyses
- **Outputs** in `output/` — tables and figures corresponding to manuscript results
- **Manuscript** in `text/manuscript/` — main text and supplementary information (for reference)

## Quick start

```bash
git clone https://github.com/noamgoldway/ketamine-pain-dissociation.git
cd ketamine-pain-dissociation

R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
R -e "renv::restore()"

export ROOT_DIR="$(pwd)"
make all
make verify
```

`make all` runs the main analysis, supplementary analyses, and supplementary figures S1–S3.

## Repository layout

```
├── README.md
├── Makefile
├── renv.lock
├── data/
├── code/
│   ├── 01_pain_ketamine_analysis_temp_covariate.r
│   ├── 02_supplementary_revision.Rmd
│   └── revision/
├── output/
│   ├── tables/
│   ├── figures/
│   └── revision/
├── text/manuscript/
└── docs/REPRODUCTION_MANIFEST.md
```

## Reproducing the analyses

| Step | Command | Output |
|------|---------|--------|
| Main analyses | `make main` | `output/tables/`, `output/figures/` |
| Supplementary analyses | `make supplementary` | `output/revision/tables/` |
| Supplementary Figure S1 | `make dose-equiv` | `output/revision/figures/` |
| Everything | `make all` | All of the above |

See [`docs/REPRODUCTION_MANIFEST.md`](docs/REPRODUCTION_MANIFEST.md) for a full mapping of scripts to manuscript tables and figures.

## Data files

| File | Description |
|------|-------------|
| `pain_ratings.csv` | Trial-level pain ratings |
| `CADSS.csv` | Dissociation (CADSS) scores by session and time point |
| `NPS.csv` | Neurological Pain Signature (NPS) scores |
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

## Citation

Goldway, Noam, Talma Hendler, Itamar Jalon, Yotam Pasternak, Roy Sar-El, Dan Mirelman, Noam Sarna, Nili Green, Yara Agbaria, and Haggai Sharon. "The Analgesic and Dissociative Properties of Ketamine are Separate and Correspond to Distinct Neural Mechanisms." *bioRxiv* (2025): 2025-07.
