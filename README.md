# The Analgesic and Dissociative Properties of Ketamine are Separate and Correspond to Distinct Neural Mechanisms

_Goldway et al. — analysis code and derived data for the NPP revision resubmission._

Repository: [github.com/noamgoldway/ketamine-pain-dissociation](https://github.com/noamgoldway/ketamine-pain-dissociation)

Submitted manuscripts (snapshots): [`text/npp_revision_2026/`](text/npp_revision_2026/)

---

This repository reproduces **Tier 1** analyses from derived CSVs (behavioral ratings, ROI betas, CADSS, NPS, connectivity summaries). fMRI preprocessing is upstream of the files in `data/`.

Layout follows [tobywise/interactive-avoidance-mental-health_public](https://github.com/tobywise/interactive-avoidance-mental-health_public): clear structure, install instructions, and **numbered reproduction steps**.

## Repository structure

```
├── README.md
├── Makefile
├── renv.lock                 # R package versions (optional but recommended)
├── requirements.txt          # Python (optional; timing scripts not in submission scope)
├── docs/
│   ├── SUBMISSION_INVENTORY.md
│   ├── REPRODUCTION_MANIFEST.md
│   └── VERIFICATION.md
├── data/                     # Analysis inputs (derived summaries, in git)
├── code/
│   ├── 00_pain_ketamine_analysis_legacy.r
│   ├── 01_pain_ketamine_analysis_temp_covariate.r
│   ├── 02_supplementary_revision.Rmd
│   ├── verify_manuscript_numbers.R
│   └── revision/
│       └── plot_dose_equivalence.R
├── output/
│   ├── tables/
│   ├── figures/
│   └── revision/             # NPP resubmission figures & Steiger tables
└── text/
    ├── Main.pdf              # Prior OSF main PDF
    └── npp_revision_2026/  # Submitted Word files (May 2026)
```

## Installing dependencies

**Run all commands from this directory** (repository root after clone).

### R (recommended: renv)

```bash
cd submission/osf   # or repo root if you cloned ketamine-pain-dissociation directly
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

### Python (optional)

Only needed for exploratory CADSS-timing scripts **not** in the submission inventory:

```bash
pip install -r requirements.txt
```

## Reproducing the analyses

Numbered steps match `docs/REPRODUCTION_MANIFEST.md`.

### 1. Main behavioral and fMRI summary analyses

```bash
export ROOT_DIR="$(pwd)"
make main
# or: Rscript code/01_pain_ketamine_analysis_temp_covariate.r
```

Writes to `output/tables/` and `output/figures/`.

Legacy pipeline (original submission, without calibration-temperature extensions):

```bash
Rscript code/00_pain_ketamine_analysis_legacy.r
```

### 2. Supplementary revision analyses (Steiger, ROI × calibration temperature)

```bash
make supplementary
```

Runs `code/02_supplementary_revision.Rmd` via `knitr::purl` (no pandoc required). Writes to `output/revision/tables/`.

### 3. Supplementary Figure S1 (weight vs ketamine post-bolus CADSS)

```bash
make dose-equiv
```

Writes `output/revision/figures/Figure_S_dose_equivalence_weight_cadss.{png,pdf}`.

### All steps + verification

```bash
make all
make verify
make check
```

See `docs/VERIFICATION.md` for statistics checked against `text/npp_revision_2026/`.

## Data

Files in `data/` (see `docs/SUBMISSION_INVENTORY.md`):

- `demog.csv`, `pain_ratings.csv`, `CADSS.csv`, `roi_beta_values_by_condition.csv`, `NPS.csv`
- `pain_calibration.xlsx`, `within_connectivity.xlsx`
- `participants.csv`, `prior_ketamine_use.csv`, `CADSS_Weight_DoseEquivalence_Data.csv`

Unlike the Wise reference repo, Tier-1 inputs are **committed in git** (small derived CSVs) rather than downloaded via script.

## What is not rebuilt here

See `docs/SUBMISSION_INVENTORY.md`. Supplementary Figure S2 is a **frozen** export. CADSS infusion-timing and arrow-task figures are out of scope.

## Citation

If you use this code, please cite the paper (preprint/manuscript as available).
