# Reproduction manifest

Run all commands from the repository root.

## Commands

| Command | Scripts | Output |
|---------|---------|--------|
| `make main` | `code/01_pain_ketamine_analysis_temp_covariate.r` | `output/tables/`, `output/figures/` |
| `make supplementary` | `code/02_supplementary.Rmd`, `code/supplementary/*.R` | `output/supplementary/` |
| `make all` | Both of the above | All outputs |
| `make verify` | `code/verify_manuscript_numbers.R` | Terminal spot-checks |

```bash
R -e "renv::restore()"   # if using renv
make all
make verify
```

## Output layout

**`output/tables/`** — mixed models, post-hoc contrasts, and model-comparison tables from the main script. This includes most supplementary statistics cited in the manuscript.

**`output/figures/`** — main-text figures from the main script.

**`output/supplementary/tables/`** — Steiger tests, ROI × calibration-temperature post-hocs, and the SPM cluster table (`supp_table_s5_ketamine_clusters.csv`, copied from `data/spm_ketamine_clusters_reference.csv`).

**`output/supplementary/figures/`** — supplementary figures. The connectivity figure is copied here from the main script.
