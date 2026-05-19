# Reproduction manifest

Run all commands from the **repository root** (`submission/osf/` after clone).

| Submitted asset | Reproducible | Script / step | Output path |
|-----------------|-------------|---------------|-------------|
| Main Figs 1–4 | Partial | `01_pain_ketamine_analysis_temp_covariate.r` + manual export to revision PDFs | `output/figures/`, `output/revision/figures/Main_Figure_*` |
| Supp Fig S1 | Yes | `code/revision/plot_dose_equivalence.R` | `output/revision/figures/Figure_S_calibtemp_pain_slopes.*` |
| Supp Fig S2 | Frozen | Submitted docx embed | `output/revision/figures/Supplementary_Figure_S2_roi_calib_temperature.png` |
| Supp Fig S3 | Yes | Step 1 → connectivity plot | `output/figures/connectivity_all_networks.*` → `output/revision/figures/Supplementary_Figure_S3_connectivity.*` |
| Steiger tables | Yes | `Rscript -e 'rmarkdown::render("code/02_supplementary_revision.Rmd")'` | `output/revision/tables/steiger_*.csv` |
| Temp-adjusted supp tables | Yes | `01_pain_ketamine_analysis_temp_covariate.r` | `output/tables/supp_table_tempadj_*`, `model_comparison_*calib*` |
| CADSS timing scatter | No | Not in Final_files S1–S3 | — |

Reference layout: [interactive-avoidance-mental-health_public](https://github.com/tobywise/interactive-avoidance-mental-health_public).
