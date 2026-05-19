# NPP revision submission inventory

Source documents (canonical copies in `text/npp_revision_2026/`):

- `Main_revision_clean.docx` (and `Main_revision.docx`)
- `Supplementary_Information_revision_clean.docx`
- `rebuttal_letter.docx`

## Main manuscript figures

| Asset | Description | Committed output |
|-------|-------------|------------------|
| Figure 1 | Study design / CADSS timing schematic | `output/revision/figures/Main_Figure_1.pdf` |
| Figure 2 | Behavioral results | `output/revision/figures/Main_Figure_2.pdf`, `.png` |
| Figure 3 | Univariate / brain–pain | `output/revision/figures/Main_Figure_3.pdf` |
| Figure 4 | NPS / connectivity summary | `output/revision/figures/Main_Figure_4.pdf` |
| Table 1 | Demographics / design (in Word) | — (manuscript only) |

## Supplementary figures (revision)

| Asset | Description | Committed output |
|-------|-------------|------------------|
| Supplementary Figure S1 | Body weight vs post-infusion CADSS (session average) | `output/revision/figures/Figure_S_dose_equivalence_weight_cadss.png` (script); legacy `Figure_S_calibtemp_pain_slopes.png` |
| Supplementary Figure S2 | ROI activation vs calibrated temperature (placebo vs ketamine) | `output/revision/figures/Supplementary_Figure_S2_roi_calib_temperature.png` (frozen from submitted docx) |
| Supplementary Figure S3 | Within-network connectivity | `output/revision/figures/Supplementary_Figure_S3_connectivity.png` (from `output/figures/connectivity_all_networks.png`) |

## Supplementary tables

Tables S1–S14+ are embedded in `Supplementary_Information_revision_clean.docx`. Reproducible statistics for temperature-adjusted models and Steiger tests are in `output/tables/` and `output/revision/tables/` (see `REPRODUCTION_MANIFEST.md`).

## Analyses cited in revision text (not separate figures)

- Steiger tests (brain–pain correlations, placebo vs ketamine): `output/revision/tables/steiger_*.csv` via `code/02_supplementary_revision.Rmd`
- Calibration-temperature mixed models: `code/01_pain_ketamine_analysis_temp_covariate.r` → `output/tables/*temp*`, `supp_table_tempadj_*`

## Explicitly out of repo scope

- Exploratory revision workspace (`revision/tables/`, timing pipelines, arrow-task figures) unless added to this inventory later
- Raw DICOM / Qualtrics rebuild from `tlv_computer/`
