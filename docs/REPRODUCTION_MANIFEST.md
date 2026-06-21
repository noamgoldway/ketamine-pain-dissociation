# Reproduction manifest

Maps manuscript tables and figures to scripts and output files. Run all commands from the **repository root**.

| Manuscript asset | Reproducible | Script / step | Output path |
|------------------|--------------|---------------|-------------|
| Main-text figures | Yes | `make main` | `output/figures/` |
| Supp Fig S1 (weight × post-bolus CADSS) | **Yes** | `make dose-equiv` | `output/revision/figures/Supplementary_Figure_S1_weight_cadss_submitted.*` |
| Supp Fig S2 (ROI × calib temp) | **Yes** | `make supp-fig-s2` | `output/revision/figures/Supplementary_Figure_S2_roi_calib_temperature.*` |
| Supp Fig S3 (connectivity) | Yes | `make main` + `make sync-s3` | `output/revision/figures/Supplementary_Figure_S3_connectivity.*` |
| Supp Tables S1–S14 | Yes | `make main` + `make supplementary` | `output/tables/`, `output/revision/tables/` |
| Steiger tests | Yes | `make supplementary` | `output/revision/tables/steiger_*.csv` |
| Temp-adjusted models | Yes | `make main` | `output/tables/supp_table_tempadj_*`, `model_comparison_*` |

## Commands

```bash
R -e "renv::restore()"   # if using renv
make all                 # main + supplementary + supplementary figures
make verify              # spot-check key statistics
```

## Notes

SPM cluster inference (Supplementary Table S5) is external to this repository; the committed table is synced from `data/spm_ketamine_clusters_reference.csv` via `make sync-s5`.
