# Reproduction manifest

Maps manuscript tables and figures to scripts and output files. Run all commands from the **repository root**.

| Manuscript asset | Command | Output path |
|------------------|---------|-------------|
| Main-text tables & figures | `make main` | `output/tables/`, `output/figures/` |
| Supplementary tables (Steiger, ROI × temperature, S5) | `make supplementary` | `output/supplementary/tables/` |
| Supplementary Figure S1 | `make supplementary` | `output/supplementary/figures/Supplementary_Figure_S1_weight_cadss_submitted.png` |
| Supplementary Figure S2 | `make supplementary` | `output/supplementary/figures/Supplementary_Figure_S2_roi_calib_temperature.png` |
| Supplementary Figure S3 | `make all` | `output/supplementary/figures/Supplementary_Figure_S3_connectivity.png` (from main connectivity plot) |
| Temp-adjusted models | `make main` | `output/tables/supp_table_tempadj_*`, `model_comparison_*` |

## Commands

```bash
R -e "renv::restore()"   # if using renv
make main                # main-text tables and figures
make supplementary       # all supplementary tables and figures (S3 needs main first)
make all                 # full reproduction
make verify              # spot-check key statistics
```

## Notes

SPM cluster inference (Supplementary Table S5) is external to this repository; the committed table is synced from `data/spm_ketamine_clusters_reference.csv` during `make supplementary`.
