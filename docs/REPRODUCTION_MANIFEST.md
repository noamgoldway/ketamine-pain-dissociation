# Reproduction manifest

Maps manuscript tables and figures to scripts and output files. Run all commands from the **repository root**.

| Manuscript asset | Reproducible | Script / step | Output path |
|------------------|--------------|---------------|-------------|
| Main Fig 1 | Partial (export) | `make main` | `output/revision/figures/Main_Figure_1.pdf` |
| Main Fig 2 | Partial | `make main` | `output/revision/figures/Main_Figure_2.{pdf,png}` |
| Main Fig 3 | Partial | `make main` | `output/revision/figures/Main_Figure_3.pdf` |
| Main Fig 4 | Partial | `make main` | `output/revision/figures/Main_Figure_4.pdf` |
| Supp Fig S1 (weight × post-bolus CADSS) | **Yes** | `make dose-equiv` | `output/revision/figures/Figure_S_dose_equivalence_weight_cadss.*` |
| Supp Fig S2 (ROI × calib temp) | **Yes** | `make supp-fig-s2` | `output/revision/figures/Supplementary_Figure_S2_roi_calib_temperature.*` |
| Supp Fig S3 (connectivity) | Yes | `make main` | `output/figures/connectivity_all_networks.*` → `Supplementary_Figure_S3_connectivity.*` |
| Supp Tables S1–S14 | Mixed | Main + supplementary scripts | `output/tables/`, `output/revision/tables/` |
| Steiger tests | Yes | `make supplementary` | `output/revision/tables/steiger_*.csv` |
| Temp-adjusted models | Yes | `make main` | `output/tables/supp_table_tempadj_*`, `model_comparison_*` |

## Commands (ordered)

```bash
R -e "renv::restore()"   # if using renv
make all                 # main + supplementary + supplementary figures
make verify              # spot-check key statistics
```

## Notes

| Asset | Note |
|-------|------|
| Main Figs 1–4 in `output/revision/figures/` | Committed PDFs/PNGs; full re-export from analysis figures may differ in styling |
| SPM cluster table (Supp Table S5) | First-level/cluster inference is external; table synced from `data/spm_ketamine_clusters_reference.csv` via `make sync-s5` |

For a full cross-check of manuscript numbers against Word tables, see [`docs/VERIFICATION.md`](VERIFICATION.md).
