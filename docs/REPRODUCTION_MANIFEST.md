# Reproduction manifest

Derived from `docs/SUBMISSION_INVENTORY.md`. **One row per submitted asset** (or analysis block cited in Final_files).

Run all commands from the **repository root** (`submission/osf/` after clone).

| Submitted asset | In git | Reproducible | Script / step | Output path |
|-----------------|--------|--------------|---------------|-------------|
| Main Fig 1 | Yes | Partial (export) | `make main` → copy/export | `output/revision/figures/Main_Figure_1.pdf` |
| Main Fig 2 | Yes | Partial | `make main` | `output/revision/figures/Main_Figure_2.{pdf,png}` |
| Main Fig 3 | Yes | Partial | `make main` | `output/revision/figures/Main_Figure_3.pdf` |
| Main Fig 4 | Yes | Partial | `make main` | `output/revision/figures/Main_Figure_4.pdf` |
| Supp Fig S1 (weight × post-bolus CADSS) | Yes | **Yes** | `make dose-equiv` | `output/revision/figures/Figure_S_dose_equivalence_weight_cadss.*` |
| Supp Fig S2 (ROI × calib temp) | Yes | **Frozen** | Docx embed only | `output/revision/figures/Supplementary_Figure_S2_roi_calib_temperature.png` |
| Supp Fig S3 (connectivity) | Yes | Yes | `make main` | `output/figures/connectivity_all_networks.*` → `Supplementary_Figure_S3_connectivity.*` |
| Supp Tables S1–S14 (embedded) | Partial | Mixed | Word tables; stats from steps below | `output/tables/`, `output/revision/tables/` |
| Steiger tests | Yes | Yes | `make supplementary` | `output/revision/tables/steiger_*.csv` |
| Temp-adjusted models | Yes | Yes | `make main` | `output/tables/supp_table_tempadj_*`, `model_comparison_*` |
| Manuscript snapshots | Yes | N/A | Copied from `revision/Final_files/` | `text/npp_revision_2026/*.docx` |

## Commands (ordered)

```bash
cd submission/osf   # repository root
R -e "renv::restore()"   # if using renv
make all      # main + supplementary + dose-equiv
make verify   # spot-check vs Final_files statistics
make check    # required output files exist
```

## Frozen / no script

| Asset | Note |
|-------|------|
| Supp Fig S2 | No rebuild script; PNG matches submitted supplementary docx |
| Main Figs 1–4 in revision folder | Committed PDFs/PNGs; full re-export from analysis figures may differ in styling |

Reference layout: [interactive-avoidance-mental-health_public](https://github.com/tobywise/interactive-avoidance-mental-health_public).
