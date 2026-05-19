# Manuscript number verification

Run from repository root:

```bash
make verify
```

## Results summary (spot-checks)

| Reported in manuscript | Source in repo | Status |
|------------------------|----------------|--------|
| Supp Fig S1: Pearson **r = −0.056**, **p = 0.777** | `CADSS_Weight_DoseEquivalence_Data.csv` → `CADSS_total_post_avg` | **Matches** when using session-averaged post-infusion CADSS |
| Supp Fig S1 caption: “post-**bolus**” CADSS | `CADSS_total_post_session2` (ketamine `_2_post` in `CADSS.csv`) | **Does not match** caption stats (r ≈ 0.04, p ≈ 0.86) |
| Main text: weight **61.6 ± 9.4 kg** | `data/participants.csv` | **Matches** (mean/sd within rounding) |
| Steiger: all ROI **p > 0.11**; dlPFC **p = 0.1139** | `output/revision/tables/steiger_roi_session_corr_diff.csv`, recomputed from `data/` | **Matches** |
| Supp Table: pain `calib_temp` **p = 0.070** | `output/tables/supp_table_tempadj_pain_model.csv` | **Matches** |

## Action taken in code

`code/revision/plot_dose_equivalence.R` uses **`CADSS_total_post_avg`** so the reproduced correlation matches the **reported** r and p. The supplementary figure caption still says “post-bolus”; consider aligning caption with the average measure or switching the analysis to ketamine-only post-bolus and updating the reported statistics.

## Not fully re-run in CI

- Full `01_pain_ketamine_analysis_temp_covariate.r` (long runtime; committed `output/tables/` are the reference)
- `02_supplementary_revision.Rmd` Steiger block is verified independently; full Rmd needs pandoc to render
- Supplementary Figure S2: frozen embed (no script)

## Steiger recomputation

Independent R code using `data/roi_beta_values_by_condition.csv`, `pain_ratings.csv`, and `psych::r.test` reproduces committed `steiger_*.csv` exactly (including dlPFC p = 0.1139).
