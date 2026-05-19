# Manuscript number verification

Run from repository root (`submission/osf`):

```bash
make verify
```

Regenerate Supplementary Figure S1 (weight vs ketamine post-bolus CADSS):

```bash
make dose-equiv
```

## Last comprehensive re-run (2026-05-19)

| Step | Command | Result |
|------|---------|--------|
| Main pipeline | `make main` | Completed (~9 s); key table MD5s **unchanged** vs pre-run |
| Supplementary | `make supplementary` (purl + source; pandoc not required) | Completed; Steiger + ROI calib-temp tables refreshed |
| S1 figure | `make dose-equiv` | Post-bolus plot; r = 0.036, p = 0.855, R² = 0.001 |
| Verification | `make verify` + `make check` | **19/19 OK** |
| Docx scan | `revision/Final_files/*.docx` | S1 r/p in supplementary; rebuttal CADSS + χ²/F; main NPS + ROI stats present |

## Audit table (Final_files docx vs repo)

| Claim | Manuscript value | Reproduced value | Status | Script / source |
|-------|------------------|------------------|--------|-----------------|
| Supp Fig S1 / dose equivalence | Pearson **r = 0.04**, **p = 0.86** (n = 28; post-bolus ketamine) | r ≈ 0.036, p ≈ 0.855, R² ≈ 0.001 (`CADSS_total_post_session2`) | **OK** | `code/revision/plot_dose_equivalence.R`; `data/CADSS_Weight_DoseEquivalence_Data.csv` |
| Supp methods: post-infusion CADSS vs weight | Qualitative (no association) | Same as above | **OK** | `code/verify_manuscript_numbers.R` |
| Demographics: body weight | **61.6 ± 9.4 kg** | mean 61.6, SD 9.4 (n = 36) | **OK** | `data/participants.csv` |
| Rebuttal: post-bolus CADSS (ketamine) | mean **24.94**, SD **10.77**, range **5–48** | 24.94, 10.77, 5–48 (n = 33) | **OK** | `data/CADSS.csv` (sum of `_2_post` subscales) |
| Main: NPS vs pain (high intensity) | Spearman **ρ = .64** (placebo), **ρ = .60** (ketamine) | 0.643, 0.604 | **OK** | `output/tables/nps_pain_corr_results.csv` ← `01_pain_ketamine_analysis_temp_covariate.r` |
| Supp Table: pain `calib_temp` | **p = 0.070** | 0.070 | **OK** | `output/tables/supp_table_tempadj_pain_model.csv` |
| Rebuttal: pain + `calib_temp` LRT | **χ²(1) = 2.66**, **p = .103** | 2.66, 0.103 | **OK** | `output/tables/model_comparison_pain_calib_temp.csv` |
| Main / rebuttal: ROI + `calib_temp` LRT | **χ²(1) = 9.11**, **p = .003** | 9.11, 0.0025 (rounded .003) | **OK** | `output/tables/model_comparison_roi_calib_temp.csv` |
| Main: intensity × session × `calib_temp` (ROI) | **F(1, 681) = 7.16**, **p = .008** | F = 7.16, p = .008 | **OK** | `output/tables/roi_m2_nice_anova.csv` |
| Main: Steiger (session-specific pain–ROI vs dissociation) | all ROI **p > 0.11**; dlPFC **p = 0.1139** | min p > 0.11; dlPFC 0.1139 | **OK** | `output/revision/tables/steiger_roi_session_corr_diff.csv`; `02_supplementary_revision.Rmd` |
| Steiger (independent recomputation) | As committed CSV | Exact match from `data/` | **OK** | `psych::r.test` on `roi_beta_values_by_condition.csv` + `pain_ratings.csv` |
| Supplementary Figure S2 | ROI activation vs calibrated temperature | Frozen PNG in `text/npp_revision_2026/` | **Frozen** | No rebuild script in repo |
| Supplementary Figure S3 | Connectivity | `output/figures/connectivity_all_networks.*` | **Committed** | `01_pain_ketamine_analysis_temp_covariate.r` |
| Arrow / timing supplementary figures | Various | Not in OSF package | **Out of scope** | Frozen or not submitted |
| Full main pipeline (~30 min) | All `output/tables/` | Spot-checked tables above | **Spot-check** | `make main` |

## Manuscript docx check (revision/Final_files)

- **Supplementary_Information_revision_clean.docx**: S1 caption uses post-bolus ketamine CADSS (**r = 0.04**, **p = 0.86**). Session-averaged post-infusion (`CADSS_total_post_avg`) is a different sensitivity (r ≈ −0.056, p ≈ 0.78).
- **Main_revision_clean.docx**: Steiger, calib-temp χ²/F, NPS ρ, demographics consistent with repo.
- **rebuttal_letter.docx**: CADSS descriptives, χ²/F, pain LRT consistent with repo.

The string `0.777` still appears in the supplementary docx inside **Supplementary Table S14** (S2(L) model **R²** column), not S1.

## Figure S1 variable

Supplementary Figure S1 uses **`CADSS_total_post_session2`** (ketamine post-bolus CADSS total). Session-averaged post-infusion (`CADSS_total_post_avg`) is not the submitted S1.

## Not fully re-run in CI

- Full `01_pain_ketamine_analysis_temp_covariate.r` (long runtime; committed tables are reference).
- `02_supplementary_revision.Rmd` render requires **pandoc ≥ 1.12.3**.
- Supplementary Figure S2: embedded asset only.

## Minor repo hygiene (optional)

- S1 primary outputs: `Figure_S_dose_equivalence_weight_cadss.{png,pdf}`; legacy copies `Figure_S_calibtemp_pain_slopes.{png,pdf}`.
- `text/npp_revision_2026/` may lag `revision/Final_files/`; sync before OSF upload.
