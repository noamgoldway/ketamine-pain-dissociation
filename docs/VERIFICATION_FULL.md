# Full verification report (R2)

Generated from R outputs in `/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/submission/osf`.

Pipeline: `make verify-all` runs **`make all`** (main + supplementary + figures + static sync) before numeric checks.

Source docx: `text/npp_revision_2026_r2/`

## Summary

- **Passed:** 26
- **Failed:** 0

### Supplementary tables + anchors

- OK: S1: 84/84 cells match connectivity_roi_list.csv
- OK: S2: 21/21 cells match supp_table_cadss_mixed_model.csv
- OK: S3: 15/15 cells match cadss_posthoc_session_contrasts.csv
- OK: S4: 60/60 cells match cadss_posthoc_timepoint_contrasts.csv
- OK: S4: row-level estimate/t/p/session checks OK
- OK: S5: 40/40 cells match supp_table_s5_ketamine_clusters.csv
- OK: S6: 21/21 cells match supp_table_roi_activation_model_anova.csv
- OK: S7: 35/35 cells match roi_diff_posthoc.csv
- OK: S8: 45/45 cells match roi_m2_nice_anova.csv
- OK: S9: 8/8 cells match posthoc_m2_calibtemp_slope_by_intensity_noadj.csv
- OK: S10: 8/8 cells match posthoc_m2_calibtemp_slope_by_session_noadj.csv
- OK: S11: 16/16 cells match posthoc_m2_calibtemp_slope_by_intensity_within_session_noadj.csv
- OK: S12: 28/28 cells match pain_roi_corr_results.csv
- OK: S13: 49/49 cells match pain_roi_dissoc_model_comparison.csv
- OK: S14: 4/4 cells match nps_dissoc_model_comparison.csv
- OK: supp_table_tempadj_pain_model: calib_temp p = 0.070 (got 0.07)
- OK: Main abstract pain interaction F(1,54)=11.22: got 11.22
- OK: Main abstract CADSS session F(1,32)=57.44: got 57.44
- OK: ROI rho ant_insula_L placebo ≈ 0.528
- OK: ROI rho ant_insula_L ketamine ≈ 0.656
- OK: ROI rho dlPFC_R ketamine ≈ 0.57
- OK: pain_calibration_stats.csv present for calibration temperature descriptives
- OK: Supp Fig S2 regenerated: placebo slope 0.331 > ketamine 0.131
- OK: Supp Fig S2 slopes consistent with submitted figure (0.331/0.131 vs 0.341/0.131)
- OK: S5 cluster table synced from data/spm_ketamine_clusters_reference.csv
- OK: verify_manuscript_numbers.R spot-checks passed

### R spot-check suite (`verify_manuscript_numbers.R`)

OK

```
OK   S1_n
OK   S1_pearson_r
OK   S1_pearson_p
OK   demog_weight_mean
OK   demog_weight_sd
OK   steiger_min_p_above_0.11
OK   steiger_dlPFC_p
OK   pain_calib_temp_p
OK   cadss_post_bolus_mean
OK   cadss_post_bolus_sd
OK   cadss_post_bolus_min
OK   cadss_post_bolus_max
OK   nps_high_pain_placebo_rho
OK   nps_high_pain_ketamine_rho
OK   roi_calib_temp_chisq
OK   roi_calib_temp_chisq_p
OK   roi_intensity_session_calib_F
OK   roi_intensity_session_calib_p
OK   pain_calib_temp_chisq
OK   pain_calib_temp_chisq_p

All manuscript spot-checks passed.
```

**Overall:** PASS
