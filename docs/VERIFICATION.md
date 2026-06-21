# Manuscript number verification

Run from the repository root.

## Quick check

```bash
make verify
```

Runs `code/verify_manuscript_numbers.R` against committed outputs (~seconds). This is the recommended check after `make all`.

## Full cross-check

```bash
make verify-all
```

This command:

1. Rebuilds all analyses (`make all`)
2. Parses numeric values from the manuscript Word files in `text/manuscript/` (Python inventory only)
3. Compares supplementary table cells and in-text statistics to R outputs (`code/verify_all_manuscript.R`)

All statistical comparisons are done in R. Python is used only to read the Word documents.

Results are printed to the terminal. Ephemeral working files are written to `output/verification/` (not committed).

Table → CSV mapping: [`TABLE_MANIFEST.yaml`](TABLE_MANIFEST.yaml).

## Regenerate individual outputs

```bash
make dose-equiv    # Supplementary Figure S1
make supp-fig-s2   # Supplementary Figure S2
make sync-s5       # Supplementary Table S5 (from data reference)
```

## Scope

| Included | Not included |
|----------|----------------|
| Behavioral mixed models, ROI summaries, connectivity | SPM first-level fMRI preprocessing |
| Steiger tests, calibration-temperature models | CADSS infusion-timing figures |
| Supplementary Figures S1–S3 | Arrow-task supplementary figures |

Supplementary Figure S1 uses **`CADSS_total_post_session2`** (ketamine post-bolus CADSS total).
