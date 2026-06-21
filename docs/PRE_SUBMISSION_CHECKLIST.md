# Pre-submission checklist (R2)

Run from `submission/osf/`:

```bash
make sync-r2-final   # copy authoritative docx from Box Final_files
make verify-all      # must exit 0 before upload
make check           # key figures/tables on disk
```

**Authoritative Word files** live on Box at:

`revision/revesion_2/Final_files/`

(full path: `/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/revision/revesion_2/Final_files`)

Edit there first; then `make sync-r2-final` before verification.

## Automated verification

| Check | Status |
|-------|--------|
| `make all` (main + supplementary + figures) | ~14 s |
| `make verify-all` (strict S1–S14 + anchors) | **PASS** |
| `verify_manuscript_numbers.R` (19 spot-checks) | **PASS** |
| R2 docx in `text/npp_revision_2026_r2/` | Synced from `revision/revesion_2/Final_files/` |

## CADSS session coding (fixed 2026-05-19)

`code/01_pain_ketamine_analysis_temp_covariate.r` now maps `_1_` → Placebo, `_2_` →
Ketamine (consistent with pain ratings and raw `CADSS.csv`). Supplementary Table S4 in
the R2 supplementary docx was updated via `make patch-s4-docx`.

## External / not in verify-all

- SPM first-level fMRI (Table S5): `data/spm_ketamine_clusters_reference.csv`
- Main Figures 1–4 PDFs (committed exports)
- CADSS timing / arrow-task figures (out of scope)
