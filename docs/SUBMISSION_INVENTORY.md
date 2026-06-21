# NPP revision submission inventory (R2)

## Authoritative submission files (Box)

**Source of truth for what gets submitted:**

`/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/revision/revesion_2/Final_files`

| File | Role |
|------|------|
| `Main_R2.docx` | Main manuscript (clean) |
| `Main_R2_tc.docx` | Main manuscript (track-changes export) |
| `Supplementary_Information_revision_clean.docx` | Supplementary Information |
| `Rebuttal_R2.docx` | Response to reviewers (round 2) |

Sync into this repo for verification:

```bash
make sync-r2-final    # copies Final_files → text/npp_revision_2026_r2/
make patch-s4-docx    # if Table S4 needs refresh from R
make verify-all
```

## Git snapshot (verification copies)

[`text/npp_revision_2026_r2/`](../text/npp_revision_2026_r2/) — working copies used by `make extract-claims` and `make verify-all`. **Not authoritative** until synced from `Final_files/` above.

Legacy R1 snapshots remain in `text/npp_revision_2026/`.

## Main manuscript figures

| Asset | Committed output |
|-------|------------------|
| Figure 1–4 | `output/revision/figures/Main_Figure_{1-4}.pdf` |

## Supplementary figures

| Asset | Reproducible | Committed output |
|-------|--------------|------------------|
| S1 (weight × post-bolus CADSS) | Yes | `output/revision/figures/Figure_S_dose_equivalence_weight_cadss.*` |
| S2 (ROI × calib temp) | Yes | `make supp-fig-s2` | `output/revision/figures/Supplementary_Figure_S2_roi_calib_temperature.*` |
| S3 (connectivity) | Yes | `output/revision/figures/Supplementary_Figure_S3_connectivity.*` |

## Supplementary tables S1–S14

See [`TABLE_MANIFEST.yaml`](TABLE_MANIFEST.yaml) and [`REPRODUCTION_MANIFEST.md`](REPRODUCTION_MANIFEST.md).

Verification: `make verify-all` (full cell + in-text checks against R2 docx).

## Out of repo scope

- Exploratory `revision/tables/` bulk exports not cited in R2 Final_files
- Raw DICOM / fMRIPrep / first-level fMRI
