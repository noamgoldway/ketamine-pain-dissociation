# Data and output inventory

Analysis inputs in `data/` and committed outputs in `output/`.

## Analysis inputs (`data/`)

| File | Description |
|------|-------------|
| `pain_ratings.csv` | Trial-level pain ratings |
| `CADSS.csv` | CADSS scores by session and time point |
| `NPS.csv` | Neurophysiological pain signature scores |
| `roi_beta_values_by_condition.csv` | ROI activation betas |
| `within_connectivity.xlsx` | Functional connectivity summaries |
| `pain_calibration.xlsx` | Calibration temperatures and pain slopes |
| `demog.csv` | Demographics |
| `participants.csv` | Participant IDs and metadata |
| `prior_ketamine_use.csv` | Prior ketamine exposure |
| `CADSS_Weight_DoseEquivalence_Data.csv` | Supplementary Figure S1 |
| `connectivity_roi_list.csv` | ROI list for connectivity analyses |
| `spm_ketamine_clusters_reference.csv` | Reference clusters for Supplementary Table S5 |

## Manuscript files (`text/manuscript/`)

Main text and supplementary information Word documents, used as reference for `make verify-all`.

## Key outputs

| Output | Script |
|--------|--------|
| `output/tables/` | `make main` |
| `output/figures/` | `make main` |
| `output/revision/tables/` | `make supplementary`, `make export-s4` |
| `output/revision/figures/` | `make dose-equiv`, `make supp-fig-s2`, `make sync-s3` |

Verification: `make verify` (spot-check) or `make verify-all` (full cross-check). See [`VERIFICATION.md`](VERIFICATION.md).
