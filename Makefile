# Run from repository root: make -C submission/osf  OR  cd submission/osf && make

ROOT_DIR := $(CURDIR)
export ROOT_DIR

.PHONY: all main supplementary dose-equiv verify check sync-s3

all: main supplementary dose-equiv sync-s3

sync-s3:
	@cp -f output/figures/connectivity_all_networks.png output/revision/figures/Supplementary_Figure_S3_connectivity.png 2>/dev/null || true
	@cp -f output/figures/connectivity_all_networks.pdf output/revision/figures/Supplementary_Figure_S3_connectivity.pdf 2>/dev/null || true

verify:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/verify_manuscript_numbers.R

main:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/01_pain_ketamine_analysis_temp_covariate.r

supplementary:
	cd code && ROOT_DIR="$(ROOT_DIR)" Rscript -e 'knitr::purl("02_supplementary_revision.Rmd", output = "_supp_run.R", quiet = TRUE); source("_supp_run.R"); file.remove("_supp_run.R")'

dose-equiv:
	Rscript code/revision/plot_dose_equivalence.R

check: verify
	@test -f output/revision/figures/Figure_S_dose_equivalence_weight_cadss.png
	@test -f output/revision/tables/steiger_roi_session_corr_diff.csv
	@echo "Key revision outputs present."
