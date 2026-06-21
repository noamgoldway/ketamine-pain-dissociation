# Run from repository root: cd ketamine-pain-dissociation && make

ROOT_DIR := $(CURDIR)
export ROOT_DIR

.PHONY: all main supplementary verify check \
	supp-tables supp-figures dose-equiv supp-fig-s2 sync-s3 sync-s5

all: main supplementary

main:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/01_pain_ketamine_analysis_temp_covariate.r

supplementary: supp-tables supp-figures

supp-tables:
	cd code && ROOT_DIR="$(ROOT_DIR)" Rscript -e 'knitr::purl("02_supplementary.Rmd", output = "_supp_run.R", quiet = TRUE); source("_supp_run.R"); file.remove("_supp_run.R")'
	@$(MAKE) sync-s5

supp-figures: dose-equiv supp-fig-s2 sync-s3

dose-equiv:
	Rscript code/supplementary/plot_dose_equivalence.R

supp-fig-s2:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/supplementary/plot_roi_calib_temperature.R

sync-s3:
	@test -f output/figures/connectivity_all_networks.png || (echo "Run make main first (connectivity figure missing)." && exit 1)
	@cp -f output/figures/connectivity_all_networks.png output/supplementary/figures/Supplementary_Figure_S3_connectivity.png

sync-s5:
	@test -f data/spm_ketamine_clusters_reference.csv || (echo "Missing data/spm_ketamine_clusters_reference.csv" && exit 1)
	cp -f data/spm_ketamine_clusters_reference.csv output/supplementary/tables/supp_table_s5_ketamine_clusters.csv

verify:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/verify_manuscript_numbers.R

check: verify
	@test -f output/supplementary/figures/Supplementary_Figure_S1_weight_cadss_submitted.png
	@test -f output/supplementary/figures/Supplementary_Figure_S2_roi_calib_temperature.png
	@test -f output/supplementary/tables/steiger_roi_session_corr_diff.csv
	@echo "Key outputs present."
