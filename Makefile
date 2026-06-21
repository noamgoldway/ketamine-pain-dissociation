# Run from repository root: cd ketamine-pain-dissociation && make

ROOT_DIR := $(CURDIR)
export ROOT_DIR

.PHONY: all main supplementary dose-equiv verify verify-all check sync-s3 sync-s5 supp-fig-s2 export-s1 export-s4

all: main supplementary dose-equiv sync-s3 sync-s5 supp-fig-s2 export-s4

export-s4:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/export_supp_table_s4.R
	ROOT_DIR="$(ROOT_DIR)" Rscript code/export_supp_table_s4_word_subset.R

export-s1:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/export_supp_table_s1.R

sync-s5:
	@test -f data/spm_ketamine_clusters_reference.csv || (echo "Missing data/spm_ketamine_clusters_reference.csv" && exit 1)
	cp -f data/spm_ketamine_clusters_reference.csv output/revision/tables/supp_table_s5_ketamine_clusters.csv

supp-fig-s2:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/revision/plot_roi_calib_temperature.R

verify-all: all
	python3 code/extract_manuscript_claims.py
	ROOT_DIR="$(ROOT_DIR)" Rscript code/export_supp_table_s1.R
	ROOT_DIR="$(ROOT_DIR)" Rscript code/verify_all_manuscript.R

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
	@test -f output/revision/figures/Supplementary_Figure_S2_roi_calib_temperature.png
	@test -f output/revision/tables/steiger_roi_session_corr_diff.csv
	@echo "Key revision outputs present."
