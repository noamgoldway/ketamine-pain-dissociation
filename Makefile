# Run from repository root: make -C submission/osf  OR  cd submission/osf && make

ROOT_DIR := $(CURDIR)
export ROOT_DIR

.PHONY: all main supplementary dose-equiv check

all: main supplementary dose-equiv

main:
	ROOT_DIR="$(ROOT_DIR)" Rscript code/01_pain_ketamine_analysis_temp_covariate.r

supplementary:
	ROOT_DIR="$(ROOT_DIR)" Rscript -e 'rmarkdown::render("code/02_supplementary_revision.Rmd", quiet = TRUE)'

dose-equiv:
	Rscript code/revision/plot_dose_equivalence.R

check:
	@test -f output/revision/figures/Figure_S_calibtemp_pain_slopes.png
	@test -f output/revision/tables/steiger_roi_session_corr_diff.csv
	@echo "Key revision outputs present."
