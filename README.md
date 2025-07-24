# ketamine-pain-dissociation
Analysis code and input data for the paper The Analgesic and Dissociative Properties of Ketamine are Separate and Correspond to Distinct Neural Mechanisms (Goldway et al.). 

📂 Contents
pain_ketamine_analysis.R: Main analysis script that:

Computes activation differences across high/low pain conditions under ketamine.

Tests whether dissociation (CADSS subscales) explains additional variance in brain activation.

Assesses the relationship between pain ratings and neural signals (e.g., NPS, DMN).

input_data/: Includes anonymized subject-level input data:

fMRI ROI activation values

Subjective pain ratings (placebo & ketamine)

CADSS dissociation questionnaire scores

🔍 Key Finding Summary
Brain activation in pain-sensitive regions is robustly predicted by subjective pain but not dissociation.

DMN disintegration correlates with dissociation, not analgesia.

📜 Requirements
R (≥ 4.0.0)

Packages: tidyverse, broom, purrr, dplyr, ggplot2, readr, tidyr

📄 License
MIT License (or your preferred license)

