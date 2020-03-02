# README

This repository contains code used in the paper
Optimal sampling in unbiased active learning
presented at AISTATS 2020
by Henrik Imberg, Johan Jonasson and Marina Axelson-Fisk.

---

## Preparations

1. Run R script install_packages.R to install the R packages used in this project.
2. Download the FASTA and GenBank files for the E. coli genome from https://www.ncbi.nlm.nih.gov/nuccore/15829254. Store the files as ec.chr01.fa and ec.chr01.gb in the Data folder.
3. Run the script create_ecoli_dataset.R to extract features (relative transition frequencies) and labels (coding vs non-coding sequence) and prepare the E. coli dataset for analysis. The resulting R data file (ecoli.RData) is stored in the Data folder.
4. Run get_uci_datasets.R to download datasets from the UCI machine learning repository and prepare for analysis. The resulting R data files are stored in the Data folder.

## Experiments

5. The R script run_experiments.R runs the experiments. Assumes that analysis-ready R data files are available and located in the Data folder. Complete results (named dataset_algorithm.RData) and summary statistics (summary_results.RData and optmod_results.RData) are stored in the Results folder. Change the number of cores (n_cores) as appropriate. To reduce computation time, you may wish to reduce the maximal training size (n_final, currently set to 500) and/or the number of experiments (n_reps, currently set to 10 000).

## Results

6. Run R script tables_and_figures.RData to generate tables and figures. Assumes that analysis-ready R data files are available and located in the Data folder (needed for Table S1), and that R data files with summary results (summary_results.RData and optmod_results.RData) are available and located in the Results folder. Output is saved as .pdf (figures) or .tex (tables) files in the Output folder.

## Miscellanea

- al.R implements various active learning algorithms. Used by run_experiments.R.
- compute_sampling_scheme.R implements and computes selection probabilities for various sampling schemes. Used by al.R.
- misc.R implements quiet and error-safe wrappers around the glm and glmnet::cv.glmnet functions, using the purrr:safely function. Used by al.R.
- print_session_info.R prints session info (OS, R and package info) to the file sessionInfo.txt.

## Additional information.

Additional details are found in the file headers.
