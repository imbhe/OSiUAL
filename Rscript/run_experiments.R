#############################################################################
#
# FILE NAME: run_experiments.R
#
# FILE DESCRIPTION: Main program. Runs experiments in 
#                   Optimal sampling in unbiased active learning
#                   presented at AISTATS 2020
#                   by Henrik Imberg, Johan Jonasson and Marina Axelson-Fisk.
#
#############################################################################

datafolder <- "./../Data/" # Input datasets stored here.
datasets <- c("abalone", "australian", "ecoli", "german", "red_wine", "white_wine") # Dataset names (.RData files, without file extension).
respath <- "./../Results/" # Output datasets with active learning performance will be stored here.


# Init --------------------------------------------------------------------

rm(list = setdiff(ls(), c("datafolder", "datasets", "respath")))

library(magrittr)
library(purrr)
library(tibble)
library(tidyr)
library(dplyr)

# Assumes that R scripts al.R, compute_sampling_scheme.R and misc.R are located in the current folder.
source("al.R")
source("compute_sampling_scheme.R")
source("misc.R")


# Prepare -----------------------------------------------------------------

# I/O.
io <- tibble(datafolder = datafolder, respath = respath, dsname = datasets) %>%
  mutate(infile = paste0(datafolder, dsname, ".RData"),
         outfile = paste0(respath, dsname, ".RData")) 

# Create dataset with execution information.
# See al.R for additional information about these parameters.
params <- crossing(n_init = 25, # Size of initial sample.
                   n_per_step = 25, # Number of instances to query in each iteration of the active learning algorithm.
                   n_final = 500, # Final sample size.
                   n_reps = 10000, # Number of experiments.
                   sampling_scheme = c("prop1", "cor1a", "cor1b", "prob_un", "det_un", "uniform"), # Sample seletion procedures. 
                   n_cores = 50, # Number of cores.
                   seed = 123456) %>% # Set seed, for reproducibility.
  crossing(io, .) %>% 
  dplyr::select(-datafolder, -respath)


# Run ---------------------------------------------------------------------

print("--- run_experiments.R executing. ---")
exec <- pmap_dbl(params, al) # exec has return value 1 upon successful completion.



# Compute and save summary statistics -------------------------------------

print("Saving summary statistics...")

# Get filenames.
files <- paste0(respath, list.files(respath, pattern="*.RData"))
names(files) <- as.character(seq_along(files))

# Load files.
all_results <- map_dfr(files, function(file) {
  load(file)
  df <- bind_cols(file = rep(file, nrow(res)), res)
  return(df)
})

# Compute average of performance metrics.
summary_results <- all_results %>%
  mutate(dsname = str_replace_all(dsname, "_", " "), 
         dsname = str_to_title(dsname),
         dsname = str_replace_all(dsname, "Ecoli", "E. coli"),
         dsname = ifelse(dsname == "German", paste(dsname, "Credit Data"), dsname),
         dsname = ifelse(dsname == "Australian", paste(dsname, "Credit Approval"), dsname)) %>%
  group_by(dsname, sampling_scheme, n) %>%
  summarise(auc = mean(auc, na.rm = TRUE),
            cs = median(cs, na.rm = TRUE),
            cil = mean(cil, na.rm = TRUE),
            misclass = 1 - mean(accuracy, na.rm = TRUE),
            nll = median(nll, na.rm = TRUE),
            rmse_pred = mean(rmse_pred, na.rm = TRUE),
            tpr = mean(tpr, na.rm = TRUE)) %>% 
  ungroup() 


# Save results for optimal model, i.e. using entire dataset for training.
optmod_results <- summary_results %>% 
  filter(is.na(sampling_scheme)) %>% 
  dplyr::select(-sampling_scheme, -n, -rmse_pred) # Drop columns.

save(optmod_results, file = paste0(respath, "optmod_results.RData"))

# Save active learning results.
summary_results <- summary_results %>% 
  filter(!is.na(sampling_scheme)) 

save(summary_results, file = paste0(respath, "summary_results.RData"))


# Clean-up ----------------------------------------------------------------

print(ifelse(all(exec == 1), "Execution completed without errors.", "Execution completed with errors."))

rm(list = ls())