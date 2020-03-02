###########################################################################
#
# FILE NAME: print_session_info.R
# 
# FILE DESCRIPTION: Write session info to file in current folder.
#
###########################################################################

# Load all packages used in this project.
library("bigstatsr")
library("cowplot")
library("doParallel")
library("dplyr")
library("ggplot2")
library("glmnet")
library("magrittr")
library("plyr")
library("pROC")
library("purrr")
library("readr")
library("stringr")
library("stringi")
library("tibble")
library("tictoc")
library("tidyr")
library("xtable")

# Write sessionInfo to file.
writeLines(capture.output(sessionInfo()), "./../sessionInfo.txt")