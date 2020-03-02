###########################################################################
# 
# FILE NAME: get_uci_datasets.R
#
# FILE DESCRIPTION: Download UCI datasets.
#                   Downloaded and computed data files will be stored in folder ./Data/.
#
###########################################################################

# Downloaded (.csv and .dat) and computed (.RData) data files will be stored here.
datafolder <- "./../Data/"


# Init --------------------------------------------------------------------

rm(list = setdiff(ls(), "datafolder"))


library(dplyr)
library(magrittr)
library(readr)


# Red wine ----------------------------------------------------------------

# Download.
download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv", 
              paste0(datafolder, "winequality-red.csv"), 
              method = "auto", 
              quiet = FALSE)

# Import.
red_wine <- read_delim(paste0(datafolder, "winequality-red.csv"), 
                        ";", 
                       escape_double = FALSE, 
                       trim_ws = TRUE) %>%
  mutate(quality = as.numeric(quality > 6))

# Create X.
X <- red_wine %>% 
  dplyr::select(-quality) %>% 
  model.matrix(~-1 + ., data = .) 

# Remove constant columns.
X <- X[, which(apply(X, 2, function(x) !all(x == x[1])))] 

# Create y.
y <- red_wine %>% dplyr::select(quality) %>% as.matrix()

# Save.
save(file = paste0(datafolder, "red_wine.RData"), X, y)

# Clean-up.
rm(red_wine, X, y)


# White wine --------------------------------------------------------------

# Download.
download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv", 
              paste0(datafolder, "winequality-white.csv"), 
              method = "auto", 
              quiet = FALSE)

# Import. 
white_wine <- read_delim(paste0(datafolder, "winequality-white.csv"), 
                         ";", 
                         escape_double = FALSE, 
                         trim_ws = TRUE) %>%
  mutate(quality = as.numeric(quality > 6))

# Create X.
X <- white_wine %>% 
  dplyr::select(-quality) %>% 
  model.matrix(~-1 + ., data = .) 

# Remove constant columns.
X <- X[, which(apply(X, 2, function(x) !all(x == x[1])))] 

# Create y.
y <- white_wine %>% dplyr::select(quality) %>% as.matrix()

# Save
save(file = paste0(datafolder, "white_wine.RData"), X, y)

# Clean-up.
rm(white_wine, X, y)


# Statlog: Australian credit approval -------------------------------------

# Download.
download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/statlog/australian/australian.dat", 
              paste0(datafolder, "australian.dat"), 
              method = "auto", 
              quiet = FALSE)

# Import.
australian <- read_table2(paste0(datafolder, "australian.dat"), 
                          col_names = FALSE, col_types = cols(X1 = col_factor(levels = c("0", "1")), 
                                                              X11 = col_factor(levels = c("1", "0")), 
                                                              X12 = col_factor(levels = c("1", "2", "3")), 
                                                              X4 = col_factor(levels = c("1", "2", "3")), 
                                                              X5 = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")), 
                                                              X6 = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9")), 
                                                              X8 = col_factor(levels = c("1", "0")), 
                                                              X9 = col_factor(levels = c("1", "0"))))

# Create X.
X <- australian %>% 
  dplyr::select(-X15) %>% 
  model.matrix(~-1 + ., data = .) 

# Remove constant columns.
X <- X[, which(apply(X, 2, function(x) !all(x == x[1])))] 

# Create y.
y <- australian %>% dplyr::select(X15) %>% as.matrix()

# Save.
save(file = paste0(datafolder, "australian.RData"), X, y)

# Clean-up.
rm(australian, X, y)


# Statlog: German credit data ---------------------------------------------

# Download.
download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/german.data-numeric", 
              paste0(datafolder, "german.dat"), 
              method = "auto", 
              quiet = FALSE)

# Import.
german <- read.table(paste0(datafolder, "german.dat"), 
                     quote = "\"", 
                     comment.char = "", 
                     stringsAsFactor = FALSE) 

# Create X.
X <- german %>% 
  dplyr::select(-V25) %>% 
  model.matrix(~-1 + ., data = .)

# Remove constant columns.
X <- X[, which(apply(X, 2, function(x) !all(x == x[1])))] 

# Create y.
y <- german %>% dplyr::select(V25) %>% as.matrix() - 1

# Save.
save(file = paste0(datafolder, "german.RData"), X, y)

# Clean-up.
rm(german, X, y)


# Abalone -----------------------------------------------------------------

# Download.
download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data", 
              paste0(datafolder, "abalone.dat"), 
              method = "auto", 
              quiet = FALSE)

# Import.
abalone <- read_csv(paste0(datafolder, "abalone.dat"), 
                    col_names = FALSE, col_types = cols(X1 = col_factor(levels = c("M", "F", "I")))) %>%
  mutate(X9 = as.numeric(X9 > 10))

# Create X.
X <- abalone %>% 
  dplyr::select(-X9) %>% 
  model.matrix(~-1 + ., data = .)

# Remove constant columns.
X <- X[, which(apply(X, 2, function(x) !all(x == x[1])))] 

# Create y.
y <- abalone %>% dplyr::select(X9) %>% as.matrix()

# Save.
save(file = paste0(datafolder, "abalone.RData"), X, y)

# Clean-up.
rm(abalone, X, y)


# Final clean-up ----------------------------------------------------------

rm(list = ls())