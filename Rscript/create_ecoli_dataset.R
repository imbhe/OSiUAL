###########################################################################
# 
# FILE NAME: create_ecoli_dataset.R
#
# FILE DESCRIPTION: Create E. coli data set.
#                   Assumes that FASTA and GenBank files ec.chr01.fa and ec.chr01.gb are available 
#                   and located in the folder ./Data/.
#                   Features (transition frequencies) and labels (coding vs non-coding sequence) 
#                   are extracted and stored in ./Data/ecoli.RData. 
#                   Necessary FASTA and GenBank files may be downloaded from 
#                   https://www.ncbi.nlm.nih.gov/nuccore/15829254.
#
###########################################################################

# Data located and stored here.
datafolder <- "./../Data/"


# Init --------------------------------------------------------------------

rm(list = setdiff(ls(), "datafolder"))

library(dplyr)
library(magrittr)
library(purrr)
library(readr)
library(stringr)
library(stringi)
library(tidyr)


# Functions ---------------------------------------------------------------

# Compute number of fist order transitions of in genomic sequence seq.
compute_transitions <- function(seq) {
  acgt <- c("A", "C", "G", "T")
  trans <-  do.call(paste0, expand.grid(acgt, acgt)) %>% str_sort()
  names(trans) <- trans
  count <- map_dfc(trans, function(x) str_count(seq, x))
  names(count) <- paste(names(trans), "count", sep = "_")
  return(count)
}


# Create E. coli dataset -----------------------------------------------------------------

# Read FASTA file with genomic sequence. Should be located in folder specified by datafolder.
fa <- read_file(paste0(datafolder, "ec.chr01.fa")) %>%
  str_remove(str_sub(., 1, str_locate(.,"\n")[1])) %>%
  str_remove_all("\n")

# Read GenBank file with annotations, extract start and stop of coding sequences (CDS), i.e. genes.
cds <- read_delim(paste0(datafolder, "ec.chr01.gb"), 
                         "§", 
                         escape_double = FALSE, 
                         col_names = FALSE, 
                         trim_ws = TRUE,
                  col_types = cols(X1 = col_character())) %>%
  rename(str = X1) %>%
  filter(str_detect(str, "CDS ") & !str_detect(str, "join")) %>% # remove spliced exons.
  mutate(str = str_remove(str, "CDS ")) %>%
  mutate(strand = ifelse(str_detect(str, "complement"), "reverse", "forward")) %>%
  mutate(str = str_remove(str, "complement"),
         str = str_remove(str, "\\("),
         str = str_remove(str, "\\)")) %>%
  separate(str, into = c("start", "stop"), sep = "\\.\\.") %>%
  mutate(start = as.numeric(start), 
         stop = as.numeric(stop)) %>%
  mutate(type = "exon")

# Find non-coding sequences (intergenes).
non_cds <- tibble(start = c(1, cds$stop + 1), stop = c(cds$start - 1,  str_length(fa)), type = "intergene", strand = "forward") %>% 
  filter(stop >= start)

# Create data in "GFF" format.
gff <- bind_rows(cds, non_cds) %>%
  arrange(start, stop) %>%
  mutate(seq = str_sub(fa, start, stop),
         seq = ifelse(type == "exon" & strand == "reverse",
                      stri_reverse(str_to_upper(str_replace_all(seq, c("T" = "a", "C" = "g", "G" = "c", "A" = "t")))), 
                      seq),
         length = str_length(seq)) %>%
  dplyr::select(type, strand, start, stop, length, seq) %>%
  filter(length > 6)

# Create data with relative frequency of transitions in coding and non-coding sequences.
transition_data <- gff %>%
  dplyr::select(seq) %>%
  map_dfr(compute_transitions) %>%
  bind_cols(gff) %>%
  dplyr::select(type, strand, start, stop, length, seq, everything()) %>%
  gather(key = "transition", value = "count", contains("count")) %>%
  arrange(start, transition) %>%
  group_by(start) %>%
  mutate(transition = str_remove(transition, "_count"),
         relfreq = count / sum(count)) %>%
  dplyr::select(-seq, -count) %>%
  spread(key = "transition", value = "relfreq") %>%
  ungroup() 

# Create X.
X <- transition_data %>%
  dplyr::select(-one_of("type", "strand", "start", "stop", "length")) %>%
  as.matrix()

# Remove redundant column.
X <- X[,-1]

# Create y.
y <- as.numeric(transition_data$type == "exon")

# Save.
save(file = paste0(datafolder, "ecoli.RData"), X, y)

# Clean-up.
rm(list = ls())