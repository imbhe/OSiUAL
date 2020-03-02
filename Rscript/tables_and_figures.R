#############################################################################
#
# FILE NAME: tables_and_figures.R
#
# FILE DESCRIPTION: Creates tables and figures in 
#                   Optimal sampling in unbiased active learning
#                   presented at AISTATS 2020
#                   by Henrik Imberg, Johan Jonasson and Marina Axelson-Fisk.
#
#############################################################################

# Input datasets stored here.
datafolder <- "./../Data/" 

# Results of active learning algorithm stored here.
respath <- "./../Results/" 

# Output (tables and figures) will be stored here.
outpath <- "./../Output/" 


  
# Init --------------------------------------------------------------------

rm(list = setdiff(ls(), c("outpath", "respath", "datafolder", "summary_results", "optmod_results")))

library(cowplot)
library(dplyr)
library(ggplot2)
library(glmnet)
library(MASS)
library(purrr)
library(RColorBrewer)
library(stringr)
library(tibble)
library(tidyr)
library(xtable)


# My ggplot theme.
theme_set(theme_bw()) 
theme_update(axis.line = element_line(colour = "black", size = 0.3), 
             axis.text = element_text(family = "serif", size = 8, colour = "black"),
             axis.ticks = element_line(colour = "black", size = 0.3),
             axis.title = element_text(size = 10),
             legend.key.height = unit(0.4, "cm"),
             legend.key.width = unit(1.2, "cm"),
             legend.margin = margin(t = -0.25, unit = 'cm'),
             legend.spacing =  unit(0, "cm"),
             legend.position = "bottom",
             legend.text = element_text(family = "serif", size = 8),
             legend.text.align = 0,
             legend.title = element_text(size = 10),
             panel.border = element_blank(),
             panel.grid = element_line(size = 0.1),
             panel.grid.major = element_line(size = 0.1),
             plot.title = element_text(hjust = 0.5, family = "serif", size = 10, colour = "black"),
             strip.background.x = element_blank(),
             strip.text = element_text(hjust = 0.5, family = "serif", size = 10, colour = "black"),
             text = element_text(family = "serif", size = 10, colour = "black"))


# Load results.
load(paste0(respath, "summary_results.RData"))
load(paste0(respath, "optmod_results.RData"))


# Figure 1 ----------------------------------------------------------------

# For reproducibility.
set.seed(1)

# Some parameters.
N <- 10^3 # Size of dataset.
beta <- c(0, rep(1, 5)) # Regression coefficients.


# Simulate data.
X <- cbind(1, mvrnorm(n = N, mu = rep(0, 5), Sigma = diag(1, 5))) # Simulate X.
y <- rbinom(n = N, size = 1, p = as.numeric(1 / (1 + exp(-X %*% beta))) ) # Simulate Y.

# Fit model and compute various quantities of interest.
mod <- glm(y~ -1 + X, family = "binomial") # Fit model.
pred <- fitted.values(mod) # Predictions.
var_y <- pred * (1 - pred) # Label uncertainty (variance of Y_i)
V <- diag(var_y) # V-matrix, variance of Y_i on diagonal. 
H <- t(X) %*% V %*% X # Hessian of total loss. Alternative: solve(vcov(mod))
Hinv <- solve(H) # Inverse Hessian. Alternative: vcov(mod)

# Probabilistic uncertainty sampling (Chu et al., 2011, Ganti & Gray, 2012)
pi1 <- pred * log(pred) + (1 - pred) * log(1 - pred)
pi1 = pi1 / sum(pi1) # Normalise: sum up to 1.

# Proposition 1: minimise anticipated variance of estimated loss.
pi2 <- pred * log(pred)^2 + (1 - pred) * log(1 - pred)^2
pi2 = pi2 / sum(pi2) # Normalise: sum up to 1.

# Corollary 1a: minimise expectation of total loss.
pi3 <- vapply(1:N, function(ix) var_y[ix] * t(X[ix,]) %*% Hinv %*% X[ix,], numeric(1))
# Alternative: pi3 <- hatvalues(mod) 
pi3 <- pi3 / sum(pi3) # Normalise: sum up to 1.

# Corollary 1b: minimise anticipated mean squared error of predictions.
L <- V %*% X %*% Hinv # Compute once.
pi4 <- vapply(1:N, function(ix) {v <- sqrt(var_y[ix]) * L %*% X[ix,]; sqrt(t(v) %*% v)}, numeric(1))
pi4 <- pi4 / sum(pi4) # Normalise: sum up to 1.

df <- tibble(pred = pred, pi1 = pi1, pi2 = pi2, pi3 = pi3, pi4 = pi4) 

# For Corollary 1: plot subset to avoid overplotting.
dfs <- df %>%
  filter(row_number() %in% sample(1:N, 100))  

labels <- c("Chu et al. (2011), Ganti & Gray (2012)", "Proposition 1", "Corollary 1a", "Corollary 1b")
colours <- brewer.pal(n = 9, "Set1")

# Probabilistic uncertainty sampling (Chu et al. ,2011; Ganti & Gray, 2012) and sampling according to Proposition 1 and Corollary 1.
fig1 <- ggplot(df, aes(x = pred)) +
  geom_line(aes(y = pi1, colour = "Chu et al. (2011), Ganti & Gray (2012)", linetype = "Chu et al. (2011), Ganti & Gray (2012)"), lwd = 0.75) + 
  geom_line(aes(y = pi2, colour = "Proposition 1", linetype = "Proposition 1"), lwd = 0.75) +
  geom_point(data = dfs, aes(y = pi3, colour = "Corollary 1a", shape = "Corollary 1a"), size = 1) +
  geom_point(data = dfs, aes(y = pi4, colour = "Corollary 1b", shape = "Corollary 1b"), size = 1) +
  expand_limits(x = c(0, 1)) + 
  scale_y_continuous(breaks = seq(0, max(c(pi1, pi2, dfs$pi3, dfs$pi4)), length.out = 3), labels = c("0", rep("", 2))) + 
  scale_colour_manual(breaks = letters[1:4], values = colours, limits = labels) + 
  labs(x = "P(Y = 1)",
       y = "Sampling probability", 
       colour = NULL,
       shape = NULL,
       linetype = NULL) +   
  guides(colour = guide_legend(ncol = 1),
         linetype = guide_legend(ncol = 1, override.aes = list(colour = colours[1:2])),
         shape = guide_legend(ncol = 1, override.aes = list(size = 1.5, colour = colours[3:4]))) +
  theme(axis.ticks.y = element_blank(),
        legend.box = "vertical",
        legend.box.just = "left",
        legend.key.width = unit(0.75, "cm"),
        legend.margin = margin(t = 0, b = 0, unit = 'cm'))



# Clean-up
rm(df, dfs, H, Hinv, L, mod, V, X, beta, N, pi1, pi2, pi3, pi4, pred, var_y, y)



# Figure 2 and S.1 - S7 -----------------------------------------------------

# Set colours, linetypes and labels.
breaks <- c("prop1", "cor1a", "cor1b", "prob_un", "det_un", "uniform")
labels <- c("Prop. 1", "Cor. 1a", "Cor. 1b", "Prob. un.", "Det. un.", "Uniform")
colours <- brewer.pal(n = 8, name = "Dark2")[2:7]
linetypes <- c("73", "42", "13", "22", "1343", "solid")
  
plotdata <- summary_results %>%
  filter(n <= 250) # Plot results for sample sizes up to n = 250.


# Figure 2 --------

# Figure 2a. Misclassification rate on Abalone, E. coli and Red Wine datasets. 
fig2a <- plotdata %>%
  filter(dsname %in% c("Abalone", "E. coli", "Red Wine")) %>% 
  ggplot(aes(x = n, y = misclass, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_hline(data = optmod_results %>% filter(dsname %in% c("Abalone", "E. coli", "Red Wine")), 
             aes(yintercept = misclass), lwd = 0.25, col = "grey60", lty = 1) + # Result using entire pool for training.
  geom_line(lwd = 0.75) + 
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  scale_y_continuous(breaks = c(seq(0, 0.15, 0.01), seq(0.15, 1, 0.02))) + 
  labs(x = "Size of labelled sample",
       y = "Misclassification rate",
       colour = NULL,
       linetype = NULL) 


# Figure 2b. negative log-likelihood of predictions on Abalone, E. coli and Red Wine datasets. 
fig2b <- plotdata %>%
  filter(dsname %in% c("Abalone", "E. coli", "Red Wine")) %>% 
  ggplot(aes(x = n, y = nll, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_line(lwd = 0.75) + 
  geom_hline(data = optmod_results %>% filter(dsname %in% c("Abalone", "E. coli", "Red Wine")), 
             aes(yintercept = nll), lwd = 0.25, col = "grey60", lty = 1) + # Result using entire pool for training.
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Size of labelled sample",
       y = "(Scaled) Negative log-likelihood",
       colour = NULL,
       linetype = NULL) 


# Combine.
fig2 <- plot_grid(fig2a + theme(axis.title.x = element_blank(), 
                                legend.position = "none"), 
                  fig2b + theme(strip.text = element_blank(), 
                                legend.position = "none"), 
                  get_legend(fig2a),
                  nrow = 3, rel_heights = c(1, 1, 0.1))

rm(fig2a, fig2b)


# Figure S1. Misclassification rate, all datasets. --------

figS1 <- ggplot(plotdata, aes(x = n, y = misclass, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_hline(data = optmod_results, aes(yintercept = misclass), lwd = 0.25, col = "grey60", lty = 1) + # Result using entire pool for training.
  geom_line(lwd = 0.75) + 
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  labs(x = "Size of labelled sample",
       y = "Misclassification rate",
       colour = NULL,
       linetype = NULL) 


# Figure S2. Proportion correctly classified minority examples, all datasets. --------

figS2 <- ggplot(plotdata, aes(x = n, y = tpr, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_hline(data = optmod_results, aes(yintercept = tpr), lwd = 0.25, col = "grey60", lty = 1) + # Result using entire pool for training.
  geom_line(lwd = 0.75) + 
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Size of labelled sample",
       y = "Proportion correctly classified minority examples",
       colour = NULL,
       linetype = NULL) 


# Figure S3. AUC, all datasets. --------

figS3 <- ggplot(plotdata, aes(x = n, y = auc, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_hline(data = optmod_results, aes(yintercept = auc), lwd = 0.25, col = "grey60", lty = 1) + # Result using entire pool for training.
  geom_line(lwd = 0.75) + 
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Size of labelled sample",
       y = "AUC",
       colour = NULL,
       linetype = NULL) 



# Figure S4. Negative log-likelihood, all datasets. --------

figS4 <- ggplot(plotdata, aes(x = n, y = nll, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_hline(data = optmod_results, aes(yintercept = nll), lwd = 0.25, col = "grey60", lty = 1) + # Result using entire pool for training.
  geom_line(lwd = 0.75) + 
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Size of labelled sample",
       y = "(Scaled) Negative log-likelihood",
       colour = NULL,
       linetype = NULL) 



# Figure S5. RMSE, all datasets. --------

figS5 <- ggplot(plotdata, aes(x = n, y = rmse_pred, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_line(lwd = 0.75) + 
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Size of labelled sample",
       y = "RMSE of predictions",
       colour = NULL,
       linetype = NULL) 



# Figure S6. Calibration-in-the-large: ratio of observed vs predicted (expected) number of minority examples, all datasets. --------

figS6 <- ggplot(plotdata, aes(x = n, y = cil, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_hline(aes(yintercept = 1), lwd = 0.25, col = "grey60", lty = 1) + # Reference value for well-calibrated model.
  geom_line(lwd = 0.75) + 
  coord_cartesian(ylim = c(0.8, 1.2)) +
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels,values = linetypes) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Size of labelled sample",
       y = "Observed / predicted number of minority examples",
       colour = NULL,
       linetype = NULL) 


# Figure S7. Calibration slope, all datasets. --------

figS7 <- ggplot(plotdata, aes(x = n, y = cs, colour = sampling_scheme, lty = sampling_scheme)) + 
  geom_hline(aes(yintercept = 1), lwd = 0.25, col = "grey60", lty = 1) + # Reference value for well-calibrated model.
  geom_line(lwd = 0.75) + 
  coord_cartesian(ylim = c(0.5, 1.5)) +
  facet_wrap(~dsname, scales = "free_y", dir = "h", ncol = 3) + 
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_colour_manual(breaks = breaks, labels = labels, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Size of labelled sample",
       y = "Calibration slope",
       colour = NULL,
       linetype = NULL) 



# Table S1 -----------------------------------------------------------------

describe_dataset <- function(filepath, dsname){
  load(file = filepath)
  
  # Recode outcome so that 1 is minority class and 0 is majority class. 
  if ( mean(y) > 0.5 ) { y <- 1 - y } 

  # Number of observations.
  N <- sprintf("%s", nrow(X)) 

  # Proportion of observations in minority class.
  n_min <- sprintf("%d", sum(y == 1))
  pct_min <- sprintf("(%.1f\\%%)", 100 * mean(y == 1))
  
  # Number of features.
  n_features = ncol(X)
  
  # Get results when using entire dataset for training.
  res <- optmod_results[which(str_detect(optmod_results$dsname, dsname)), ]

  # Combine in dataset.
  summary <- tibble(vals = c(N, n_min, pct_min, n_features, "", 
                             sprintf("%.2f", 1 - res$misclass), 
                             sprintf("%.2f", res$auc)))
  names(summary) <- dsname
 
  return(summary)
  
}

# Data info.
data_info <- tibble(datafolder = datafolder, 
                    filename = c("abalone.RData","australian.RData", "ecoli.RData", "german.RData", "red_wine.RData", "white_wine.RData"),
                    dsname = c("Abalone", "Australian", "E. coli", "German", "Red Wine", "White Wine")) %>% 
  mutate(filepath = paste0(datafolder, filename))

# Row labels of Table 1.
row_labels <- tibble(Feature = c("Number of records",
                                 "n (\\%) in minority class",
                                 "",
                                 "Number of predictors$^a$",                           
                                 "Performance under optimal model$^b$",
                                 "\\hspace{0.3cm} Accuracy$^c$",
                                 "\\hspace{0.3cm} AUC")) 

# Compute dataset summary statistics.
tabS1 <- data_info %>% 
  dplyr::select(dsname, filepath) %>% 
  pmap_dfc(describe_dataset) %>% 
  bind_cols(row_labels, .)


# Clean-up
rm(data_info, describe_dataset, row_labels)



# Table S2 and Figure S8 - S11 --------------------------------------------


# Set colours, linetypes and labels.
breaks <- c("prop1", "cor1a", "cor1b", "prob_un", "det_un", "uniform")
labels <- c("Prop. 1", "Cor. 1a", "Cor. 1b", "Prob. un.", "Det. un.", "Uniform")
colours <- brewer.pal(n = 8, name = "Dark2")[2:7]
linetypes <- c("73", "42", "13", "22", "1343", "solid")


compute_label_complexity <- function(dataset, method, metric, dir, n_ref) {
  
  # Get active and passive learning results.
  al <- summary_results %>% 
    filter(dsname == dataset & sampling_scheme == method)
  
  pl <- summary_results %>%
    filter(dsname == dataset & sampling_scheme == "uniform")
  
  # Sample sizes.
  nseq_al <- al$n
  nseq_pl <- pl$n
  
  # Get performance.
  al_perf <- as.numeric(unlist(al[, metric]))
  pl_perf <- as.numeric(unlist(pl[, metric]))
  
  # Interpolate results for sample sizes between the ones actually used.
  ipol_al <- approx(x = nseq_al, y = al_perf, n = diff(range(nseq_al)) + 1)
  ipol_pl <- approx(x = nseq_pl, y = pl_perf, n = diff(range(nseq_pl)) + 1)
  
  n_al <- vapply(n_ref, function(n_pl) {
    
    # Performance using passive learning of sample size n_pl.
    y <- ipol_pl$y[which(ipol_pl$x == n_pl)] 
    
    # Return smallest sample size that produces better performance with active learning.
    return(ipol_al$x[which(get(dir)(ipol_al$y, y))[1]])
    
  }, numeric(1))
  
  
  res <- tibble(dsname = dataset, 
                sampling_scheme = method,
                metric = metric,
                n_pl = n_ref, 
                n_al = n_al,
                ss_ratio = n_al / n_pl)
  
  return(res)
  
}

params <- tibble(metric = c("auc", "misclass", "rmse_pred", "nll"), 
                 dir = c(">=", "<=", "<=", "<=")) %>% 
  crossing(., dataset = unique(summary_results$dsname), 
           method = setdiff(unique(summary_results$sampling_scheme), "uniform"),
           n_ref = seq(25, 500, 25))


label_complexity <- params %>% 
  pmap_dfr(compute_label_complexity) %>% 
  mutate(n_al = ifelse(is.na(n_al), max(n_pl + 1), n_al)) %>% # %>% # If label complexity for active learning missing (out of range): set to maximal passive learning sample size + 1)
  mutate(ss_ratio = ifelse(is.na(ss_ratio), sprintf("$>$%.2f", max(summary_results$n) / n_pl), sprintf("%.2f", ss_ratio))) # If label complexity for active learning missing (out of range): set ratio to > maximal sample size / passive learning sample size.
  

plotdata <- label_complexity %>%
  mutate(dsname = ifelse(dsname == "German Credit Data", "German", dsname),
         dsname = ifelse(dsname == "Australian Credit Approval", "Australian", dsname)) %>% 
  filter(n_pl <= 250) # Plot results for passive learning up to size n = 250. 


# Figure S8. Label complexity of active vs. passive learning, misclassification rate. ----
figS8 <- ggplot(plotdata %>% filter(metric == "misclass"), aes(x = n_pl, y = n_al, colour = sampling_scheme, lty = sampling_scheme)) +
  geom_abline(aes(intercept = 0, slope = 1), lwd = 0.25, col = "grey60") +
  geom_line(lwd = 0.75) + 
  coord_equal(xlim = c(25, 250), ylim = c(0, 250)) +
  facet_wrap(~dsname, dir = "h", ncol = 3) + 
  scale_colour_manual(labels = labels, breaks = breaks, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  labs(x = "Size of labelled sample, passive learning",
       y = "Size of labelled sample, active learning",
       colour = NULL,
       linetype = NULL) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) 


# Figure S9. Label complexity of active vs. passive learning, AUC. ----
figS9 <- ggplot(plotdata %>% filter(metric == "auc"), aes(x = n_pl, y = n_al, colour = sampling_scheme, lty = sampling_scheme)) +
  geom_abline(aes(intercept = 0, slope = 1), lwd = 0.25, col = "grey60") +
  geom_line(lwd = 0.75) + 
  coord_equal(xlim = c(25, 250), ylim = c(0, 250)) +
  facet_wrap(~dsname, dir = "h", ncol = 3) + 
  scale_colour_manual(labels = labels, breaks = breaks, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  labs(x = "Size of labelled sample, passive learning",
       y = "Size of labelled sample, active learning",
       colour = NULL,
       linetype = NULL) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) 



# Figure S10. Label complexity of active vs. passive learning, negative log-likelihood of predictions. ----
figS10 <- ggplot(plotdata %>% filter(metric == "nll"), aes(x = n_pl, y = n_al, colour = sampling_scheme, lty = sampling_scheme)) +
  geom_abline(aes(intercept = 0, slope = 1), lwd = 0.25, col = "grey60") +
  geom_line(lwd = 0.75) + 
  coord_equal(xlim = c(25, 250), ylim = c(0, 250)) +
  facet_wrap(~dsname, dir = "h", ncol = 3) + 
  scale_colour_manual(labels = labels, breaks = breaks, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  labs(x = "Size of labelled sample, passive learning",
       y = "Size of labelled sample, active learning",
       colour = NULL,
       linetype = NULL) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) 



# Figure S11. Label complexity of active vs. passive learning, RMSE of predictions. ----
figS11 <- ggplot(plotdata %>% filter(metric == "rmse_pred"), aes(x = n_pl, y = n_al, colour = sampling_scheme, lty = sampling_scheme)) +
  geom_abline(aes(intercept = 0, slope = 1), lwd = 0.25, col = "grey60") +
  geom_line(lwd = 0.75) + 
  coord_equal(xlim = c(25, 250), ylim = c(0, 250)) +
  facet_wrap(~dsname, dir = "h", ncol = 3) + 
  scale_colour_manual(labels = labels, breaks = breaks, values = colours) +
  scale_linetype_manual(breaks = breaks, labels = labels, values = linetypes) +
  labs(x = "Size of labelled sample, passive learning",
       y = "Size of labelled sample, active learning",
       colour = NULL,
       linetype = NULL) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) 



# Table S2. Label complexity of active vs. passive learning with n = 250 training examples.
tabS2 <- label_complexity %>% 
  filter(n_pl == 250) %>%
  dplyr::select(-n_al) %>%
  group_by(dsname, metric) %>%
  spread(key = "sampling_scheme", value = "ss_ratio") %>% 
  ungroup() %>%
  mutate(metric = str_replace_all(metric, "auc", "AUC"),
         metric = str_replace_all(metric, "nll", "Negative log-likelihood"),
         metric = str_replace_all(metric, "misclass", "Misclassification rate"),
         metric = str_replace_all(metric, "rmse_pred", "RMSE of predictions")) %>% 
  arrange(dsname, metric, n_pl) %>% 
  group_by(dsname) %>% 
  mutate(Dataset = ifelse(row_number() == 1, dsname, "")) %>% 
  group_by(dsname, metric) %>% 
  mutate(`Performance metric` = ifelse(row_number() == 1, metric, "")) %>% 
  ungroup() %>% 
  dplyr::select(Dataset, `Performance metric`, prop1, cor1a, cor1b, prob_un, det_un) %>% 
  dplyr::rename("Prop. 1" = prop1,
                "Cor. 1a" = cor1a,
                "Cor. 1b" = cor1b,
                "Prob. un." = prob_un,
                "Det. un." = det_un) 


rm(plotdata, label_complexity, params, compute_label_complexity)


# Save --------------------------------------------------------------------

ggsave(plot = fig1, paste0(outpath, "Figure1.pdf"), width = 8, height = 6.5, units = "cm")

ggsave(plot = fig2, paste0(outpath, "Figure2.pdf"), width = 16, height = 11, units = "cm")

ggsave(plot = figS1, paste0(outpath, "FigureS1.pdf"), width = 16, height = 10, units = "cm")

ggsave(plot = figS2, paste0(outpath, "FigureS2.pdf"), width = 16, height = 10, units = "cm")

ggsave(plot = figS3, paste0(outpath, "FigureS3.pdf"), width = 16, height = 10, units = "cm")

ggsave(plot = figS4, paste0(outpath, "FigureS4.pdf"), width = 16, height = 9, units = "cm")

ggsave(plot = figS5, paste0(outpath, "FigureS5.pdf"), width = 16, height = 9, units = "cm")

ggsave(plot = figS6, paste0(outpath, "FigureS6.pdf"), width = 16, height = 9, units = "cm")

ggsave(plot = figS7, paste0(outpath, "FigureS7.pdf"), width = 16, height = 9, units = "cm")

ggsave(plot = figS8, paste0(outpath, "FigureS8.pdf"), width = 16, height = 10, units = "cm")

ggsave(plot = figS9, paste0(outpath, "FigureS9.pdf"), width = 16, height = 10, units = "cm")

ggsave(plot = figS10, paste0(outpath, "FigureS10.pdf"), width = 16, height = 10, units = "cm")

ggsave(plot = figS11, paste0(outpath, "FigureS11.pdf"), width = 16, height = 10, units = "cm")


print(xtable(tabS1, 
             caption = "Descriptive statistics of benchmark datasets.",
             label = "tab:datasets"), 
      add.to.row = list(pos = list(7), command = "\\hline \\multicolumn{7}{l}{\\small $^{a}$ After re-coding of categorical predictors and removal of redundant variables.} \\\\ \\multicolumn{7}{l}{\\small $^{b}$ Using $L_2$-penalised logistic regression on the entire data set.} \\\\ \\multicolumn{7}{l}{\\small $^{c}$ Using 50\\% probability cut-off.} \\\\ \\multicolumn{7}{l}{\\small AUC, area under the receiver operating characteristic curve.} \\\\"),
      caption.placement = "top",
      include.rownames = FALSE, 
      sanitize.text.function = identity,
      file = paste0(outpath, "TableS1.tex"))

print(xtable(tabS2, 
             align = "lllccccc",
             caption = "Label complexity of active vs. passive learning, presented as the relative increase (ratio $>1$) or decrease (ratio $<1$) in the sample size needed for active learning to achieve equal performance as passive learning with $n = 250$ training examples.",
             label = "tab:label_complexity"),
      add.to.row = list(pos = list(-1, 24), command = c("\\multicolumn{2}{c}{ } & \\multicolumn{5}{c}{Ratio of sample sizes using active vs. passive learning} \\\\ \\cmidrule(lr){3-7}", "\\hline \\multicolumn{7}{l}{\\small AUC, area under the receiver operating characteristic curve; RMSE, root mean squared error.} \\\\")),
      caption.placement = "top",
      floating.environment = "table*",
      hline.after = seq(0, 24, 4),
      include.rownames = FALSE, 
      sanitize.text.function = identity,
      file = paste0(outpath, "TableS2.tex"))

# Clean-up
rm(list = ls())
