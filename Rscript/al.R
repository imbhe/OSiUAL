#############################################################################
#
# FILE NAME: al.R
#
# FILE DESCRIPTION: Function used to run active learning algorithm.
#
# INPUT:
#  - infile: Path to .RData file data containing data matrix X and binary response vector y taking values 0 and 1. 
#  - outfile: Filepath (without file extension) where results will be stored. Results will be stored as .RData file.
#  - dsname: Name of dataset. 
#  - n_init: Size of initial sample, selected using simple random sampling.
#  - n_per_step: Number of instances to query per iteration.
#  - n_final: Total number of instances to query, including initial sample.
#  - n_reps: Number of repetitions, i.e. the algorithm runs n_reps times.
#  - sampling_scheme: Sampling scheme. Valid options are 
#      * "prop1" (optimised to minimise the anticipated variance of the estimated loss).
#      * "cor1a" (optimised to minimise the anticipated generalisation error in terms of the total loss of the active learning algorithm).
#      * "cor1b" (optimised to minimise the anticipated mean squared error of the predictions).
#      * "prob_un" (probabilistic uncertainty sampling).
#      * "det_un" (deterministic unceratinty sampling).
#      * "uniform" (uniform or simple random samplingpassive learning, i.e. passive learning).
#  - n_cores: Number of cores. Repeated instances of the algorithm are run in parallell using forach loop provided by doParallel package.
#  - seed: Seed for random number generator. For reproducibility.
#
# OUTPUT: Dataset storing the result of the active learning algorithnm.
#
# REQUIREMENTS: bigstatsr, doParallel, glmnet, pROC, stringr and tictoc packages.
#
#############################################################################


al <- function(infile, 
               outfile, 
               dsname,
               n_init, 
               n_per_step, 
               n_final, 
               n_reps, 
               sampling_scheme,
               n_cores, 
               seed = 123456) {
  
  require(bigstatsr)
  require(doParallel)
  require(glmnet)
  require(pROC)
  require(stringr)
  require(tictoc)
  
  # Print execution info.
  print(sprintf("Running sampling scheme %s on %s", sampling_scheme, infile))
  
  # Start timer.
  tic()
  
  
  # Init ------------------------------------------------------------------
  
  # For reproducibility.
  set.seed(seed) 
  
  # Load and prepare data.
  load(file = infile)
  X <- scale(X)  # Scale features.
  if ( mean(y) > 0.5 ) { y <- 1 - y } # Recode outcome so that 1 is minority class and 0 is majority class. 
  y <- factor(y, levels = c(0, 1), labels = c("a", "b")) # Recode as factor. Necessary for glmnet.

  # For parallel computing.
  doParallel::registerDoParallel(cores = n_cores)
  
  # Make case insensitive.
  sampling_scheme <- tolower(sampling_scheme)
  
  # Some constants.
  N <- nrow(X) # Number of observations.
  n_iter <- (n_final - n_init) / n_per_step + 1 # Number of iterations.
  n_per_step <- c(n_init, rep(n_per_step, n_iter - 1)) # Number of instances to query in each iteration.
  nseq <- cumsum(n_per_step) # Cumulative number of instances queried.
  
  
  # Allocate matrices to store performance metrics.
  accuracy <- bigstatsr::FBM(nrow = n_iter, ncol = n_reps, init = NA) # Accuracy.
  auc <- bigstatsr::FBM(nrow = n_iter, ncol = n_reps, init = NA) # AUC.
  cil <- bigstatsr::FBM(nrow = n_iter, ncol = n_reps, init = NA) # Calibration-in-the-large, computed as the ratio of observed vs predicted number of instances in minority class.
  cs <- bigstatsr::FBM(nrow = n_iter, ncol = n_reps, init = NA) # Calibration slope. 
  nll <- bigstatsr::FBM(nrow = n_iter, ncol = n_reps, init = NA) # Negative log-likelihood.
  rmse_pred <- bigstatsr::FBM(nrow = n_iter, ncol = n_reps, init = NA) # Root mean squared error of predictions.
  tpr <- bigstatsr::FBM(nrow = n_iter, ncol = n_reps, init = NA) # True positive rate, i.e. proportion correctly classified minority examples. (Minority class coded as positive).
  
  
  # Performance using L2 penalised logistic regression on entire dataset. 
  # Penalty parameter chosen using 10-fold cross validation.
  fit0 <- glmnet::cv.glmnet(x = X, y = y, family = "binomial", alpha = 0, grouped = FALSE,
                            standardize = FALSE, intercept = TRUE) # Optimal fit.
  pred0 <-  as.numeric(predict(fit0, X, type = "response", s =  fit0$lambda.min)) # Predictions.
  accuracy0 <- mean((y == levels(y)[2]) == (pred0 > 0.5)) # Accuracy.
  auc0 <- as.numeric(pROC::auc(pROC::roc(y, pred0, levels = c("a", "b"), direction = "<"))) # AUC.
  cil0 <- sum(y == "b") / sum(pred0) # Calibration-in-the-large, computed as the ratio of observed vs predicted number of instances in minority class.
  cs0 <- compute_cs(y, pred0) # Calibration slope. 
  nll0 <- - mean((y == "b") * log(pred0) + (y == "a") * log(1 - pred0)) # Negative log-likelihood (scaled by a factor 1/N).
  tpr0 <- sum(pred0[y == "b"] > 0.5) / sum(y == "b") # True positive rate, i.e. proportion correctly classified minority examples. (Minority class coded as positive).
  
  
  # Active learning algorithm -------------------------------------------------------------
  
  # Parallellised active learning algorithm.
  # Results written to bigstatsr::FBM matrix initialised above.
  # Foreach loop gives return value 1 upon successful completion.
  tmp <- foreach(irep = 1:n_reps) %dopar% { 
    # tmp <- for(irep in 1:n_reps) { # For debugging.
    
    # Initialise.
    w <- rep(0, N) # Set sampling weights to zero.
    labelled <- integer(0) # Set collection of labelled instances to empty set.
    fit <- NULL # Model not yet fitted, set to NULL.
    
    # Iterate.
    for (iter in 1:n_iter) {
      
      nt <- n_per_step[iter] # Number to query in current iteration.
      bw <- nt / cumsum(n_per_step)[iter] # Batch-weight, proportional to batch size.
      
      if (is.null(fit)) {
        # If fit is NULL, i.e. model has not yet been fitted.
        
        # Use simple random sampling among unlabelled instances.
        new_instances <- sample(setdiff(1:N, labelled), nt)
        
        labelled <- c(labelled, new_instances) # Update pool of labelled and unlabelled instances.
        w <- rep(0, N)
        w[labelled] <- N / nseq[iter] # Update weights.
        
      } else { 
        # If fit is not NULL, i.e. model has sucessfully been fitted in earlier iteration.
        
        if (sampling_scheme %in% c("prop1", "cor1a", "cor1b", "prob_un")) {
          # If random sampling but not passive learning.
          
          pi <- compute_sampling_scheme(X = X, fit = fit, sampling_scheme = sampling_scheme)         
          
          # Sample nt instances from the pool of unlabelled instances.
          pi[labelled] <- 0 
          pi <- pi / sum(pi)
          new_wt <- rmultinom(n = 1, size = nt, prob = pi) / (nt * pi)
          
          # Re-query already labelled instances with probability 1. Comes with no additional cost since labels are already known. Needed to retain unbiasedness.
          new_wt[labelled] <- 1
          
          # Update weights.
          w <- w  + bw * (new_wt - w)
          
          # Update pool of labelled and unlabelled instances.
          labelled <- which(w > 0)
          
        } else if (sampling_scheme == "det_un") {
          # Deterministic uncertainty sampling among unlabelled instances.
          
          entropy <- - (pred * log(pred) + (1 - pred) * log(1 - pred)) # Compute entropy of label distribution.
          entropy[labelled] <- 0 # Set to zero for labelled instances to avoid re-queries.
          new_instances <- sort.list(entropy, decreasing = TRUE)[1:nt] # Sort according to label uncertainty, select most uncertain instances.
          
          labelled <- c(labelled, new_instances) # Update pool of labelled and unlabelled instances.
          w <- rep(0, N)
          w[labelled] <- 1 # Update weights.
          
        } else if (sampling_scheme == "uniform") {
          # Passive learning, i.e. simple random sampling among unlabelled instances.
          
          new_instances <- sample(setdiff(1:N, labelled), nt)
          
          labelled <- c(labelled, new_instances) # Update pool of labelled and unlabelled instances.
          w <- rep(0, N)
          w[labelled] <- 1 # Update weights.
          
        }
        
      } # End if fit is not NULL (line 97).
      
      # Retrieve features, labels and weights of labelled instances.
      wt <- w[labelled]
      Xt <- X[labelled, ]
      yt <- y[labelled]
      
      
      # Update model using L2 penalised logistic regression. Penalty parameter chosen using LOOCV or 10-fold cross validation.
      # Returns NULL if model cannot be fitted.
      new_fit <- safe_cv_glmnet(x = Xt, y = yt, weights = wt,
                                family = "binomial", alpha = 0, 
                                standardize = FALSE, intercept = TRUE, grouped = FALSE,
                                nfolds = ifelse(length(labelled) <= 50, length(labelled), 10))
      
      # If new_fit is not NULL, i.e model fitting succeeded without errors.
      if(!is.null(new_fit)) { 
        fit <- new_fit # Update model.
      } 
      
      # If fit is not NULL, i.e. model has been fitted in current or earlier iterations.
      if (!is.null(fit)) {
        
        # Compute predictions.
        pred <- as.numeric(predict(fit, X, type = "response", s = fit$lambda.min))
        
        # Evaluate performance. Results in current learning process stored in FBM matrices column-wise.
        accuracy[iter, irep] <- mean((y == levels(y)[2]) == (pred > 0.5)) # Accuracy
        auc[iter, irep] <- as.numeric(pROC::auc(pROC::roc(y, pred, levels = c("a", "b"), direction = "<"))) # AUC.
        cs[iter, irep] <- compute_cs(y, pred) # Calibration slope.
        cil[iter, irep] <- sum(y == "b") / sum(pred) # Calibration-in-the-large, computed as the ratio of observed vs predicted number of instances in minority class.
        nll[iter, irep] <- - mean((y == "b") * log(pred) + (y == "a") * log(1 - pred)) # Negative log-likelihood (scaled by a factor 1/N).
        rmse_pred[iter, irep] <- sqrt(mean((pred - pred0)^2)) # Root mean squared error of predictions.
        tpr[iter, irep] <- sum(pred[y == "b"] > 0.5) / sum(y == "b") # True positive rate, i.e. proportion correctly classified minority examples. (Minority class coded as positive).
        
      } 
      
    } # End current iteration.
    
    return(1) # Return value, foreach.
    
  } # End foreach.
  
  doParallel::stopImplicitCluster() # Stop cluster.
  
  
  # Save results ------------------------------------------------------------
  
  # Results of active learning algorithm.
  res <- tibble(dsname = dsname,
                N = N, 
                n_final = n_final,
                n_init = n_init,
                n_per_step = n_per_step[2],
                sampling_scheme = sampling_scheme,
                repetition = rep(1:n_reps, each = n_iter),
                n = rep(nseq, n_reps), 
                accuracy = as.numeric(accuracy[]), # Convert from FBM to vector, column-wise.
                auc = as.numeric(auc[]),
                cil = as.numeric(cil[]),
                cs = as.numeric(cs[]),
                nll = as.numeric(nll[]),
                rmse_pred = as.numeric(rmse_pred[]),
                tpr = as.numeric(tpr[])) 
  
  save(res, file = stringr::str_replace(outfile, ".RData", paste0("_", sampling_scheme, ".RData")))
  
  
  # Results of optimal model, i.e. using the entire dataset for training.
  res <- tibble(dsname = dsname, N = N, accuracy = accuracy0, auc = auc0, cs = cs0, cil = cil0, nll = nll0, tpr = tpr0)
  
  save(res, file = stringr::str_replace(outfile, ".RData", paste0("_optmod", ".RData")))
  
  
  # Stop timer.
  t <- toc(quiet = TRUE)
  print(sprintf("%.2f min elapsed.", as.numeric(t$toc - t$tic) / 60))
  
  # al function completed execution.
  return(1)
  
}
