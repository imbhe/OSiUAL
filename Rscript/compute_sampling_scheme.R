#############################################################################
#
# FILE NAME: compute_sampling_scheme.R
# 
# FILE DESCRIPTION: Function that computes selection probabilities for various sampling schemes.
#
# INPUT:
#   - X: Data matrix, intercept excluded.
#   - fit: L2 penalised (ridge) logistic regresion model fitted using cv.glmnet.
#   - sampling_scheme: Sampling scheme. Valid options are 
#       * "prop1" (optimised to minimise the anticipated variance of the estimated loss).
#       * "cor1a" (optimised to minimise the anticipated generalisation error in terms of the total loss of the active learning algorithm).
#       * "cor1b" (optimised to minimise the anticipated mean squared error of the predictions).
#       * "prob_un" (probabilistic uncertainty sampling).
#       * "det_un" (deterministic unceratinty sampling).
#       * "uniform" (uniform or simple random samplingpassive learning, i.e. passive learning).
#
# OUTPUT: Sampling probability for each observation (row) in X.
#
# REQUIREMENTS: glmnet package.
#
#############################################################################


compute_sampling_scheme <- function(X, fit, sampling_scheme) {
  
  sampling_scheme <- tolower(sampling_scheme) # Make case insensitive.
  lambda <- fit$lambda.min # Regularisation parameter.
  N <- nrow(X)
  
  # Predictions.
  if (!is.matrix(X)) { X <- as.matrix(X) }
  pred <- as.numeric(predict(fit, X, s = fit$lambda.min, type = "response"))
  
  # Add column of ones to X (corresponding to intercept).
  X <- cbind(1, X) 
  
  # Compute selection probabilities.
  if (sampling_scheme == "prop1") {
    
    pi <- sqrt(pred * log(pred)^2 + (1 - pred) * log(1 - pred)^2)
    
  } else if (sampling_scheme == "cor1a") {
    
    var_y <- pred * (1 - pred)
    V <- diag(var_y)
    XVX <- t(X) %*% V %*% X
    
    # Inverse of Hessian of total loss.
    # No penatly on intercept parameter. 
    H_inv <- solve(XVX / N + diag(c(0, rep(lambda, ncol(X) - 1)))) 
    
    pi <- vapply(1:N, function(ix) sqrt(var_y[ix] * t(X[ix, ]) %*% H_inv %*% X[ix, ]), numeric(1))
    
  } else if (sampling_scheme == "cor1b") {
    
    var_y <- pred * (1 - pred)
    V <- diag(var_y)
    XVX <- t(X) %*% V %*% X
    
    # Inverse of Hessian of total loss.
    # No penatly on intercept parameter. 
    H_inv <- solve(XVX / N + diag(c(0, rep(lambda, ncol(X) - 1)))) 
    
    HXV2XH <- H_inv %*% t(X) %*% V^2 %*% X %*% H_inv # pxp matrix.
    pi <- vapply(1:N, function(ix) sqrt(var_y[ix] * t(X[ix, ]) %*% HXV2XH %*% X[ix, ]), numeric(1))
    
  } else if (sampling_scheme == "prob_un") {
    
    pi <- -(pred * log(pred) + (1 - pred) * log(1 - pred))
    
  } 
  
  # Ad hoc adjustment if computations produce NA or infinite values.
  pi[is.na(pi)] <- min(pi[is.finite(pi)], na.rm = TRUE)
  pi[is.infinite(pi)] <- min(pi[is.finite(pi)], na.rm = TRUE)
  
  # Normalise sampling probabilities to sum to 1.
  pi <- pi / sum(pi)
  
  return(pi)
}