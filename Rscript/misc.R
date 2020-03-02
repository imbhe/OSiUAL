#############################################################################
#
# FILE NAME: misc.R
# 
# FILE DESCRIPTION: Miscellaneous functions.
#
# FUNCTIONS:
#
#    - safe_cv_glmnet: Quiet and error safe wrapper around glmnet::cv.glmnet. 
#
#    - compute_cs: Computes calibration slope, using a quiet and error safe wrapper around the glm function.
#
# REQUIREMENTS: glmnet and purrr packages.
#
#############################################################################


# Fit penalised glm using glmnet::cv.glmnet.
# Input: See ?glmnet::cv.glmnet.
# Output: See ?glmnet::cv.glmnet. 
# Returns NULL if erorrs are encountered in glmnet::cv.glmnet.
safe_cv_glmnet <- function(...) { 
  purrr::safely(glmnet::cv.glmnet, otherwise = NULL, quiet = TRUE)(...)[["result"]]
}


# Compute calibration slope. 
# Input: binary response vector y and predicted probabilities pred.
# Output: calibration slope.
# Returns NA if errors are encountered in glm.
compute_cs <- function(y, pred) {
  purrr::safely(function(y, pred) {
    coef(glm(y ~ I(log(pred / (1  - pred))), family = "binomial"))[2]
  }, 
  otherwise = NA, 
  quiet = TRUE)(y, pred)$result
}