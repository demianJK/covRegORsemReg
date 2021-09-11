#' simulate data for a MIMIC model
#'
#' @param p_pred number of predictors
#' @param p_pred_zero amount of regression weights being zero (remaining regression weights are 0.2, 0.5 and 0.8; equal parts)
#' @param p_pred_r predictor correlation (equals predictor covariance)
#' @param N sample size
#' @param seed seed for random number generation
#'
#' @import lavaan
#'
#' @return simulated data set

simulateMIMIC <- function(p_pred, p_pred_zero, p_pred_r, N, seed)
  {
  # sanity checks:
  # if (p_pred <= 0) {
  #   stop("number of predictors (p_pred) should be positive")
  # }  else if (N <= 0) {
  #   stop("sample size (N) should be positive")
  # }  else if (seed <= 0) {
  #   stop("(seed) should be positive")
  # }
  if (p_pred_zero == "40%"){
    np <- c(0.4*p_pred, 0.6*p_pred, 0.8*p_pred)
  } else if (p_pred_zero == "70%"){
    np <- c(0.7*p_pred, 0.8*p_pred, 0.9*p_pred)
  }

  # create MIMIC population model:
  factor_loadings <- c("f1 =~ 1*y1 + 0.8*y2 + 0.8*y3 + 0.8*y4 + 0.5*y5 + 0.5*y6") # number of indicators is fixed
  factor_variance <- c("f1 ~~ 1*f1") # fix variance to 1 (for model identification), actually already done in default simulateData(..., std.lv=T)
  # possible regression weights are fixed: 0, 0.2, 0.5, 0.8
  # relative number of predictors with certain regression weights (np)
  predictors <- paste("f1 ~ ",
                      paste0(paste0("0*", "x", 1:np[1], collapse="+"), # predictors with 0 weights
                            "+",
                            paste0("0.2*", "x", (np[1]+1):np[2], collapse="+"), # predictors with 0.2 weights
                            "+",
                            paste0("0.5*", "x", (np[2]+1):np[3], collapse="+"), # predictors with 0.5 weights
                            "+",
                            paste0("0.8*", "x", (np[3]+1):p_pred, collapse="+") # predictors with 0.8 weights
                            ))
  variances <- paste("x", 1:p_pred, "~~", "1*", "x", 1:p_pred, sep="", collapse=";") # variances are set to unity by default in simulateData()
  # since predictor variances are fixed to 1, predictor correlation = covariance
  # when no fixed values are specified --> set to zero (e.g. means!)

  covariances = list()
  count=0
  for(i in 1:p_pred){   ## Kovarianzen in der Form: x1 ~~ 0*x2
    for(j in 1:p_pred){
      if(i != j & j > i){
        count = count+1
        covariances[count] = paste("x", i, "~~", p_pred_r, "*","x", j, sep="")
      }
    }
  }
  covariances = paste(covariances, collapse=";")
  pop.model <- paste(factor_loadings, factor_variance, predictors, variances, covariances, sep="\n")

  # simulate
  dat <- lavaan::simulateData(pop.model, sample.nobs=N, seed=seed, model.type="lavaan")

  return(dat)
}
