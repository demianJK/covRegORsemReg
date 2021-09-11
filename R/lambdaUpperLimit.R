#' determine maximal lambda for Bien & Tiebshirani (2011) covLASSO estimator or van Wieringen & Wessel (2016) covRidge estimator
#'
#' @param data raw data set
#' @param lambdaMax guess for maximal lambda
#' @param stop stop criterion (default: 80% of off-diagonal values, i.e., covariances, shrunken to zero)
#' @param type covLASSO or covRidge estimator
#'
#' @import covglasso
#' @import matrixcalc
#' @import rags2ridges
#'
#' @return congruence relative percentage of congruence between sample covariance matrix and target matrix
#' @return lambdaMax attained lambdaMax


compareL <- function(S_L, N, lambdaMax){
  covLASSO <- covglasso::covglasso(S=S_L, n=N, lambda=lambdaMax)$sigma # compute regularized matrix with lambdaMax
  compL <- covLASSO[lower.tri(covLASSO)] == 0 # compare all off-diagonal elements of lower triangle (symmetric matrix) with maximally regularized condition (i.e., zeros)
  # compute relative percentage of congruence:
  abs <- 0
  nelm <- length(compL) # number of unique off-diagonale elements in covmat
  for (i in 1:nelm){
    if (compL[i]){
      abs <- abs + 1
    }
  }
  rel <- round(abs/nelm, 1)
  # if (all(compL)){ # only TRUEs
  #   t <- unname(table(compL)[1])
  #   f <- 0
  # } else if (length(table(compL)) == 1){ # also the case in all(compL) but this was tested in condition before and there is no other way to check if all are FALSE
  #   t <- 0
  #   f <- unname(table(compL)[1])
  # } else { # TRUEs and FALSEs
  #   t <- unname(table(compL)[2])
  #   f <- unname(table(compL)[1])
  # }
  # rel <- round(t / (f + t), 1) # accuracy: first decimal place
  return(rel)
}

compareR <- function(S_R, T_RP, lambdaMax){
  covRidgeP <- rags2ridges::ridgeP(S=S_R, lambda=lambdaMax) # compute regularized matrix with lambdaMax
  compR <- round(covRidgeP, 1) == round(T_RP, 1) # compare regularized matrix with target (both in precision terms); liberal rounding
  # compute relative percentage of congruence:
  if (all(compR)){ # only TRUEs
    t <- unname(table(compR)[1])
    f <- 0
  } else if (length(table(compR)) == 1){ # only FALSEs (opposite of: at least one TRUE)
    t <- 0
    f <- unname(table(compR)[1])
  } else { # TRUEs and FALSEs
    t <- unname(table(compR)[2])
    f <- unname(table(compR)[1])
  }
  rel <- round(t / (f + t), 1) # accuracy: first decimal place
  return(rel)
}

lambdaUpperLimit <- function(data, lambdaMax, stop=0.8, type="covLASSO"){
  # sanity checks
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop("data must be a matrix or data frame")
  }
  else if (lambdaMax <= 0) {
    stop("lambdaMax should be positive")
  }
  else if (stop < 0 || stop > 1) {
    stop("stop criterion should be a decimal number between 0 and 1")
  }
  else if (!(type %in% c("covLASSO", "covRidge"))){
    stop("type should be one of {'covLASSO', 'covRidge'}")
  }

  N <- nrow(data)
  if (type == "covLASSO"){
    data <- scale(data, center = TRUE, scale = FALSE) # centering necessary
    S_L <- ((N-1)/N)*cov(data) # biased cov estim
    if (!matrixcalc::is.positive.definite(S_L)){
      S_L <- ridgeOption(S_L)
      warning("ridge option evoked to attain positive definite matrix for covLASSO estimation")
    }
    rel <- compareL(S_L, N, lambdaMax)

    if (rel < stop){ # lamda_max too small
      return(list(congruence=rel, lambdaMax=lambdaMax))
      stop("Please supply a bigger lambdaMax value.")
    } else if (rel == stop){ # lambdaMax approximately fine
      return(list(congruence=rel, lambdaMax=lambdaMax))
    } else { # lambdaMax too big
      while (rel > stop){ # "pre"-loop (to determine first lower bound of lambdaMax)
        lambdaMax_test <- lambdaMax/2    # lambdaMax_test = 5   lambdaMax = 10
        rel <- compareL(S_L, N, lambdaMax_test)
        lambdaMax_higher <- lambdaMax # buffer variable      # lambdaMax_higher = 10
        lambdaMax <- lambdaMax_test # for next iteration      # lambdaMax = 5     (am Ende: lambdaMax = 5   lambda_test=5   lambda_higher = 10)
      }
      # actual loop:
      while (rel != stop){
        lambdaMax_lower <- 0 # used as determination criterion
        if (rel < stop){ # after pre-loop first (rel$rel1 < stop) is always true
          lambdaMax_test <- (lambdaMax_higher + lambdaMax)/2
          rel <- compareL(S_L, N, lambdaMax_test)
          lambdaMax_lower <- lambdaMax # buffer variable
          lambdaMax <- lambdaMax_test # for next iteration
        } else { # rel > stop
          lambdaMax_test <- (lambdaMax_lower + lambdaMax)/2
          rel <- compareL(S_L, N, lambdaMax_test)
          lambdaMax_higher <- lambdaMax
          lambdaMax <- lambdaMax_test
        }
        if (round(lambdaMax_lower,1) == round(lambdaMax_higher,1)){ # determination criterion
          return(list(congruence=rel, lambdaMax=lambdaMax, determined="TRUE"))
        }
      }
      return(list(congruence=rel, lambdaMax=lambdaMax))
    }
  } else { # type == "covRidge"
    data <- scale(data, center = TRUE, scale = FALSE) # centering necessary
    S_R <- rags2ridges::covML(as.matrix(data))
    T_RP <- rags2ridges::default.target(S_R) # data will be centered and covmat computed via cross-prod approach
    rel <- compareR(S_R, T_RP, lambdaMax) # compute regularized matrix and relative percentage of zero values in off-diagonal triangle

    if (rel < stop){ # lamda_max too small
      return(list(congruence=rel, lambdaMax=lambdaMax))
    } else if (rel == stop){ # lambdaMax approximately fine
      return(list(congruence=rel, lambdaMax=lambdaMax))
    } else { # lambdaMax too big
      # "pre"-loop (to determine first lower bound of lambdaMax):
      while (rel > stop){
        lambdaMax_test <- lambdaMax/2
        rel <- compareR(S_R, T_RP, lambdaMax_test)
        lambdaMax_higher <- lambdaMax
        lambdaMax <- lambdaMax_test
      }
      # actual loop:
      while (rel != stop){
        lambdaMax_lower <- 0 # used as determination criterion
        if (rel < stop){ # after pre-loop first (rel$rel1 < stop) is always true
          lambdaMax_test <- (lambdaMax_higher + lambdaMax)/2
          rel <- compareR(S_R, T_RP, lambdaMax_test)
          lambdaMax_lower <- lambdaMax
          lambdaMax <- lambdaMax_test
        } else { # rel$rel1 > stop
          lambdaMax_test <- (lambdaMax_lower + lambdaMax)/2
          rel <- compareR(S_R, T_RP, lambdaMax_test)
          lambdaMax_higher <- lambdaMax
          lambdaMax <- lambdaMax_test
        }
        if (round(lambdaMax_lower,1) == round(lambdaMax_higher,1)){ # determination criterion
          return(list(congruence=rel, lambdaMax=lambdaMax, determined="TRUE"))
        }
      }
      return(list(congruence=rel, lambdaMax=lambdaMax))
    }
  }
}
