#' ridge option: add small positive constant to diagonal of covariance matrix until it is positive definite
#'
#' @param S covariance matrix estimate
#' @param const small positive number added to diagonal elements of S in every iteration until matrix is pd
#'

ridgeOption <- function(S, const=1e-5, type="ex"){
  # 1e-5 ist scheinbar Wert aus lavaan Ridge option
  # https://github.com/yrosseel/lavaan/blob/master/R/lav_options.R
  # if (matrixcalc::is.positive.definite(S)){
  #   stop("matrix is positive-definite; ridge option does not need to be applied")
  # }

  if (!(type %in% c("ex", "all"))){
    stop("type should be one of {'ex, 'all'}")
  }

  if (type == "ex"){
    diag(S)[6:ncol(S)] <- diag(S)[6:ncol(S)] + const # only ridging exogenuous variables once
    # see:
    # https://github.com/yrosseel/lavaan/blob/master/R/lav_samplestats_icov.R
    return(S)
  } else {
    diag(S)<- diag(S)+ const
    return(S)
  }
  return(S)
}


# diag_const <- diag(const, ncol = ncol(S), nrow = nrow(S))
# while(any(eigen(S)$values < 0)){
#   S <- S + diag_const
# }

# alternative loop:
# while(!is.positive.definite(covMat_ridgeReg)){
#   covMat_ridgeReg <- covMat_ridgeReg + diag_const
# }
## beide Alternativen scheinen äquivalent (weil gleiche cond numb cond_ridgeReg)
