LL <- function(S, P){ # aus rags2ridges:::.LL
  ##############################################################################
  # - Function that computes the value of the (negative) log-likelihood
  # - S > sample covariance matrix of validation set
  # - P > regularized covariance matrix of complement of validation set (i.e. training set)
  #
  #    --> fixed simulated data set, fixed data splitting, fixed lambda
  ##############################################################################
  LL <- -log(det(P)) + sum(S*(P^(-1))) #.trace(S %*% P)
  return(LL)
}


cvLL <- function(lambda, data, folds){ # aus rags2ridges:::.kcvl
  ##############################################################################
  # - Function that calculates a cross-validated negative log-likelihood score
  #   for single penalty value
  # - lambda > value penalty parameter
  # - Y      > (raw) Data matrix, variables in columns
  # - folds  > cross-validation sample splits
  #
  #    --> fixed data set, variable data splittings, fixed lambda 
  ##############################################################################
  
  dataC  <- scale(data, scale = FALSE) # Y centered (respective mean subtracted from X); rags2ridges::covML()
  
  cv_LL <- 0
  for (f in 1:length(folds)){
    ## P
    # dataset splitting Ai_complement (TRAIN)
    dataC_P <- dataC[-folds[[f]], , drop=FALSE] # drop=FALSE damit auch single row Output Dataframe bleibt
    n_P <- nrow(dataC_P)
    S_P <- ((n_P-1)/n_P)*cov(dataC_P) # biased estimator weil /n und nicht /(n-1)
    if (matrixcalc::is.positive.definite(S_P) == TRUE){
      P <- covglasso(S=S_P, n=n_P, lambda=lambda,
                     # Defaults: crit="bic", penalize.diag = FALSE, start=diag(cov(dat))
      )$sigma
    } else {
      S_ridgeReg <- S_P
      while(any(eigen(S_ridgeReg)$values<0)){
        S_ridgeReg <- S_ridgeReg + diag(1e-5, ncol = ncol(S_ridgeReg), nrow = nrow(S_ridgeReg))
      }
      P <- covglasso(S=S_ridgeReg, n=n_P, lambda=lambda)$sigma
    }
    
    ## S
    # dataset splitting Ai (TEST)
    Ai <- dataC[folds[[f]], , drop=FALSE]
    n <- nrow(Ai)
    S_Ai <- ((n-1)/n)*cov(Ai) # biased estimator
    
    cv_LL <- cv_LL + LL(P, S_Ai / length(folds[[f]])) # / length(folds[[f]]) um für n_k zu kontrollieren???
  } 
  
  return(cv_LL / length(folds))
}


lamda_Optim_cL <- function(data, # data set
                           lambdas, # vector with lambda values to be compared in terms of CV or BIC
                           criterion="CV", # CV or BIC
                           folds # if CV
                           )
{
  if (!(schoolmath::is.positive(lambdas))){
    stop("all values of lambda must be positive")
  }
  else if (!(criterion %in% c("CV", "BIC"))){
    stop("type should be one of {'CV', 'BIC'}")
  } else if (criterion == "CV"){
    if (length(fold) != 1)
    { stop("Input (fold) must be a scalar") }                      # neu
    if ((class(fold) != "numeric") & (class(fold) != "integer"))   # abgeändert
    { stop("Input (fold) is of wrong class") }
    if ((fold <=  1) | (fold > nrow(data)))
    { stop("Input (fold) out of range") } 
    if (fold %% 1 != 0)                                              # neu
    { stop("Input (fold) has to be a whole number")}
    if ((nrow(data) == 3) & (fold == 2))                                # neu
    { stop("Training set cannot consist of one observation. Input (fold) has to be 3 in that special case.")}
  }
}
# function selects optimal value for lambda (out of user-specified set of given lambda values)
# lambda_max can be determined with lambda_Upperlimit(..., type="covLASSO")


for (i in length(lambdas)){
  
}