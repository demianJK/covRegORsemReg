############################################## Funktionsdefinitionen CV Bien & Tiebshirani ##############################################

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


kCVL <- function(lambda, data, folds){ # aus rags2ridges:::.kcvl
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
  
  cvLL <- 0
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
    
    cvLL <- cvLL + LL(P, S_Ai / length(folds[[f]])) # / length(folds[[f]]) um für n_k zu kontrollieren???
  } 
  
  return(cvLL / length(folds))
}


optLambda.covLASSO <- function(data, 
                               lambdaMin, lambdaMax, 
                               lambdaInit = (lambdaMin + lambdaMax)/2,
                               fold = as.integer(nrow(data)/2) # abgerundet auf ganze Zahl
                               ){ # aus rags2ridges::optPenalty.kCVauto
  ##############################################################################
  # - Function that performs an 'automatic' search for the optimal penalty parameter 
  #   for covglasso() or spcov() call by employing Brent's method to the calculation 
  #   of a cross-validated negative log-likelihood score.
  # - data       > Data matrix. Variables assumed to be represented by columns.
  # - lambdaMin	 > A numeric giving the minimum value for the penalty parameter.
  # - lambdaMax	 > A numeric giving the maximum value for the penalty parameter.
  # - lambdaInit > A numeric giving the initial (starting) value for the penalty parameter.
  # - fold       > A numeric or integer specifying the number of folds to apply in the cross-validation.
  #
  #    --> fixed data set, variable data splittings, variable lambdas
  ##############################################################################
  # Formel aus Bien & Tiebshirani, S. 812, Abschnitt 4 alpha_cv
  
  # input checks
  if (is.data.frame(data))                      # neu
  {data <- as.matrix(data)}
  if (!is.matrix(data)) 
  { stop("Input (data) should be a matrix or a data frame") }
  if (class(lambdaMin) != "numeric")
  { stop("Input (lambdaMin) is of wrong class") }
  if (length(lambdaMin) != 1)
  { stop("Input (lambdaMin) must be a scalar") }
  if (lambdaMin <= 0)                                          # 0 <= lambda <= unendlich
  { stop("Input (lambdaMin) must be positive") }
  if (class(lambdaMax) != "numeric")
  { stop("Input (lambdaMax) is of wrong class") }
  if (length(lambdaMax) != 1)
  { stop("Input (lambdaMax) must be a scalar") }
  if (lambdaMax <= lambdaMin)
  { stop("Input (lambdaMax) must be larger than lambdaMin") }
  if (class(lambdaInit) != "numeric")
  { stop("Input (lambdaInit) is of wrong class") }
  if (length(lambdaInit) != 1)
  { stop("Input (lambdaInit) must be a scalar") }
  if (lambdaInit <= lambdaMin)
  { stop("Input (lambdaInit) must be larger than lambdaMin") }
  if (lambdaInit > lambdaMax)
  { stop("Input (lambdaInit) must be smaller than lambdaMax") }
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
  # der eine Fall, wenn training set aus nur einem Fall besteht; dann kann keine CovMat berechnet werden
  
  # make k-folds as list:
  #fold <- max(min(ceiling(fold), nrow(data)), 2) #rausgenommen
  ## aufgedröselt:
  # ceiling(fold): rundet auf nächst größere Integer, aber kann man nicht in checks oben angeben, dass ganze Zahl sein soll
  # min(ceiling(fold), nrow(data)): gibt doch immer fold aus, da fold < nrow(data) sein muss...
  # max(min(ceiling(fold), nrow(data)), 2): das fold mind. 2 sein soll ist in checks oben
  fold <- rep(1:fold, ceiling(nrow(data)/fold))[1:nrow(data)] # fold: Größe samples
  # rep(1:fold, ceiling(nrow(data)/fold))[1:nrow(data)] und rep(1:fold, ceiling(nrow(data)/fold)) gibt nur selben Output wenn kein Vektorrecycling
  shuffle <- sample(1:nrow(data), nrow(data))
  folds <- split(shuffle, as.factor(fold)) # folds: finale subsamples (Matching Zeilenindizes mit k)
  # offene Fragen:
  # wie geht Funktion damit um wenn ungleich große subsamples? ...
  # sollte man die Größe der Samples dann noch in kCVL aufnehmen? --> scheinbar nicht bedacht, siehe Formel Bien & Tiebsh. p.812
  
  # determine optimal value of penalty parameter
  optLambda <- optim(par = lambdaInit, 
                     fn = kCVL, 
                     method = "Brent", 
                     lower = lambdaMin,
                     upper = lambdaMax, 
                     data = data, # für kCVL
                     folds = folds # für kCVL
  )$par # value 'par': "best set of parameters found"
  
  # Return
  return(optLambda)
}
