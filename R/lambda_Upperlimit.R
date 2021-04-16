compute <- function(S, n, lambda_max)
  { 
  # compute covLASSO (Bien & Tiebshirani) with lambda_max
  if (is.positive.definite(S) == TRUE){
    S <- covglasso(S=S, n=n, lambda=lambda_max)$sigma
  } else { # if matrix is npd apply ridge option till it is pd
    while(any(eigen(S)$values<0)){
      S <- S + diag(1e-5, ncol = ncol(S), nrow = nrow(S))
    }
    S <- covglasso(S=S, n=n, lambda=lambda_max)$sigma
  }
}

#####

compare <- function(S, n, lambda_max, S_max)
  {
  # compute regularized matrix with lambda_max
  S <- compute(S, n, lambda_max)
  
  # compare all off-diagonal elements of lower triangle (symmetric matrix)
  comp <- S[lower.tri(S)] == S_max[lower.tri(S_max)]
  # compute relative percentage of congruence
  rel <- table(comp)[2] / (table(comp)[1] + table(comp)[2])
  rel1 <- round(rel, 1) # liberal rounding
  rel2 <- round(rel, 2) # conservative rounding (for output)
  
  return(list(rel1, rel2))
}

#####

lambda_upperLimit <- function(S, # ML estimator (/n)
                              n, # sample size
                              lambda_max, # guess for lambda max
                              stop=0.8, # stop criterion; default: 80% of off-diagonal values (covariances) shrunken to zero
                              type
                              )
{
  if (!isSymmetric(S)) {
    stop("S should be a symmetric matrix")
  }
  else if (n <= 0) {
    stop("sample size (n) should be positive")
  }
  else if (lambda_max <= 0) {
    stop("lambda_max should be positive")
  }
  else if (stop < 0 | stop > 1) {
    stop("stop criterion should be a value between 0 and 1")
  }
  else if (!(type %in% c("covLASSO", "covRidge"))){
    stop("type should be one of {'covLASSO', 'covRidge'}")
  }
  
  ########################################################
  
  
  # construct maximally regularized matrix ######################################## only for covLASSO?????
  S_max <- matrix(0, nrow=nrow(S), ncol=ncol(S))
  diag(S_max) <- diag(S)
  # compute regularized matrix and relative percentage of zero values in off-diagonal triangle
  rel <- compare(S, n, lambda_max, S_max) # list with two values
    
    if (rel$rel1 < stop){ # lamda_max too small
      return(list(rel$rel2, lambda_max))
      stop("This is the relative percentage of zero values in one off-diagonal triangle. Consider indicating a higher value for (lambda_max)")
    } else if (rel$rel1 == stop){ # lambda_max approximately fine
      return(list(rel$rel2, lambda_max))
      stop("This is the relative percentage of zero values in one off-diagonal triangle. If that's close enough to your stop criterion, then you've found your upper limit for lambda.")
    } else { # lambda_max too big
      # "pre"-loop (to determine first lower bound of lambda_max):
        while (rel$rel1 > stop){
          lambda_max_test <- lambda_max/2    # lambda_max_test = 5   lambda_max = 10  
          rel <- compare(S, n, lambda_max_test, S_max)
          # buffer variable:
          lambda_max_higher <- lambda_max      # lambda_max_higher = 10
          # for next iteration:
          lambda_max <- lambda_max_test      # lambda_max = 5     (am Ende: lambda_max = 5   lambda_test=5   lambda_higher = 10)
        }
        # actual loop:
        while (rel$rel1 != stop){
          ######################### buffer ??????? ####################
          if (rel$rel1 < stop){ # after pre-loop first (rel$rel1 < stop) is always true
            lambda_max_test <- (lambda_max_higher + lambda_max)/2    # test=7.5   higher=10    max=5
            rel <- compare(S, n, lambda_max_test, S_max)
            # buffer variables:
            # lambda_max_higher muss bleiben für (rel$rel1 < stop)
            # lambda_max_lower für (rel$rel1 > stop)
            lambda_max_lower <- lambda_max  # lower=5
            # for next iteration:
            lambda_max <- lambda_max_test   # max=7.5
          } else { # rel$rel1 > stop
            lambda_max_test <- (lambda_max_lower + lambda_max)/2   # test=6.25    lower=5    max=7.5
            rel <- compare(S, n, lambda_max_test, S_max)
            # buffer variables:
            # lambda_max_lower muss bleiben für (rel$rel1 > stop)
            # lambda_max_higher für (rel$rel1 < stop)
            lambda_max_higher <- lambda_max  # higher=7.5
            # for next iteration:
            lambda_max <- lambda_max_test    # max=6.25
          }
        }
      return(list(rel$rel2, lambda_max))
      stop("This is the relative percentage of zero values in one off-diagonal triangle. If that's close enough to your stop criterion, then you've found your upper limit for lambda.")
    }
}
      