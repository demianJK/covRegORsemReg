simulate <- function(srep, pred_coll, N, repl, nsReps){
  
  ### create population model:
  factor_loadings <- c("f1 =~ y1 + 0.8*y2 + 0.8*y3 + 0.8*y4 + 0.5*y5 + 0.5*y6")
  factor_variances <- c("f1 ~~ 1*f1") # Varianz auf 1 fixiert?
  predictors <- paste("f1 ~ ",
                      paste(paste("0*", "x", 1:70, sep = "", collapse="+"), # 70 predictors with 0 
                            paste("0.2*", "x", 71:80, sep = "", collapse="+"), # 10 predictors with 0.2 
                            paste("0.5*", "x", 81:90, sep = "",  collapse="+"), # 10 predictors with 0.5
                            paste("0.8*", "x", 91:100, sep = "",  collapse="+"), # 10 predictors with 0.8
                            sep=" + ")
  )
  variances <- paste(" x", 1:100, "~~", "1*", "x", 1:100, sep="", collapse=";") # Varianzen auf 1 fixiert???
  covariances = list()
  count=0
  for(i in 1:100){   ## Kovarianzen in der Form: x1 ~~ 0*x2
    for(j in 1:100){
      if(i != j & j > i){
        count = count+1
        covariances[count] = paste("x", i, "~~", pred_coll, "*","x", j, sep="")
      }
    }
  }
  covariances = paste(covariances, collapse=";")
  pop.model <- paste(factor_loadings, factor_variances, predictors, variances, covariances, sep="\n")
  
  gc(verbose=FALSE) # garbage collection to free some RAM; silent i.e. without printing memory statistics
  
  trueRep <- nsReps*(repl - 1) + srep # seed for simulation ################################## why his way?
  dat <- lavaan::simulateData(pop.model, sample.nobs=N, seed=trueRep)
  
  return(dat) 
}