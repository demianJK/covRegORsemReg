analyse <- function(dat#, # outputobjekt von simulation()
                     #srep, pred_coll, N, repl, nsReps
                    )
  {
  # auch im Original auskommentiert
  # dat[,7:106] <- scale(dat[,7:106]) # zentrieren der Prädiktoren (?)
  
  ################################################################################
  ### from multivariate normal distribution derived estimator for covariance matrix (i.e. without bessel correction)
  ################################################################################  

  n <- nrow(dat)
  covMat <- ((n-1)/n)*cov(scale(dat, scale = FALSE)) # create ML covariance matrix of centered data (necessary for unbiased MLE)
  
  ################################################################################
  ### Model Specification (lavaan)
  ################################################################################  
  
  fac_loads <- c("f1 =~ NA*y1 + y2 + y3 + y4 + y5 + y6")
  fac_varcs <- c("f1 ~~ 1*f1")
  preds <- paste("f1 ~ ", paste("x", 1:100, sep = "", collapse="+"))
  analyse.model_lavaan <- paste(fac_loads, fac_varcs, preds, sep="\n")
  
  #################################################################################
  ### Matrix Properties
  #################################################################################
  
  ## PD
  # default: tol=1e-8
  if (is.positive.definite(covMat) == TRUE) {
    pd <- TRUE
  } else {
    pd <- FALSE
  }
  
  ## condition number 
  # well-conditioned später beurteilen --> nochmal recherchieren
  cond_before <- kappa(covMat) 
  
  ###################################################################################
  ### Data Analysis
  ###################################################################################
  
  ################### RegSEM ######################################################
  # Jacobucci et al. (2016)
  # LASSO + Ridge
  # nutzt ganzen Datensatz, aber nutzt automatisch RidgeOpt wenn NPD nach Dok
  
  #### 'normales' SEM:
  
  sem <- try(lavaan::sem(analyse.model_lavaan, # user specified model
                         data=dat,
                         fixed.x=TRUE), # TRUE: exogenous ‘x’ covariates are considered fixed variables and the means, variances and covariances of these variables are fixed to their sample values
             silent=FALSE) # report of error messages not suppressed
  
  if(inherits(sem, "try-error")){
    if (pd == TRUE){
      sem <- NA
    } else {
      sem <- lavaan::sem(analyse.model_lavaan,
                         data=dat,
                         fixed.x=TRUE,
                         do.fit=FALSE) # model is not fit and current starting values of the model parameters are preserved
    }
  }
  
  ################### RegSEM Lasso ######################
  
  ####### lslx ?????? besser ####################################################################################
  
  semLASSO <- try(cv_regsem(sem, # lavaan output object (which is input of this function)
                             n.lambda=30, # number of penalization values to test ############ less for tetsing
                             jump=.01, # amount to increase penalization each iteration
                             type="lasso", # penalty type
                             pars_pen=c(7:106), # parameter indicators to penalize (here: alle RegCoeffs)
                             optMethod="rsolnp", # solver; rsolp is a nonlinear solver that doesn't rely on gradient information
                             fit.ret=c("BIC","rmsea","AIC","CAIC","EBIC.5","EBIC.25"), # which fit indices to return (6)
                             fit.ret2="train", # fit indices of which data set
                             warm.start=FALSE), # whether start values are based on previous iteration (not recommmended)
                   silent=TRUE)
  if(inherits(semLASSO,"try-error")){
    semLASSO <- NA
  }
  
  ################### RegSEM Ridge ######################
  
  semRidge <- try(cv_regsem(sem, # lavaan output object (which is input of this function)
                             n.lambda=30, # number of penalization values to test
                             jump=.01, # amount to increase penalization each iteration
                             type="ridge", # penalty type
                             pars_pen=c(7:106), # parameter indicators to penalize (here: alle RegCoeffs)
                             optMethod="rsolnp", # solver; rsolp is a nonlinear solver that doesn't rely on gradient information
                             fit.ret=c("BIC","rmsea","AIC","CAIC","EBIC.5","EBIC.25"), # which fit indices to return (6)
                             fit.ret2="train", # fit indices of which data set (?)
                             warm.start=FALSE), # whether start values are based on previous iteration (not recommmended)
                   silent=TRUE)
  if(inherits(semRidge,"try-error")){
    semLASSO <- NA
  }
  
  ################### RegcovMat ##############################################################################################################
  
  ################### RegcovMat Shrink ######################
  ### Touloumis (2015)
  ## Wahl Target Matrix (und optimal penalty)
  # !!:   # t(dat) weil functions expect variables to correspond to rows and subjects to columns
  # Routine, um Wahl Target Matrix zu automatisieren für Loop:
  dat_mt <- t(as.matrix(dat))
  input_target <- targetselection(dat_mt)
  target_lambdas <- c(input_target[1][[1]], input_target[2][[1]], input_target[3][[1]])
  
  # Eigentlich ist die Wahl der Targetmatrix noch dezidierter, aber teils schwammige Beschreibung in Paper
  # within func: cov(t(dat)) --> /n-1
  if (max(target_lambdas) == target_lambdas[1]) {
    covShrink <- shrinkcovmat.equal(dat_mt)
    target_type <- "equal"
  } else if (max(target_lambdas) == target_lambdas[2]) {
    covShrink <- shrinkcovmat.identity(dat_mt)
    target_type <- "identity"
  } else if (max(target_lambdas) == target_lambdas[3]) {
    covShrink <- shrinkcovmat.unequal(dat_mt)
    target_type <- "unequal"
  }
  
  covShrink$target_type <- target_type
  
  ############################################################## noch speichern welche Target-Matrix genutzt?
  
  # Extrahieren der regularisierten KovMat
  covShrink_sigmaHat <- covShrink$Sigmahat
  
  # well-conditioned (später beurteilen; erstmal nur condition number speichern):
  cond_Shrink <- kappa(covShrink_sigmaHat)
  
  lav.out_covShrink <- try(lavaan::sem(analyse.model_lavaan, # user specified model
                                       sample.cov = covShrink_sigmaHat,
                                       sample.nobs = N,
                                       fixed.x=TRUE), # TRUE: exogenous ‘x’ covariates are considered fixed variables and the means, variances and covariances of these variables are fixed to their sample values
                           silent=FALSE)
  
  if(inherits(lav.out_covShrink,"try-error")){   ## wenn es Fehler gab
    lav.out_covShrink <- lavaan::sem(analyse.model_lavaan,
                                     sample.cov = covShrink,
                                     sample.nobs = N,
                                     fixed.x=TRUE,
                                     do.fit=FALSE)
  }
  
  
  ################### RegcovMat Lasso ######################
  ###### Bien & Tiebshirani (2011)
  ### covglasso
  # braucht Ridge Regression when npd (S.811f)
  
  if (pd == TRUE){
    ## empirical determination of maximum value of lambda
    #lambdaMax <- 
    
    ## find optimal lambda value (numerical approach)
    lambda_covLASSO <- optLambda.covLASSO(data=dat, 
                                          lambdaMin=0.001, 
                                          lambdaMax=0.1 #lambdaMax
                                          #fold = as.integer(nrow(dat)/2)  # Default
    )
    # whichRep = 6 --> 36 warnings
    # In optimize(function(par) fn(par, ...)/con$fnscale, lower = lower,  ... :
    # NA/Inf durch größte positive Zahl ersetzt
    
    # estimation of sigma hat
    covLASSO <- covglasso(data=dat, # very similar results to S=covMAt, n=n
                          lambda=lambda_covLASSO
                          # start=diag(cov(dat)) Default (starting matrix for estimation algorithm)
                          # crit = bic  Default
                          # penalize.diag = FALSE   Default
                          )
  } else {
    ## empirical determination of maximum value of lambda
    #lambdaMax <- 
    
    ## find optimal lambda value (numerical approach)
    lambda_covLASSO <- optLambda.covLASSO(data=dat, 
                                          lambdaMin=0.001, 
                                          lambdaMax=0.1 #lambdaMax
                                          #fold = as.integer(nrow(dat)/2)  # Default
    )
    # whichRep = 6 --> 36 warnings
    # In optimize(function(par) fn(par, ...)/con$fnscale, lower = lower,  ... :
    # NA/Inf durch größte positive Zahl ersetzt
    ## apply Ridge Option to covMat (add little constant to diagonal)
    covMat_ridgeReg <- covMat # biased covmat
    const <- diag(1e-5, ncol = ncol(covMat_ridgeReg), nrow = nrow(covMat_ridgeReg))
    while(any(eigen(covMat_ridgeReg)$values < 0)){
      covMat_ridgeReg <- covMat_ridgeReg + const
    }
    # 1e-5 ist scheinbar Wert aus lavaan Ridge option 
    # https://github.com/yrosseel/lavaan/blob/master/R/lav_options.R
    # https://github.com/yrosseel/lavaan/blob/master/R/lav_samplestats_icov.R
    # alternativ:
    # while(is.positive.definite(covMat_ridgeReg) == FALSE){
    #   covMat_ridgeReg <- covMat_ridgeReg + const
    # }
    ## beide Alternativen scheinen äquivalent (weil gleiche cond numb cond_ridgeReg)
    
    # well-conditioned (später beurteilen; erstmal nur condition number speichern):
    cond_ridgeReg <- kappa(covMat_ridgeReg)
    
    covLASSO <- covglasso(S=covMat_ridgeReg,
                          n=n,
                          lambda=lambda_covLASSO
                          )
    
    covLASSO$cond_ridgeReg <- cond_ridgeReg
  } 
  
  covLASSO_sigmaHat <- covLASSO$sigma
  
  ### Modellschätzung:
  lav.out_covLASSO <- try(lavaan::sem(analyse.model_lavaan, # user specified model
                                      sample.cov = covLASSO_sigmaHat, # braucht row und colnames
                                      sample.nobs = N,
                                      fixed.x=TRUE), # TRUE: exogenous ‘x’ covariates are considered fixed variables and the means, variances and covariances of these variables are fixed to their sample values
                          silent=FALSE)
  
  if(inherits(lav.out_covLASSO,"try-error")){   ## wenn es Fehler gab
    lav.out_covLASSO <- lavaan::sem(analyse.model_lavaan,
                                    sample.cov = covLASSO,
                                    sample.nobs = N,
                                    fixed.x=TRUE,
                                    do.fit=FALSE)
  }
  
  ################### RegcovMat Ridge (Precision Matrix) ######################
  ##### van Wieringen & Wessel (2016)
  ## https://github.com/CFWP/rags2ridges/blob/master/R/rags2ridges.R
  ### rags2ridges
  #  "(5) is always p.d. when λa 2 (0; 1) ... the estimate is not necessarily well-conditioned ...
  # To obtain a well-conditioned estimate in such situations, one should choose λa not too close to zero" (S.4)
  # "ridge penalty deals with singularity and ill-conditioning through shrinkage of the eigenvalues of S^−1" (S.9)
  # hier nur covMat nutzen
  
  ### LAMBDA & covMat ESTIMATION
  # choice trough cross-validation or information criteria (see section 3.5, p 7ff)
  # here: automatic search via Brent's method  to the calculation of a cross-validated negative log-likelihood score
  # k-fold CV (Formel aus van Wieringen & Wessel, S. 8, die über Formel 9)
  covRidge <- optPenalty.kCVauto(Y=as.matrix(dat),                                
                                        lambdaMin=.001, # aus Dok
                                        lambdaMax=30, # aus Dok
                                        #lambdaInit = (lambdaMin + lambdaMax)/2,
                                        fold = as.numeric(as.integer(nrow(dat)/2)) # abgerundet auf ganze zahl,  # gleicher Wert wie bei covLASSO
                                        #cor = FALSE,
                                        #target = default.target(covML(Y)),
                                        #type = "Alt"     ## it seems NOT formula 9 (p.8) is evoked, but formulae before (unnumbered)
  )
  # output: optimal lambda and estimated precision matrix
  covRidge_sigmaHat <- solve(covRidge$optPrec) # invert precision matrix estimate
  
  ### TARGET MATRIX
  # choice trought default function: default.target(S)
  # data-driven in the sense that the input matrix S provides the information for the diagonal entries
  # lead to rotation equivariant alternative and archetypal Type I ridge estimators
  # only assumption: has to be pd, therefore covmat has to be pd
  
  ### Model ESTIMATION
  lav.out_covRidge <- try(lavaan::sem(analyse.model_lavaan, # user specified model
                                      sample.cov = covRidge_sigmaHat, #
                                      sample.nobs = N,
                                      fixed.x=TRUE), # TRUE: exogenous ‘x’ covariates are considered fixed variables and the means, variances and covariances of these variables are fixed to their sample values
                          silent=FALSE)
  
  if(inherits(lav.out_covRidge,"try-error")){   ## wenn es Fehler gab
    lav.out_covRidge <- lavaan::sem(analyse.model_lavaan, 
                                    sample.cov = covRidge, 
                                    sample.nobs = N, 
                                    fixed.x=TRUE, 
                                    do.fit=FALSE)
  }
  
  ################### FIML ######################
  ### openmx

  ### RAM Model erstellen:
  mx_model <- mxModel("model_RAM", type="RAM",
                      manifestVars = names(dat),
                     latentVars   = c("f1"),
                     mxFitFunctionML(rowDiagnostics=TRUE),
                     mxPath(from = c("f1"), to = names(dat)[1:6]), # Factor loadings
                     mxPath(from = names(dat), arrows = 2), # manifest residuals
                     mxPath(from = c("f1"), arrows = 2, free = F, values = 1), # latents fixed@1
                     mxPath(from = "one", to = names(dat), arrows = 1), # manifest means
                     mxData(dat, type = "raw")
                     )

  ### Model Estimation:
  fiml  <- mxAutoStart(mx_model) 
  fiml <- mxTryHard(mx_model)
  if (fiml$output$status$code > 1){ # dann Modell nicht konvergiert (nur konvergiert wenn 0 und 1)
    # Quelle mx codes: https://openmx.ssri.psu.edu/wiki/errors
    fiml$converged <- "no"
  }
  
  
  #####################################
  
  out <- list(
    dat,
    covMat,
    pd,
    cond_before,
    cond_ridgeReg,
    # cond_shrink = cond_shrink,
    # mehr condition numbers von anderen Verfahren?
    pred_coll,
    N,
    trueRep,
    sem, # "normales" SEM
    semLASSO, # regsem lasso
    semRidge, # regsem ridge
    covShrink, # shrinkcovmat 
    lav.out_covShrink, 
    covLASSO, # spcov
    lav.out_covLASSO,
    covRidge, # rags2ridges
    lav.out_covRidge,
    fiml # openmx
  )
  
  return(out) ########################################
}