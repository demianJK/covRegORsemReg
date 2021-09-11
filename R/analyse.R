#' analyse data with 7 different methods (covLASSO, covRidge, covShrink, semLASSO, semRidge, regular ML sem, regular FIML sem)
#'
#' @param dat data set
#'
#' @import covglasso
#' @import lavaan
#' @import matrixcalc
#' @import OpenMx
#' @import rags2ridges
#' @import regsem
#' @import ShrinkCovMat
#'
#' @return model estimations (and covariance estimations) of all 7 approaches


########### einfach alle param estims standardisieren? weil regsem benötigt standard. input (ach, das ist ja was unterschiedliches)


analyse <- function(dat){

  # sanity checks:
  if (!is.matrix(dat) & !is.data.frame(dat)) {
    stop("data must be a matrix or data frame")
  }

  ########################################################################################################################
  ### Data Preparation & Properties
  ########################################################################################################################

  N <- nrow(dat)

  # for semReg:
  datS <- scale(dat, center=TRUE, scale=TRUE) # necessary for regsem (see Jacobucci et al., 2016, p.7)
  S_biasedS <- ((N-1)/N)*cov(datS) # normal-based ML estimator of covariance matrix (i.e. without bessel correction)
  pdS <- matrixcalc::is.positive.definite(S_biasedS)
  S_ridgeS <- ridgeOption(S_biasedS)
  pd_ridgeS <- matrixcalc::is.positive.definite(S_ridgeS)

  # for covReg:
  datC <- scale(dat, center=TRUE, scale=FALSE)
  S_biasedC <- ((N-1)/N)*cov(datC) # normal-based ML estimator of covariance matrix (i.e. without bessel correction)
  pdC <- matrixcalc::is.positive.definite(S_biasedC)
  S_ridgeC <- ridgeOption(S_biasedC)
  pd_ridgeC <- matrixcalc::is.positive.definite(S_ridgeC)

  data_all <- list(data=dat,
                   dataS=datS, S_biasedS=S_biasedS, pdS=pdS, S_ridgeS=S_ridgeS, pd_ridgeS=pd_ridgeS,
                   dataC=datC, S_biasedC=S_biasedC, pdC=pdC, S_ridgeC=S_ridgeC, pd_ridgeC=pd_ridgeC)

  ########################################################################################################################
  ### Data Analysis
  ########################################################################################################################

  ########### regular SEM ################################################################################################

  ##### 1) ML ############################################################################################################

  # model specification
  fac_loads <- c("f1 =~ NA*y1 + y2 + y3 + y4 + y5 + y6") # NA for free estimation of first factor loading (default: fixed to 1)
  fac_vars <- c("f1 ~~ 1*f1") # alternatively sem(..., std.lv = TRUE); then auto.fix.first = F (opposite is default)
  preds <- paste0("f1 ~ ", paste0(names(dat)[7:ncol(dat)], collapse=" + "))
  analyse.model_lavaan <- paste(fac_loads, fac_vars, preds, sep="\n")

  # model estimation
  sem_start <- Sys.time()
  if (pdS == TRUE){
    sem  <- try(lavaan::sem(analyse.model_lavaan, # user specified model
                           data=datS,
                           fixed.x=TRUE), # TRUE: exogenous ‘x’ covariates are considered fixed variables and the means, variances and covariances of these variables are fixed to their sample values (and are not estimated)
    silent=TRUE) # report of error messages not suppressed

    if(inherits(sem, "try-error")){
      sem<- try(lavaan::sem(analyse.model_lavaan,
                         data=datS,
                         fixed.x=TRUE,
                         do.fit=FALSE), # model is not fit and current starting values of the model parameters are preserved
                silent=TRUE)
    }

    if(inherits(sem,"try-error")){
      sem <- NA
    }

  } else { # when npd then supply through ridge option attained pd covmat instead of data (since regsem extracts data or covmat from sem output; but ridge option covmat not in sem output)
    sem <- try(lavaan::sem(analyse.model_lavaan, # user specified model
                           sample.cov = S_ridgeS,
                           sample.nobs = N,
                           sample.cov.rescale=FALSE,
                           missing = "listwise",
                           fixed.x=TRUE), # TRUE: exogenous ‘x’ covariates are considered fixed variables and the means, variances and covariances of these variables are fixed to their sample values (and are not estimated)
    silent=TRUE) # report of error messages not suppressed

    if(inherits(sem, "try-error")){
      sem <- try(lavaan::sem(analyse.model_lavaan,
                         sample.cov = S_ridgeS,
                         sample.nobs = N,
                         sample.cov.rescale=FALSE,
                         fixed.x=TRUE,
                         do.fit=FALSE), # model is not fit and current starting values of the model parameters are preserved
                 silent=TRUE)
    }

    if(inherits(sem,"try-error")){
      sem <- NA
    }

  }
  sem_end <- Sys.time()
  sem_time <- as.numeric(difftime(sem_end, sem_start, units = "secs"))


  ###### 2) FIML ############################################################################################################

  # model specification
  mx_model <- OpenMx::mxModel(model="model_RAM", type="RAM",
                              manifestVars = names(dat), latentVars = c("f1"),
                              #mxFitFunctionML(rowDiagnostics=TRUE),
                              mxData(observed=datS, type = "raw"), # raw and rowDiagnostics --> FIML
                              mxPath(from = c("f1"), to = names(dat)[1:6]#, values = rep(1, 6)
                                     ), # Factor loadings
                              mxPath(from = names(dat), arrows = 2, #values = rep(1, ncol(dat))
                                     ), # residuals (indicators and predictors)
                              mxPath(from = c("f1"), arrows = 2, free = F, values = 1), # fix latent var to one
                              mxPath(from = "one", to = names(dat)), # manifest means
                              mxPath(from = "one", to = "f1", free=F, values=c(0)), # fix latent mean to zero
                              mxPath(from = names(dat)[7:ncol(dat)], to = c("f1")) # predictors
  )
  # free=T (default) und values=... dann frei geschätzt mit starting values
  # wenn free=F dann fixiert auf diesem Wert

  # model estimation
  fiml_start <- Sys.time()
  fiml <- OpenMx::mxRun(mx_model)
  fiml_end <- Sys.time()
  fiml_time <- as.numeric(difftime(fiml_end, fiml_start, units = "secs"))

  # if (fiml$output$status$code > 1){ # if not converged use tryhard; for status codes see https://openmx.ssri.psu.edu/wiki/errors
  #   fimlTH_start <- Sys.time()
  #   fiml <- OpenMx::mxTryHard(mx_model)
  #   fimlTH_end <- Sys.time()
  #   fimlTH_time <- as.numeric(difftime(fimlTH_end, fimlTH_start, units = "secs"))
  # } else {
  #   fimlTH_time <- 0
  # }
  # all_fiml_time <- fiml_time + fimlTH_time


  #fiml  <- OpenMx::mxAutoStart(mx_model) # nicht nehmen, weil Model zu ULS oder DWLS geändert wird!
  #fiml <- OpenMx::mxTryHard(mx_model, exhaustive = T)
  # if (fiml$output$status$code > 1){ # dann Modell nicht konvergiert (nur konvergiert wenn 0 und 1)
  #   fiml$convergence <- FALSE
  # }

  # im vergleich zu sem: residuen x1 bis xn geschätzt und mittelwerte indikatoren (aber sonst sehr ähnliche schätzungen)


  ########## semReg ##########################################################################################################

  ##### 1) semLASSO #######################################################################################################

  semLASSO_start <- Sys.time()
  semLASSO <- try(regsem::cv_regsem(sem,
                                    n.lambda = 30,
                                    jump = .01,
                                    #lambda.start=0,
                                    type = "lasso",
                                    pars_pen = "regressions",
                                    fit.ret = c("BIC", "rmsea", "AIC", "CAIC", "EBIC.5", "EBIC.25"),
                                         ),
                       silent=TRUE)

  if(inherits(semLASSO,"try-error")){
    semLASSO <- NA
  }

  semLASSO_end <- Sys.time()
  semLASSO_time <- as.numeric(difftime(semLASSO_end, semLASSO_start, units = "secs"))

  ##### 2) semRidge #######################################################################################################

  semRidge_start <- Sys.time()
  semRidge <- try(regsem::cv_regsem(sem,
                                    n.lambda = 30,
                                    jump = .01,
                                    #lambda.start=0,
                                    type = "ridge",
                                    pars_pen = "regressions",
                                    fit.ret = c("BIC", "rmsea", "AIC", "CAIC", "EBIC.5", "EBIC.25"),
                                    ),
                       silent=TRUE)

  if(inherits(semRidge,"try-error")){
    semRidge <- NA
  }

  semRidge_end <- Sys.time()
  semRidge_time <- as.numeric(difftime(semRidge_end, semRidge_start, units = "secs"))

  ########## covReg ##############################################################################################################
  ##### 1) covShrink ##########################################################################################################
  # uses raw data

  # 1) covariance matrix estimation
  covShrink_start <- Sys.time()
  datCT <- t(datC) # parameter "data" expects rows of the data matrix data correspond to variables and the columns to subjects
  targetSelect <- ShrinkCovMat::targetselection(datCT, centered=TRUE)

  # whether lambdas differ significantly
  if (round(targetSelect[[1]], 1) != round(targetSelect[[2]], 1) || round(targetSelect[[1]], 1) != round(targetSelect[[3]], 1) || round(targetSelect[[2]], 1) != round(targetSelect[[3]], 1)){
    # if optimal shrinkage intensities differ significantly, then chose the one with the max shrinkage intensity
    optShrinks <- c(targetSelect[[1]], targetSelect[[2]], targetSelect[[3]])
    maxShrink <- max(optShrinks)
    if (maxShrink == 1){
      covShrink <- ShrinkCovMat::shrinkcovmat.equal(data=datCT, centered=TRUE)
      covShrink$TargetType <- "equal"
    } else if (maxShrink == 2){
      covShrink <- ShrinkCovMat::shrinkcovmat.unequal(data=datCT, centered=TRUE)
      covShrink$TargetType <- "unequal"
    } else {
      covShrink <- ShrinkCovMat::shrinkcovmat.equal(data=datCT, centered=TRUE)
      covShrink$TargetType <- "equal"
    }
    # whether average sample cov is close to 1 (i.e., < 0.96)
  } else if (!(round(targetSelect[[5]], 1) < 1)){
    covShrink <- ShrinkCovMat::shrinkcovmat.identity(data=datCT, centered=TRUE)
    covShrink$TargetType <- "identity"
    # whether r is large (say more than one unit so as to account for the sampling variability)
  } else if (targetSelect[[4]] >= 1){
    covShrink <- ShrinkCovMat::shrinkcovmat.unequal(data=datCT, centered=TRUE)
    covShrink$TargetType <- "unequal"
    # otherwise use vIp
  } else {
    covShrink <- ShrinkCovMat::shrinkcovmat.equal(data=datCT, centered=TRUE)
    covShrink$TargetType <- "equal"
  } # see p. 258 & 259 for decision criteria
  # output: Sigmahat (Tippfehler in Dok!), lambdahat, sigmasample, Target

  covShrink_end <- Sys.time()
  covShrink_time <- as.numeric(difftime(covShrink_end, covShrink_start, units = "secs"))

  rownames(covShrink$Sigmahat) <- colnames(dat)
  colnames(covShrink$Sigmahat) <- colnames(dat)

  # 2) model estimation
  sem_covShrink_start <- Sys.time()
  sem_covShrink <- try(lavaan::sem(analyse.model_lavaan,
                                       sample.cov = covShrink$Sigmahat,
                                       sample.cov.rescale = FALSE, # sample.cov would be rescaled with (n/n-1)
                                       sample.nobs = N,
                                       fixed.x=TRUE),
                           silent=TRUE)

  if(inherits(sem_covShrink,"try-error")){
    sem_covShrink <- try(lavaan::sem(analyse.model_lavaan,
                                     sample.cov = covShrink$Sigmahat,
                                     sample.cov.rescale = FALSE,
                                     sample.nobs = N,
                                     fixed.x=TRUE,
                                     do.fit=FALSE),
                         silent=TRUE)
  }

  if(inherits(sem_covShrink,"try-error")){
    sem_covShrink <- NA
  }

  sem_covShrink_end <- Sys.time()
  sem_covShrink_time <- as.numeric(difftime(sem_covShrink_end, sem_covShrink_start, units = "secs"))
  all_covShrink_time <- covShrink_time + sem_covShrink_time

  ################### covLASSO ######################
  # needs centered data (see p. 808)

  # 1) covariance matrix estimation
  covLASSO_start <- Sys.time()
  # find maximum value of lambda
  covLASSO_lambdaMax <- lambdaUpperLimit(data=dat, lambdaMax=30, # value from documentation optPenalty.kCVauto() (i.e. covRidge)
                                       type="covLASSO") # stop: when 80% of off-diagonal values are zero (default)

  # to find best lambda value vector of values is supplied
  covLASSO_lambdaVector <- seq(0, covLASSO_lambdaMax$lambdaMax, # empirically determined
                              .005) # step size

  if (pdC == TRUE){
    covLASSO <- covglasso::covglasso(data=datC,
                                         lambda=covLASSO_lambdaVector,
                                         crit="bic" # model selection criterion
                                         # start=diag(diag(S)) (starting matrix for estimation algorithm)
                                         # penalize.diag = FALSE
                                         # regularize all off-diagonal elements (default: diagonale of all ones except diagonal)
                                         )
  } else { # if matrix is npd, ridge option is employed
    covLASSO <- covglasso::covglasso(S=S_ridgeC,
                                         n=N,
                                         lambda=covLASSO_lambdaVector,
                                         crit="bic")
  } # output: sigma, ..., log-likelihood, ..., bic, BIC, ... lambda
  covLASSO_end <- Sys.time()
  covLASSO_time <- as.numeric(difftime(covLASSO_end, covLASSO_start, units = "secs"))
  # add best lambda value to result list (bc only sigma hat with best lambda saved in output):
  covLASSO$bestLambda <- covLASSO_lambdaVector[which(covLASSO$BIC == covLASSO$bic)]

  # 2) Model estimation
  sem_covLASSO_start <- Sys.time()
  sem_covLASSO <- try(lavaan::sem(analyse.model_lavaan,
                                          sample.cov = covLASSO$sigma,
                                          sample.cov.rescale = FALSE,
                                          sample.nobs = N,
                                          fixed.x=TRUE),
                              silent=TRUE)

  if(inherits(sem_covLASSO,"try-error")){
    sem_covLASSO <- try(lavaan::sem(analyse.model_lavaan,
                                        sample.cov = covLASSO$sigma,
                                        sample.cov.rescale = FALSE,
                                        sample.nobs = N,
                                        fixed.x=TRUE,
                                        do.fit=FALSE),
                        silent=TRUE)
  }

  if(inherits(sem_covLASSO,"try-error")){
    sem_covLASSO <- NA
  }

  # sem_covLASSO <- tryCatch(lavaan::sem(analyse.model_lavaan,
  #                                 sample.cov = covLASSO$sigma,
  #                                 sample.cov.rescale = FALSE,
  #                                 sample.nobs = N,
  #                                 fixed.x=TRUE),
  #                     finally = tryCatch(lavaan::sem(analyse.model_lavaan,
  #                                                    sample.cov = covLASSO$sigma,
  #                                                    sample.cov.rescale = FALSE,
  #                                                    sample.nobs = N,
  #                                                    fixed.x=TRUE,
  #                                                    do.fit=FALSE),
  #                                        finally = NA))

  sem_covLASSO_end <- Sys.time()
  sem_covLASSO_time <- as.numeric(difftime(sem_covLASSO_end, sem_covLASSO_start, units = "secs"))
  all_covLASSO_time <- covLASSO_time + sem_covLASSO_time

  ##### 3) covRidge ########################################################################################################
  # needs centered data (see p. 1)

  # 1) covariance matrix estimation
  covRidge_start <- Sys.time()
  # determine maximal value for lambda
  covRidge_lambdaMax <- lambdaUpperLimit(data=datC, lambdaMax=30, # value from documentation optPenalty.kCVauto()
                                       type="covRidge") # stop: when 80% of regularized matrix elements equal target matrix (liberal rounding first number after decimal point)
  # numerical approximation of optimal lambda
  covRidge <- rags2ridges::optPenalty.kCVauto(Y=datC, # cross-prod approach via rags2ridges::covML()
                                                  lambdaMin=0.001, # value from documentation optPenalty.kCVauto(); also lambdaMax
                                                  lambdaMax=covRidge_lambdaMax$lambdaMax,
                                                  #lambdaInit = (lambdaMin + lambdaMax)/2,
                                                  fold = as.numeric(as.integer(N/2)) # rounded on whole number; same choice as with covLASSO
                                                  #target = default.target(covML(Y)),
                                                  #type = "Alt"
                                                  )  # output: optimal lambda and estimated precision matrix
  covRidge$sigma <- solve(covRidge$optPrec) # invert precision matrix estimate
  covRidge_end <- Sys.time()
  covRidge_time <- as.numeric(difftime(covRidge_end, covRidge_start, units = "secs"))

  # 2) Model estimation
  sem_covRidge_start <- Sys.time()
  sem_covRidge <- try(lavaan::sem(analyse.model_lavaan,
                                          sample.cov = covRidge$sigma,
                                          sample.cov.rescale = FALSE,
                                          sample.nobs = N,
                                          fixed.x=TRUE),
                              silent=TRUE)

  if(inherits(sem_covRidge,"try-error")){
    sem_covRidge <- try(lavaan::sem(analyse.model_lavaan,
                                        sample.cov = covRidge$sigma,
                                        sample.cov.rescale = FALSE,
                                        sample.nobs = N,
                                        fixed.x=TRUE,
                                        do.fit=FALSE),
                        silent=TRUE)
  }

  if(inherits(sem_covRidge,"try-error")){
    sem_covRidge <- NA
  }

  sem_covRidge_end <- Sys.time()
  sem_covRidge_time <- as.numeric(difftime(sem_covRidge_end, sem_covRidge_start, units = "secs"))
  all_covRidge_time <- covRidge_time + sem_covRidge_time


  #####################################

  out <- list(
    data=data_all,
    sem=sem,
    fiml=fiml,
    semLASSO=semLASSO,
    semRidge=semRidge,
    covShrink=list(cov=covShrink, sem=sem_covShrink),
    covLASSO=list(cov=covLASSO, sem=sem_covLASSO),
    covRidge=list(cov=covRidge, sem=sem_covRidge),
    times=list(sem=sem_time, fiml=fiml_time,
               semLASSO=semLASSO_time, semRidge=semRidge_time,
               covShrink=covShrink_time, sem_covShrink=sem_covShrink_time, all_covShrink=all_covShrink_time,
               covLASSO=covLASSO_time, sem_covLASSO=sem_covLASSO_time, all_covLASSO=all_covLASSO_time,
               covRidge=covRidge_time, sem_covRidge=sem_covRidge_time, all_covRidge=all_covRidge_time)
  )

  return(out) ########################################
}
