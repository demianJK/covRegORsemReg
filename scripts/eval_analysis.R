### preparation  ########################################################################################################################################
gc()
library(xtable) # for LaTex tables
library(ggplot2)

# select conditions (currently only data with r_pred = 0 analysed)
grid <- read.csv("data/grid.csv")
selectedConds <- c() # used for all later methods
for (i in 1:nrow(grid)){
  if (grid[i, "r_pred"] == 0){
    selectedConds <- c(selectedConds, i)
  }
}

# create new grid with only analyzed conditions (r_pred = 0)
grid <- grid[selectedConds, ]

# read data
for (i in selectedConds){
  assign(paste0("analyseList_cond_", i), readRDS(paste0("ana/analyseList_cond_", i)))
}

allConds_covShrink <- readRDS("ana/covShrink/allConds_covShrink") # covShrink was re-analysed due to a wrong decision criterion in the choice of the target matrix

### Convergence Rate ####################################################################################################################################

convergence <- as.data.frame(matrix(nrow=7, ncol=18))
rownames(convergence) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
colnames(convergence) <- 1:18
allConv <- vector("list", 18)

# test 3 possbible conditions that are regarded as nonconverged
#' 1) anstatt Modell nur NA (error in comp); not needed with fiml
#' 2) min ein NA in param estims (trifft hier auf keinen Fall zu)
#' 3) model not converged

for (i in 1:18){
  cond <- selectedConds[i]
  # create temporary vectors
  ML <- c()
  FIML <- c()
  semLASSO <- c()
  semRidge <- c()
  covShrink <- c()
  covLASSO <- c()
  covRidge <- c()
  for (nrep in 1:100){
    ### TRUEs and FALSEs
    # ML
    if (is.atomic(get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]) || anyNA(get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x)){
      # if error occured while computing or at least one param estim is NA
      ML <- append(ML, FALSE)
    } else {
      ML <- append(ML, get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@converged)
    }

    # FIML
    if (anyNA(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]]) || get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["status"]][["code"]] > 1){
      # if at least one param estim is NA or if model not converged according to status code (see https://openmx.ssri.psu.edu/wiki/errors)
      FIML <- append(FIML, FALSE)
    } else {
      FIML <- append(FIML, TRUE)
    }

    # semLASSO
    if (is.atomic(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]]) || anyNA(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])){ # no param estims (weil tlw converged, aber keine param estims)
      # if error occured while computing or at least one param estim is NA
      semLASSO <- append(semLASSO, FALSE)
    } else {
      ind <- which.min(abs(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["fits"]][,"BIC"]))
      if (get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["fits"]][ind,"conv"] == 0){ # see R doc
        semLASSO <- append(semLASSO, TRUE)
      } else {
        semLASSO <- append(semLASSO, FALSE)
      }
    }

    # semRidge
    if (is.atomic(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]]) || anyNA(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])){ # no param estims (weil tlw converged, aber keine param estims)
      # if error occured while computing or at least one param estim is NA
      semRidge <- append(semRidge, FALSE)
    } else {
      ind <- which.min(abs(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["fits"]][,"BIC"]))
      if (get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["fits"]][ind,"conv"] == 0){ # see R doc
        semRidge <- append(semRidge, TRUE)
      } else {
        semRidge <- append(semRidge, FALSE)
      }
    }

    # covShrink
    if (anyNA(allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x)){
      # 100% convergence, therefore first statement superfluous, hence only checked if at least one param estim is NA
      covShrink <- append(covShrink, FALSE)
    } else {
      covShrink <- append(covShrink, allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@converged)
    }

    # covLASSO
    if (is.atomic(get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]) || anyNA(get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x)){
      # if error occured while computing or at least one param estim is NA
      covLASSO <- append(covLASSO, FALSE)
    } else {
      covLASSO <- append(covLASSO, get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@converged)
    }

    # covRidge
    if (is.atomic(get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]) || anyNA(get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x)){
      # if error occured while computing or at least one param estim is NA
      covRidge <- append(covRidge, FALSE)
    } else {
      covRidge <- append(covRidge, get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@converged)
    }
  }

  allConv[[i]] <- data.frame(ML, FIML, semLASSO, semRidge, covShrink, covLASSO, covRidge)

  ## relative convergence rate
  # ML
  if (all(ML)){ # only TRUEs
    t <- unname(table(ML)[1])
    f <- 0
  } else if (length(table(ML)) == 1){ # only FALSEs (opposite of: at least one TRUE)
    t <- 0
    f <- unname(table(ML)[1])
  } else { # TRUEs and FALSEs
    t <- unname(table(ML)[2])
    f <- unname(table(ML)[1])
  }
  convergence["ML", i] <- (t / (f + t))

  # FIML
  if (all(FIML)){ # only TRUEs
    t <- unname(table(FIML)[1])
    f <- 0
  } else if (length(table(FIML)) == 1){ # only FALSEs (opposite of: at least one TRUE)
    t <- 0
    f <- unname(table(FIML)[1])
  } else { # TRUEs and FALSEs
    t <- unname(table(FIML)[2])
    f <- unname(table(FIML)[1])
  }
  convergence["FIML", i] <- (t / (f + t))

  # semLASSO
  if (all(semLASSO)){ # only TRUEs
    t <- unname(table(semLASSO)[1])
    f <- 0
  } else if (length(table(semLASSO)) == 1){ # only FALSEs (opposite of: at least one TRUE)
    t <- 0
    f <- unname(table(semLASSO)[1])
  } else { # TRUEs and FALSEs
    t <- unname(table(semLASSO)[2])
    f <- unname(table(semLASSO)[1])
  }
  convergence["semLASSO", i] <- (t / (f + t))

  # semRidge
  if (all(semRidge)){ # only TRUEs
    t <- unname(table(semRidge)[1])
    f <- 0
  } else if (length(table(semRidge)) == 1){ # only FALSEs (opposite of: at least one TRUE)
    t <- 0
    f <- unname(table(semRidge)[1])
  } else { # TRUEs and FALSEs
    t <- unname(table(semRidge)[2])
    f <- unname(table(semRidge)[1])
  }
  convergence["semRidge", i] <- (t / (f + t))

  # covShrink
  if (all(covShrink)){ # only TRUEs
    t <- unname(table(covShrink)[1])
    f <- 0
  } else if (length(table(covShrink)) == 1){ # only FALSEs (opposite of: at least one TRUE)
    t <- 0
    f <- unname(table(covShrink)[1])
  } else { # TRUEs and FALSEs
    t <- unname(table(covShrink)[2])
    f <- unname(table(covShrink)[1])
  }
  convergence["covShrink", i] <- (t / (f + t))

  # covLASSO
  if (all(covLASSO)){ # only TRUEs
    t <- unname(table(covLASSO)[1])
    f <- 0
  } else if (length(table(covLASSO)) == 1){ # only FALSEs (opposite of: at least one TRUE)
    t <- 0
    f <- unname(table(covLASSO)[1])
  } else { # TRUEs and FALSEs
    t <- unname(table(covLASSO)[2])
    f <- unname(table(covLASSO)[1])
  }
  convergence["covLASSO", i] <- (t / (f + t))

  # covRidge
  if (all(covRidge)){ # only TRUEs
    t <- unname(table(covRidge)[1])
    f <- 0
  } else if (length(table(covRidge)) == 1){ # only FALSEs (opposite of: at least one TRUE)
    t <- 0
    f <- unname(table(covRidge)[1])
  } else { # TRUEs and FALSEs
    t <- unname(table(covRidge)[2])
    f <- unname(table(covRidge)[1])
  }
  convergence["covRidge", i] <- (t / (f + t))
}

saveRDS(convergence, file = "ana/objects/convergence")
saveRDS(allConv, file = "ana/objects/allConv")

convergence <- readRDS("ana/objects/convergence")
allConv <- readRDS("ana/objects/allConv")

#### inflated standard errors ###################################################################################################################
# note that for semReg, no se are computed

# inflatErrParams <- vector("list", length(selectedConds)) # relative numbers
# inflatErrModels <- as.data.frame(matrix(nrow=7, ncol=length(selectedConds)))
# rownames(inflatErrModels) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge") # see note above
# colnames(inflatErrModels) <- 1:18
#
# for (i in 1:length(allConv)){
#   cond <- selectedConds[i]
#   # create temporary vectors for inflatErrModels
#   ML <- 0
#   FIML <- 0
#   covShrink <- 0
#   covLASSO <- 0
#   covRidge <- 0
#     for (nrep in 1:100){
#       # ML
#       if (allConv[[i]][["ML"]][[nrep]]){ # if converged
#         infl <- get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@ParTable[["se"]] > 1000
#         if (all(infl)){ # all params inflated
#           inflatErrParams[[i]][["ML"]][[nrep]] <- 1
#           ML <- ML + 1
#         } else if (length(table(infl)) == 1){ # no params inflated
#           inflatErrParams[[i]][["ML"]][[nrep]] <- 0
#         } else { # some params inflated
#           t <- unname(table(infl)[2])
#           f <- unname(table(infl)[1])
#           inflatErrParams[[i]][["ML"]][[nrep]] <- round(t / (f + t), 5)
#           ML <- ML + 1
#         }
#       } else {
#         inflatErrParams[[i]][["ML"]][[nrep]] <- NA
#       }
#
#       # FIML
#       if (allConv[[i]][["FIML"]][[nrep]]){ # if converged
#         infl <- get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["standardErrors"]] > 1000
#         if (all(infl)){ # all params inflated
#           inflatErrParams[[i]][["FIML"]][[nrep]] <- 1
#           FIML <- FIML + 1
#         } else if (length(table(infl)) == 1){ # no params inflated
#           inflatErrParams[[i]][["FIML"]][[nrep]] <- 0
#         } else { # some params inflated
#           t <- unname(table(infl)[2])
#           f <- unname(table(infl)[1])
#           inflatErrParams[[i]][["FIML"]][[nrep]] <- round(t / (f + t), 5)
#           FIML <- FIML + 1
#         }
#       } else {
#         inflatErrParams[[i]][["FIML"]][[nrep]] <- NA
#       }
#
#       # covShrink
#       if (allConv[[i]][["covShrink"]][[nrep]]){ # if converged [note that all models are converged]
#         infl <- allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@ParTable[["se"]] > 1000
#         if (all(infl)){ # all params inflated
#           inflatErrParams[[i]][["covShrink"]][[nrep]] <- 1
#           covShrink <- covShrink + 1
#         } else if (length(table(infl)) == 1){ # no params inflated
#           inflatErrParams[[i]][["covShrink"]][[nrep]] <- 0
#         } else { # some params inflated
#           t <- unname(table(infl)[2])
#           f <- unname(table(infl)[1])
#           inflatErrParams[[i]][["covShrink"]][[nrep]] <- round(t / (f + t), 5)
#           covShrink <- covShrink + 1
#         }
#       } else {
#         inflatErrParams[[i]][["covShrink"]][[nrep]] <- NA
#       }
#
#       # covLASSO
#       if (allConv[[i]][["covLASSO"]][[nrep]]){ # if converged
#         infl <- get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@ParTable[["se"]] > 1000
#         if (all(infl)){ # all params inflated
#           inflatErrParams[[i]][["covLASSO"]][[nrep]] <- 1
#           covLASSO <- covLASSO + 1
#         } else if (length(table(infl)) == 1){ # no params inflated
#           inflatErrParams[[i]][["covLASSO"]][[nrep]] <- 0
#         } else { # some params inflated
#           t <- unname(table(infl)[2])
#           f <- unname(table(infl)[1])
#           inflatErrParams[[i]][["covLASSO"]][[nrep]] <- round(t / (f + t), 5)
#           covLASSO <- covLASSO + 1
#         }
#       } else {
#         inflatErrParams[[i]][["covLASSO"]][[nrep]] <- NA
#       }
#
#       # covRidge
#       if (allConv[[i]][["covRidge"]][[nrep]]){ # if converged
#         infl <- get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@ParTable[["se"]] > 1000
#         if (all(infl)){ # all params inflated
#           inflatErrParams[[i]][["covRidge"]][[nrep]] <- 1
#           covRidge <- covRidge + 1
#         } else if (length(table(infl)) == 1){ # no params inflated
#           inflatErrParams[[i]][["covRidge"]][[nrep]] <- 0
#         } else { # some params inflated
#           t <- unname(table(infl)[2])
#           f <- unname(table(infl)[1])
#           inflatErrParams[[i]][["covRidge"]][[nrep]] <- round(t / (f + t), 5)
#           covRidge <- covRidge + 1
#         }
#       } else {
#         inflatErrParams[[i]][["covRidge"]][[nrep]] <- NA
#       }
#     }
#
#   # if 0% convergence, then put NA there
#   if (convergence["ML", i] == 0){
#     inflatErrModels["ML", i] <- NA
#   } else {
#     inflatErrModels["ML", i] <- paste0(round(ML/(convergence["ML", i]*100), 2))
#   }
#
#   if (convergence["FIML", i] == 0){
#     inflatErrModels["FIML", i] <- NA
#   } else {
#     inflatErrModels["FIML", i] <- paste0(round(FIML/(convergence["FIML", i]*100), 2))
#   }
#
#   inflatErrModels["semLASSO", i] <- NA
#   inflatErrModels["semRidge", i] <- NA
#
#   if (convergence["covShrink", i] == 0){
#     inflatErrModels["covShrink", i] <- NA
#   } else {
#     inflatErrModels["covShrink", i] <- paste0(round(covShrink/(convergence["covShrink", i]*100), 2))
#   }
#
#   if (convergence["covLASSO", i] == 0){
#     inflatErrModels["covLASSO", i] <- NA
#   } else {
#     inflatErrModels["covLASSO", i] <- paste0(round(covLASSO/(convergence["covLASSO", i]*100), 2))
#   }
#
#   if (convergence["covRidge", i] == 0){
#     inflatErrModels["covRidge", i] <- NA
#   } else {
#     inflatErrModels["covRidge", i] <- paste0(round(covRidge/(convergence["covRidge", i]*100), 2))
#   }
# }
#
# # beautify table
# for (i in 1:7){
#   for (j in 1:18){
#     if (!is.na(inflatErrModels[i,j]) && inflatErrModels[i,j] == "0"){
#       inflatErrModels[i,j] <- "-"
#     }
#   }
# }
#
# saveRDS(inflatErrModels, file = "ana/objects/inflatErrModels")
# saveRDS(inflatErrParams, file = "ana/objects/inflatErrParams")

inflatErrModels <- readRDS("ana/objects/inflatErrModels")
inflatErrParams <- readRDS("ana/objects/inflatErrParams")


### improper solutions (negative variances) #################################################################################################################

impropParams <- vector("list", length(selectedConds)) # relative numbers
impropModels <- as.data.frame(matrix(nrow=7, ncol=length(selectedConds)))
rownames(impropModels) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
colnames(impropModels) <- 1:18

for (i in 1:length(allConv)){
  cond <- selectedConds[i]
  # create temporary vectors for impropModels
  ML <- 0
  FIML <- 0
  semLASSO <- 0
  semRidge <- 0
  covShrink <- 0
  covLASSO <- 0
  covRidge <- 0
  for (nrep in 1:100){
    # prep
    len <- length(get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x) # number of params
    # in all lavaan and resem models, last six elements are resid vars
    # ML
    if (allConv[[i]][["ML"]][[nrep]]){ # if converged
      negVar <- get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x[(len-5):len] < 0
      if (all(negVar)){ # all variances negative
        impropParams[[i]][["ML"]][[nrep]] <- 1
        ML <- ML + 1
      } else if (length(table(negVar)) == 1){ # no params inflated
        impropParams[[i]][["ML"]][[nrep]] <- 0
      } else { # some params inflated
        t <- unname(table(negVar)[2])
        f <- unname(table(negVar)[1])
        impropParams[[i]][["ML"]][[nrep]] <- round(t / (f + t), 5)
        ML <- ML + 1
      }
    } else {
      impropParams[[i]][["ML"]][[nrep]] <- NA
    }

    # FIML
    if (allConv[[i]][["FIML"]][[nrep]]){ # if converged
      # bei FIML: erst pred weights, dann factor loadings, dann resid vars (dann means)
      negVar <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(len-5):len] < 0
      if (all(negVar)){ # all params inflated
        impropParams[[i]][["FIML"]][[nrep]] <- 1
        FIML <- FIML + 1
      } else if (length(table(negVar)) == 1){ # no params inflated
        impropParams[[i]][["FIML"]][[nrep]] <- 0
      } else { # some params inflated
        t <- unname(table(negVar)[2])
        f <- unname(table(negVar)[1])
        impropParams[[i]][["FIML"]][[nrep]] <- round(t / (f + t), 5)
        FIML <- FIML + 1
      }
    } else {
      impropParams[[i]][["FIML"]][[nrep]] <- NA
    }

    # semLASSO
    if (allConv[[i]][["semLASSO"]][[nrep]]){ # if converged
      negVar <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]][(len-5):len]) < 0
      if (all(negVar)){ # all params inflated
        impropParams[[i]][["semLASSO"]][[nrep]] <- 1
        semLASSO <- semLASSO + 1
      } else if (length(table(negVar)) == 1){ # no params inflated
        impropParams[[i]][["semLASSO"]][[nrep]] <- 0
      } else { # some params inflated
        t <- unname(table(negVar)[2])
        f <- unname(table(negVar)[1])
        impropParams[[i]][["semLASSO"]][[nrep]] <- round(t / (f + t), 5)
        semLASSO <- semLASSO + 1
      }
    } else {
      impropParams[[i]][["semLASSO"]][[nrep]] <- NA
    }

    # semRidge
    if (allConv[[i]][["semRidge"]][[nrep]]){ # if converged
      negVar <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]][(len-5):len]) < 0
      if (all(negVar)){ # all params inflated
        impropParams[[i]][["semRidge"]][[nrep]] <- 1
        semRidge <- semRidge + 1
      } else if (length(table(negVar)) == 1){ # no params inflated
        impropParams[[i]][["semRidge"]][[nrep]] <- 0
      } else { # some params inflated
        t <- unname(table(negVar)[2])
        f <- unname(table(negVar)[1])
        impropParams[[i]][["semRidge"]][[nrep]] <- round(t / (f + t), 5)
        semRidge <- semRidge + 1
      }
    } else {
      impropParams[[i]][["semRidge"]][[nrep]] <- NA
    }

    # covShrink
    if (allConv[[i]][["covShrink"]][[nrep]]){ # if converged
      negVar <- allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x[(len-5):len] < 0
      if (all(negVar)){ # all params inflated
        impropParams[[i]][["covShrink"]][[nrep]] <- 1
        covShrink <- covShrink + 1
      } else if (length(table(negVar)) == 1){ # no params inflated
        impropParams[[i]][["covShrink"]][[nrep]] <- 0
      } else { # some params inflated
        t <- unname(table(negVar)[2])
        f <- unname(table(negVar)[1])
        impropParams[[i]][["covShrink"]][[nrep]] <- round(t / (f + t), 5)
        covShrink <- covShrink + 1
      }
    } else {
      impropParams[[i]][["covShrink"]][[nrep]] <- NA
    }

    # covLASSO
    if (allConv[[i]][["covLASSO"]][[nrep]]){ # if converged
      negVar <- get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x[(len-5):len] < 0
      if (all(negVar)){ # all params inflated
        impropParams[[i]][["covLASSO"]][[nrep]] <- 1
        covLASSO <- covLASSO + 1
      } else if (length(table(negVar)) == 1){ # no params inflated
        impropParams[[i]][["covLASSO"]][[nrep]] <- 0
      } else { # some params inflated
        t <- unname(table(negVar)[2])
        f <- unname(table(negVar)[1])
        impropParams[[i]][["covLASSO"]][[nrep]] <- round(t / (f + t), 5)
        covLASSO <- covLASSO + 1
      }
    } else {
      impropParams[[i]][["covLASSO"]][[nrep]] <- NA
    }


    # covRidge
    if (allConv[[i]][["covRidge"]][[nrep]]){ # if converged
      negVar <- get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x[(len-5):len] < 0
      if (all(negVar)){ # all params inflated
        impropParams[[i]][["covRidge"]][[nrep]] <- 1
        covRidge <- covRidge + 1
      } else if (length(table(negVar)) == 1){ # no params inflated
        impropParams[[i]][["covRidge"]][[nrep]] <- 0
      } else { # some params negVarated
        t <- unname(table(negVar)[2])
        f <- unname(table(negVar)[1])
        impropParams[[i]][["covRidge"]][[nrep]] <- round(t / (f + t), 5)
        covRidge <- covRidge + 1
      }
    } else {
      impropParams[[i]][["covRidge"]][[nrep]] <- NA
    }
  }

  # if 0% convergence, then put NA there
  if (convergence["ML", i] == 0){
    impropModels["ML", i] <- NA
  } else {
    impropModels["ML", i] <- paste0(round(ML/(convergence["ML", i]*100), 2))
  }
  if (convergence["FIML", i] == 0){
    impropModels["FIML", i] <- NA
  } else {
    impropModels["FIML", i] <- paste0(round(FIML/(convergence["FIML", i]*100), 2))
  }
  if (convergence["semLASSO", i] == 0){
    impropModels["semLASSO", i] <- NA
  } else {
    impropModels["semLASSO", i] <- paste0(round(semLASSO/(convergence["semLASSO", i]*100), 2))
  }
  if (convergence["semRidge", i] == 0){
    impropModels["semRidge", i] <- NA
  } else {
    impropModels["semRidge", i] <- paste0(round(semRidge/(convergence["semRidge", i]*100), 2))
  }
  if (convergence["covShrink", i] == 0){
    impropModels["covShrink", i] <- NA
  } else {
    impropModels["covShrink", i] <- paste0(round(covShrink/(convergence["covShrink", i]*100), 2))
  }
  if (convergence["covLASSO", i] == 0){
    impropModels["covLASSO", i] <- NA
  } else {
    impropModels["covLASSO", i] <- paste0(round(covLASSO/(convergence["covLASSO", i]*100), 2))
  }
  if (convergence["covRidge", i] == 0){
    impropModels["covRidge", i] <- NA
  } else {
    impropModels["covRidge", i] <- paste0(round(covRidge/(convergence["covRidge", i]*100), 2))
  }
}

# NA for conditions with zero convergence
for (i in 1:7){
  for (j in 1:18){
    if (!is.na(impropModels[i,j]) && impropModels[i,j] == "0"){
      impropModels[i,j] <- "-"
    }
  }
}

saveRDS(impropModels, file = "ana/objects/impropModels")
saveRDS(impropParams, file = "ana/objects/impropParams")

impropModels <- readRDS("ana/objects/impropModels")
impropParams <- readRDS("ana/objects/impropParams")

### final count: VALID repetitions (conv, inflatErr and negVar) #############################################################################################

## allConv, inflatErrParams, impropParams
# allUsed <- allConv
# used <- as.data.frame(matrix(nrow=7, ncol=length(selectedConds)))
# rownames(used) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
# colnames(used) <- 1:18
#
# for (i in 1:18){
#   for (nrep in 1:100){
#     # ML
#     if (allConv[[i]][["ML"]][[nrep]]){
#       if (inflatErrParams[[i]][["ML"]][[nrep]]){ # has none
#         allUsed[[i]][["ML"]][[nrep]] <- FALSE
#       }
#       if (impropParams[[i]][["ML"]][[nrep]]){
#         allUsed[[i]][["ML"]][[nrep]] <- FALSE
#       }
#     }
#
#     # FIML
#     if (allConv[[i]][["FIML"]][[nrep]]){
#       if (inflatErrParams[[i]][["FIML"]][[nrep]]){
#         allUsed[[i]][["FIML"]][[nrep]] <- FALSE
#       }
#       if (impropParams[[i]][["FIML"]][[nrep]]){
#         allUsed[[i]][["FIML"]][[nrep]] <- FALSE
#       }
#     }
#
#     # semLASSO
#     if (allConv[[i]][["semLASSO"]][[nrep]]){
#       if (impropParams[[i]][["semLASSO"]][[nrep]]){ # has none
#         allUsed[[i]][["semLASSO"]][[nrep]] <- FALSE
#       }
#     }
#
#       # semRidge
#     if (allConv[[i]][["semRidge"]][[nrep]]){
#       if (impropParams[[i]][["semRidge"]][[nrep]]){ # has none
#         allUsed[[i]][["semRidge"]][[nrep]] <- FALSE
#       }
#     }
#
#     # covShrink
#     if (allConv[[i]][["covShrink"]][[nrep]]){
#       if (inflatErrParams[[i]][["covShrink"]][[nrep]]){ # has none
#         allUsed[[i]][["covShrink"]][[nrep]] <- FALSE
#       }
#       if (impropParams[[i]][["covShrink"]][[nrep]]){ # has none
#         allUsed[[i]][["covShrink"]][[nrep]] <- FALSE
#       }
#     }
#
#     # covLASSO
#     if (allConv[[i]][["covLASSO"]][[nrep]]){
#       if (inflatErrParams[[i]][["covLASSO"]][[nrep]]){
#         allUsed[[i]][["covLASSO"]][[nrep]] <- FALSE
#       }
#       if (impropParams[[i]][["covLASSO"]][[nrep]]){
#         allUsed[[i]][["covLASSO"]][[nrep]] <- FALSE
#       }
#     }
#
#     # covRidge
#     if (allConv[[i]][["covRidge"]][[nrep]]){
#       if (inflatErrParams[[i]][["covRidge"]][[nrep]]){
#         allUsed[[i]][["covRidge"]][[nrep]] <- FALSE
#       }
#       if (impropParams[[i]][["covRidge"]][[nrep]]){
#         allUsed[[i]][["covRidge"]][[nrep]] <- FALSE
#       }
#     }
#   }
#   # ML
#   if (all(allUsed[[i]][["ML"]])){ # all TRUE
#     used["ML", i] <- 1
#   } else if (length(table(allUsed[[i]][["ML"]])) == 1){ # all FALSE
#     used["ML", i] <- 0
#   } else {
#     t <- unname(table(allUsed[[i]][["ML"]])[2])
#     f <- unname(table(allUsed[[i]][["ML"]])[1])
#     used["ML", i] <- t / (f + t)
#   }
#
#   # FIML
#   if (all(allUsed[[i]][["FIML"]])){ # all TRUE
#     used["FIML", i] <- 1
#   } else if (length(table(allUsed[[i]][["FIML"]])) == 1){ # all FALSE
#     used["FIML", i] <- 0
#   } else {
#     t <- unname(table(allUsed[[i]][["FIML"]])[2])
#     f <- unname(table(allUsed[[i]][["FIML"]])[1])
#     used["FIML", i] <- t / (f + t)
#   }
#
#   # semLASSO
#   if (all(allUsed[[i]][["semLASSO"]])){ # all TRUE
#     used["semLASSO", i] <- 1
#   } else if (length(table(allUsed[[i]][["semLASSO"]])) == 1){ # all FALSE
#     used["semLASSO", i] <- 0
#   } else {
#     t <- unname(table(allUsed[[i]][["semLASSO"]])[2])
#     f <- unname(table(allUsed[[i]][["semLASSO"]])[1])
#     used["semLASSO", i] <- t / (f + t)
#   }
#
#   # semRidge
#   if (all(allUsed[[i]][["semRidge"]])){ # all TRUE
#     used["semRidge", i] <- 1
#   } else if (length(table(allUsed[[i]][["semRidge"]])) == 1){ # all FALSE
#     used["semRidge", i] <- 0
#   } else {
#     t <- unname(table(allUsed[[i]][["semRidge"]])[2])
#     f <- unname(table(allUsed[[i]][["semRidge"]])[1])
#     used["semRidge", i] <- t / (f + t)
#   }
#
#   # covShrink
#   if (all(allUsed[[i]][["covShrink"]])){ # all TRUE
#     used["covShrink", i] <- 1
#   } else if (length(table(allUsed[[i]][["covShrink"]])) == 1){ # all FALSE
#     used["covShrink", i] <- 0
#   } else {
#     t <- unname(table(allUsed[[i]][["covShrink"]])[2])
#     f <- unname(table(allUsed[[i]][["covShrink"]])[1])
#     used["covShrink", i] <- t / (f + t)
#   }
#
#   # covLASSO
#   if (all(allUsed[[i]][["covLASSO"]])){ # all TRUE
#     used["covLASSO", i] <- 1
#   } else if (length(table(allUsed[[i]][["covLASSO"]])) == 1){ # all FALSE
#     used["covLASSO", i] <- 0
#   } else {
#     t <- unname(table(allUsed[[i]][["covLASSO"]])[2])
#     f <- unname(table(allUsed[[i]][["covLASSO"]])[1])
#     used["covLASSO", i] <- t / (f + t)
#   }
#
#   # covRidge
#   if (all(allUsed[[i]][["covRidge"]])){ # all TRUE
#     used["covRidge", i] <- 1
#   } else if (length(table(allUsed[[i]][["covRidge"]])) == 1){ # all FALSE
#     used["covRidge", i] <- 0
#   } else {
#     t <- unname(table(allUsed[[i]][["covRidge"]])[2])
#     f <- unname(table(allUsed[[i]][["covRidge"]])[1])
#     used["covRidge", i] <- t / (f + t)
#   }
# }
#
# # test if everything went fine
# for (i in 1:7){
#   for (j in 1:18){
#     if (used[i,j] <= convergence[i, j]){
#       print(TRUE)
#     }
#   }
# }
#
# saveRDS(allUsed, file = "ana/objects/allUsed")
# saveRDS(used, file = "ana/objects/used")

allUsed <- readRDS("ana/objects/allUsed")
used <- readRDS("ana/objects/used")

#xtable(used)

# ### final TABLE: conv, inflatErr and negVar #############################################################################################
# empty <- rep("", 18)
#
# valid <- rbind(empty, convergence["ML", ], inflatErrModels["ML", ], impropModels["ML", ], used["ML", ],
#                 empty, convergence["FIML", ], inflatErrModels["FIML", ], impropModels["FIML", ], used["FIML", ],
#                 empty, convergence["semLASSO", ], inflatErrModels["semLASSO", ], impropModels["semLASSO", ], used["semLASSO", ],
#                 empty, convergence["semRidge", ], inflatErrModels["semRidge", ], impropModels["semRidge", ], used["semRidge", ],
#                 empty, convergence["covShrink", ], inflatErrModels["covShrink", ], impropModels["covShrink", ], used["covShrink", ],
#                 empty, convergence["covLASSO", ], inflatErrModels["covLASSO", ], impropModels["covLASSO", ], used["covLASSO", ],
#                 empty, convergence["covRidge", ], inflatErrModels["covRidge", ], impropModels["covRidge", ], used["covRidge", ])
#
# xtable(valid)
#change rownames and other features of the table in LaTex

# saveRDS(valid, file = "ana/objects/valid")

# rm(convergence)
# rm(inflatErrModels)
# rm(inflatErrParams)
# rm(impropModels)
# rm(impropParams)
# rm(used)
# gc()

### RMSE alt ######################################################################################################################################

# note that param estim and true param are exchanged in the computations here but since the deviations are squared later on, the results are the same

# RMSE_single <- vector("list", length(selectedConds)) # for every param in every condition for every repetition
# RMSE_mean <- vector("list", length(selectedConds)) # for every parameter in every condition MEAN across repetitions
# RMSE_alt <- as.data.frame(matrix(ncol=18, nrow=7)) # table for RMSE with mean of ALL param estims in every condition
# rownames(RMSE_alt) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
# colnames(RMSE_alt) <- 1:18
#
# for (i in 1:length(allUsed)){
#   cond <- selectedConds[i]
#   # create population parameters
#   p_pred <- grid[i,"p_pred"]
#   p_pred_zero <- grid[i,"p_pred_zero"]
#   if (p_pred_zero == "40%"){
#     n0 <- 0.4
#     np <- 0.2 #c(0.4*p_pred, 0.6*p_pred, 0.8*p_pred) # each 20%
#   } else if (p_pred_zero == "70%"){
#     n0 <- 0.7
#     np <- 0.1 # c(0.7*p_pred, 0.8*p_pred, 0.9*p_pred) # each 10%
#   }
#   pars <- c(1,0.8,0.8,0.8,0.5,0.5, # 6 factor loadings
#             rep(0,n0*p_pred), # preds with 0 weights
#             rep(0.2,np*p_pred), # preds with 0.2 weight
#             rep(0.5,np*p_pred), # preds with 0.5 weight
#             rep(0.8,np*p_pred), # preds with 0.8 weight
#             rep(1,6)) # resid variances (see Jacobucci 2019 script)
#   # create temporary vectors
#   ML <- as.data.frame(matrix(ncol=100, nrow=length(pars)))
#   FIML <- as.data.frame(matrix(ncol=100, nrow=length(pars)))
#   semLASSO <- as.data.frame(matrix(ncol=100, nrow=length(pars)))
#   semRidge <- as.data.frame(matrix(ncol=100, nrow=length(pars)))
#   covShrink <- as.data.frame(matrix(ncol=100, nrow=length(pars)))
#   covLASSO <- as.data.frame(matrix(ncol=100, nrow=length(pars)))
#   covRidge <- as.data.frame(matrix(ncol=100, nrow=length(pars)))
#   for (nrep in 1:100){
#
#     if (allUsed[[i]][["ML"]][[nrep]]){ # if model converged
#       # analyseList_cond_1[[1]][["sem"]]@Fit@x has same ordering of params as pars (check with lavaan::summary())
#       ML[,nrep] <- as.vector(sqrt((t((pars - get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x)) * (pars - get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x))))
#     } else {
#       ML[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["FIML"]][[nrep]]){ # if model converged
#       # bei FIML: erst pred weights, dann factor loadings, dann resid vars (dann means)
#       fimlParams <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:length(pars)]
#       fimlParams <- c(fimlParams[(p_pred+1):(p_pred+6)], # factor loadings
#                       fimlParams[1:(p_pred)], # pred weights
#                       fimlParams[(p_pred+7):length(pars)]) # resid vars
#       FIML[,nrep] <- as.vector(sqrt((t((pars - fimlParams)) * (pars - fimlParams))))
#     } else {
#       FIML[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["semLASSO"]][[nrep]]){ # if model converged
#       # reihenfolge resem: factor loadings, reg weights, resid vars
#       semLASSO[,nrep] <- as.vector(sqrt((t((pars - unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]]))) * (pars - unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])))))
#     } else {
#       semLASSO[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["semRidge"]][[nrep]]){ # if model converged
#       semRidge[,nrep] <- as.vector(sqrt((t((pars - unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]]))) * (pars - unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])))))
#     } else {
#       semRidge[,nrep] <- NA
#     }
#
#     # all covShrink models converged
#     covShrink[,nrep] <- as.vector(sqrt((t((pars - get(paste0("analyseList_cond_", cond))[[nrep]][["covShrink"]][["sem"]]@Fit@x)) * (pars - get(paste0("analyseList_cond_", cond))[[nrep]][["covShrink"]][["sem"]]@Fit@x))))
#
#     if (allUsed[[i]][["covLASSO"]][[nrep]]){ # if model converged
#       covLASSO[,nrep] <- as.vector(sqrt((t((pars - get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x)) * (pars - get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x))))
#     } else {
#       covLASSO[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["covRidge"]][[nrep]]){ # if model converged
#       covRidge[,nrep] <- as.vector(sqrt((t((pars - get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x)) * (pars - get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x))))
#     } else {
#       covRidge[,nrep] <- NA
#     }
#   }
#   # save results
#   # single
#   RMSE_single[[i]][["ML"]] <- ML
#   RMSE_single[[i]][["FIML"]] <- FIML
#   RMSE_single[[i]][["semLASSO"]] <- semLASSO
#   RMSE_single[[i]][["semRidge"]] <- semRidge
#   RMSE_single[[i]][["covShrink"]] <- covShrink
#   RMSE_single[[i]][["covLASSO"]] <- covLASSO
#   RMSE_single[[i]][["covRidge"]] <- covRidge
#   # mean across repetitions
#   RMSE_mean[[i]][["ML"]] <- rowMeans(ML, na.rm=T)
#   RMSE_mean[[i]][["FIML"]] <- rowMeans(FIML, na.rm=T)
#   RMSE_mean[[i]][["semLASSO"]] <- rowMeans(semLASSO, na.rm=T)
#   RMSE_mean[[i]][["semRidge"]] <- rowMeans(semRidge, na.rm=T)
#   RMSE_mean[[i]][["covShrink"]] <- rowMeans(covShrink, na.rm=T)
#   RMSE_mean[[i]][["covLASSO"]] <- rowMeans(covLASSO, na.rm=T)
#   RMSE_mean[[i]][["covRidge"]] <- rowMeans(covRidge, na.rm=T)
#   # .. and mean across all params in model
#   RMSE_alt["ML", i] <- paste0(round(mean(RMSE_mean[[i]][["ML"]], na.rm=T), 2), " (", used["ML", i], ")")
#   RMSE_alt["FIML", i] <- paste0(round(mean(RMSE_mean[[i]][["FIML"]], na.rm=T), 2), " (", used["FIML", i], ")")
#   RMSE_alt["semLASSO", i] <- paste0(round(mean(RMSE_mean[[i]][["semLASSO"]], na.rm=T), 2), " (", used["semLASSO", i], ")")
#   RMSE_alt["semRidge", i] <- paste0(round(mean(RMSE_mean[[i]][["semRidge"]], na.rm=T), 2), " (", used["semRidge", i], ")")
#   RMSE_alt["covShrink", i] <- paste0(round(mean(RMSE_mean[[i]][["covShrink"]], na.rm=T), 2), " (", used["covShrink", i], ")")
#   RMSE_alt["covLASSO", i] <- paste0(round(mean(RMSE_mean[[i]][["covLASSO"]], na.rm=T), 2), " (", used["covLASSO", i], ")")
#   RMSE_alt["covRidge", i] <- paste0(round(mean(RMSE_mean[[i]][["covRidge"]], na.rm=T), 2), " (", used["covRidge", i], ")")
# }
#
# # NA for conditions with zero convergence [used]
# for (i in 1:7){
#   for (j in 1:18){
#     if (RMSE_alt[i,j] == "NaN (0)"){
#       RMSE_alt[i,j] <- ""
#     }
#   }
# }

### MSE ######################################################################################################################################

# MSE_nreps <- vector("list", length(selectedConds)) # mse for every model param in every condition for every repetition
# MSE <- as.data.frame(matrix(ncol=18, nrow=7)) # table for MSE with mean of ALL param estims in every condition
# rownames(MSE) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
# colnames(MSE) <- 1:18
#
# for (i in 1:18){
#   cond <- selectedConds[i]
#
#   # create population parameters
#   p_pred <- grid[i,"p_pred"]
#   p_pred_zero <- grid[i,"p_pred_zero"]
#   if (p_pred_zero == "40%"){
#     n0 <- 0.4
#     np <- 0.2 #c(0.4*p_pred, 0.6*p_pred, 0.8*p_pred) # each 20%
#   } else if (p_pred_zero == "70%"){
#     n0 <- 0.7
#     np <- 0.1 # c(0.7*p_pred, 0.8*p_pred, 0.9*p_pred) # each 10%
#   }
#   pars <- c(1,0.8,0.8,0.8,0.5,0.5, # 6 factor loadings
#             rep(0,n0*p_pred), # preds with 0 weights
#             rep(0.2,np*p_pred), # preds with 0.2 weight
#             rep(0.5,np*p_pred), # preds with 0.5 weight
#             rep(0.8,np*p_pred), # preds with 0.8 weight
#             rep(1,6)) # resid variances (see Jacobucci 2019 script)
#
#   # create temporary vectors
#   ML <- matrix(ncol=100, nrow=length(pars))
#   FIML <- matrix(ncol=100, nrow=length(pars))
#   semLASSO <- matrix(ncol=100, nrow=length(pars))
#   semRidge <- matrix(ncol=100, nrow=length(pars))
#   covShrink <- matrix(ncol=100, nrow=length(pars))
#   covLASSO <- matrix(ncol=100, nrow=length(pars))
#   covRidge <- matrix(ncol=100, nrow=length(pars))
#
#   for (nrep in 1:100){
#
#     if (allUsed[[i]][["ML"]][[nrep]]){ # if model converged
#       # analyseList_cond_1[[1]][["sem"]]@Fit@x has same ordering of params as pars (check with lavaan::summary())
#       ML[,nrep] <- ((get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x - pars)^2)
#     } else {
#       ML[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["FIML"]][[nrep]]){ # if model converged
#       # bei FIML: erst pred weights, dann factor loadings, dann resid vars (dann means)
#       fimlParams <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:length(pars)]
#       fimlParams <- c(fimlParams[(p_pred+1):(p_pred+6)], # factor loadings
#                       fimlParams[1:(p_pred)], # pred weights
#                       fimlParams[(p_pred+7):length(pars)]) # resid vars
#       FIML[,nrep] <- (fimlParams - pars)^2
#     } else {
#       FIML[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["semLASSO"]][[nrep]]){ # if model converged
#       # reihenfolge resem: factor loadings, reg weights, resid vars
#       semLASSO[,nrep] <- (unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]] - pars))^2
#     } else {
#       semLASSO[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["semRidge"]][[nrep]]){ # if model converged
#       semRidge[,nrep] <- (unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]] - pars))^2
#     } else {
#       semRidge[,nrep] <- NA
#     }
#
#     # all covShrink models converged
#     covShrink[,nrep] <- (allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x - pars)^2
#
#     if (allUsed[[i]][["covLASSO"]][[nrep]]){ # if model converged
#       covLASSO[,nrep] <- (get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x - pars)^2
#     } else {
#       covLASSO[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["covRidge"]][[nrep]]){ # if model converged
#       covRidge[,nrep] <- (get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x - pars)^2
#     } else {
#       covRidge[,nrep] <- NA
#     }
#   }
#   # save results
#   # MSE for every model param in every successfull repetition
#   MSE_nreps[[i]][["ML"]] <- ML
#   MSE_nreps[[i]][["FIML"]] <- FIML
#   MSE_nreps[[i]][["semLASSO"]] <- semLASSO
#   MSE_nreps[[i]][["semRidge"]] <- semRidge
#   MSE_nreps[[i]][["covShrink"]] <- covShrink
#   MSE_nreps[[i]][["covLASSO"]] <- covLASSO
#   MSE_nreps[[i]][["covRidge"]] <- covRidge
#   # .. and mean across all params in model
#   MSE["ML", i] <- paste0(round(mean(ML, na.rm=T), 2), " (", used["ML", i], ")")
#   MSE["FIML", i] <- paste0(round(mean(FIML, na.rm=T), 2), " (", used["FIML", i], ")")
#   MSE["semLASSO", i] <- paste0(round(mean(semLASSO, na.rm=T), 2), " (", used["semLASSO", i], ")")
#   MSE["semRidge", i] <- paste0(round(mean(semRidge, na.rm=T), 2), " (", used["semRidge", i], ")")
#   MSE["covShrink", i] <- paste0(round(mean(covShrink, na.rm=T), 2), " (", used["covShrink", i], ")")
#   MSE["covLASSO", i] <- paste0(round(mean(covLASSO, na.rm=T), 2), " (", used["covLASSO", i], ")")
#   MSE["covRidge", i] <- paste0(round(mean(covRidge, na.rm=T), 2), " (", used["covRidge", i], ")")
# }
#
# # NA for conditions with zero convergence [used]
# for (i in 1:7){
#   for (j in 1:18){
#     if (MSE[i,j] == "NaN (0)"){
#       MSE[i,j] <- ""
#     }
#   }
# }
#
# xtable(MSE)
#
# saveRDS(MSE_nreps, file = "ana/objects/MSE_nreps")
# saveRDS(MSE, file = "ana/objects/MSE")

MSE_nreps <- readRDS("ana/objects/MSE_nreps")
MSE <- readRDS("ana/objects/MSE")

### Bias  ####################################################################################################################

# bias_nreps <- vector("list", length(selectedConds)) # for every param in every condition for every repetition
# bias <- as.data.frame(matrix(ncol=18, nrow=7)) # table for bias with mean of ALL param estims in every condition
# rownames(bias) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
# colnames(bias) <- 1:18
#
# for (i in 1:length(allUsed)){
#   cond <- selectedConds[i]
#   # create population parameters
#   p_pred <- grid[i,"p_pred"]
#   p_pred_zero <- grid[i,"p_pred_zero"]
#   if (p_pred_zero == "40%"){
#     n0 <- 0.4
#     np <- 0.2 #c(0.4*p_pred, 0.6*p_pred, 0.8*p_pred) # each 20%
#   } else if (p_pred_zero == "70%"){
#     n0 <- 0.7
#     np <- 0.1 # c(0.7*p_pred, 0.8*p_pred, 0.9*p_pred) # each 10%
#   }
#   pars <- c(1,0.8,0.8,0.8,0.5,0.5, # 6 factor loadings
#             rep(0,n0*p_pred), # preds with 0 weights
#             rep(0.2,np*p_pred), # preds with 0.2 weight
#             rep(0.5,np*p_pred), # preds with 0.5 weight
#             rep(0.8,np*p_pred), # preds with 0.8 weight
#             rep(1,6)) # resid variances (see Jacobucci 2019 script)
#   # create temporary vectors
#   ML <- matrix(ncol=100, nrow=length(pars))
#   FIML <- matrix(ncol=100, nrow=length(pars))
#   semLASSO <- matrix(ncol=100, nrow=length(pars))
#   semRidge <- matrix(ncol=100, nrow=length(pars))
#   covShrink <- matrix(ncol=100, nrow=length(pars))
#   covLASSO <- matrix(ncol=100, nrow=length(pars))
#   covRidge <- matrix(ncol=100, nrow=length(pars))
#
#   for (nrep in 1:100){
#
#     if (allUsed[[i]][["ML"]][[nrep]]){ # if model converged
#       # analyseList_cond_1[[1]][["sem"]]@Fit@x has same ordering of params as pars (check with lavaan::summary())
#       ML[,nrep] <- get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x - pars
#     } else {
#       ML[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["FIML"]][[nrep]]){ # if model converged
#       # bei FIML: erst pred weights, dann factor loadings, dann resid vars (dann means)
#       fimlParams <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:length(pars)]
#       fimlParams <- c(fimlParams[(p_pred+1):(p_pred+6)], # factor loadings
#                       fimlParams[1:(p_pred)], # pred weights
#                       fimlParams[(p_pred+7):length(pars)]) # resid vars
#       FIML[,nrep] <- fimlParams - pars
#     } else {
#       FIML[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["semLASSO"]][[nrep]]){ # if model converged
#       # reihenfolge resem: factor loadings, reg weights, resid vars
#       semLASSO[,nrep] <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]]) - pars
#     } else {
#       semLASSO[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["semRidge"]][[nrep]]){ # if model converged
#       semRidge[,nrep] <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]]) - pars
#     } else {
#       semRidge[,nrep] <- NA
#     }
#
#     # all covShrink models converged
#     covShrink[,nrep] <- allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x - pars
#
#     if (allUsed[[i]][["covLASSO"]][[nrep]]){ # if model converged
#       covLASSO[,nrep] <- get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x - pars
#     } else {
#       covLASSO[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["covRidge"]][[nrep]]){ # if model converged
#       covRidge[,nrep] <- get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x - pars
#     } else {
#       covRidge[,nrep] <- NA
#     }
#   }
#   # save results
#   # bias for every model param in every successfull repetition
#   bias_nreps[[i]][["ML"]] <- ML
#   bias_nreps[[i]][["FIML"]] <- FIML
#   bias_nreps[[i]][["semLASSO"]] <- semLASSO
#   bias_nreps[[i]][["semRidge"]] <- semRidge
#   bias_nreps[[i]][["covShrink"]] <- covShrink
#   bias_nreps[[i]][["covLASSO"]] <- covLASSO
#   bias_nreps[[i]][["covRidge"]] <- covRidge
#   # .. and mean across all params in model
#   bias["ML", i] <- paste0(round(mean(ML, na.rm=T), 2), " (", used["ML", i], ")")
#   bias["FIML", i] <- paste0(round(mean(FIML, na.rm=T), 2), " (", used["FIML", i], ")")
#   bias["semLASSO", i] <- paste0(round(mean(semLASSO, na.rm=T), 2), " (", used["semLASSO", i], ")")
#   bias["semRidge", i] <- paste0(round(mean(semRidge, na.rm=T), 2), " (", used["semRidge", i], ")")
#   bias["covShrink", i] <- paste0(round(mean(covShrink, na.rm=T), 2), " (", used["covShrink", i], ")")
#   bias["covLASSO", i] <- paste0(round(mean(covLASSO, na.rm=T), 2), " (", used["covLASSO", i], ")")
#   bias["covRidge", i] <- paste0(round(mean(covRidge, na.rm=T), 2), " (", used["covRidge", i], ")")
# }
#
# # NA for conditions with zero convergence [used]
# for (i in 1:7){
#   for (j in 1:18){
#     if (bias[i,j] == "NaN (0)"){
#       bias[i,j] <- ""
#     }
#   }
# }
#
# xtable(bias)

# saveRDS(bias_nreps, file = "ana/objects/bias_nreps")
# saveRDS(bias, file = "ana/objects/bias")

bias_nreps <- readRDS("ana/objects/bias_nreps")
bias <- readRDS("ana/objects/bias")

### Variance  ####################################################################################################################

# var_nreps <- vector("list", length(selectedConds)) # for every parameter in every condition across repetitions
# var <- as.data.frame(matrix(ncol=18, nrow=7)) # table for variance with mean of ALL param estims in every condition
# rownames(var) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
# colnames(var) <- 1:18
#
# for (i in 1:length(allUsed)){
#
#   cond <- selectedConds[i]
#   npars <- grid[i,"p"] + 6
#   p_pred <- grid[i,"p_pred"]
#
#   # create temporary matrices
#   ML <- matrix(ncol=100, nrow=npars)
#   FIML <- matrix(ncol=100, nrow=npars)
#   semLASSO <- matrix(ncol=100, nrow=npars)
#   semRidge <- matrix(ncol=100, nrow=npars)
#   covShrink <- matrix(ncol=100, nrow=npars)
#   covLASSO <- matrix(ncol=100, nrow=npars)
#   covRidge <- matrix(ncol=100, nrow=npars)
#
#   for (nrep in 1:100){
#
#     if (allUsed[[i]][["ML"]][[nrep]]){ # if model converged
#       # analyseList_cond_1[[1]][["sem"]]@Fit@x has same ordering of params as pars (check with lavaan::summary())
#       ML[,nrep] <- get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x
#     } else {
#       ML[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["FIML"]][[nrep]]){ # if model converged
#       # bei FIML: erst pred weights, dann factor loadings, dann resid vars (dann means)
#       fimlParams <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:npars]
#       FIML[,nrep] <- c(fimlParams[(p_pred+1):(p_pred+6)], # factor loadings
#                       fimlParams[1:(p_pred)], # pred weights
#                       fimlParams[(p_pred+7):npars]) # resid vars
#     } else {
#       FIML[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["semLASSO"]][[nrep]]){ # if model converged
#       # reihenfolge resem: factor loadings, reg weights, resid vars
#       semLASSO[,nrep] <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])
#     } else {
#       semLASSO[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["semRidge"]][[nrep]]){ # if model converged
#       semRidge[,nrep] <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])
#     } else {
#       semRidge[,nrep] <- NA
#     }
#
#     # all covShrink models converged
#     covShrink[,nrep] <- allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x
#
#     if (allUsed[[i]][["covLASSO"]][[nrep]]){ # if model converged
#       covLASSO[,nrep] <- get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x
#     } else {
#       covLASSO[,nrep] <- NA
#     }
#
#     if (allUsed[[i]][["covRidge"]][[nrep]]){ # if model converged
#       covRidge[,nrep] <- get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x
#     } else {
#       covRidge[,nrep] <- NA
#     }
#   }
#   # save results
#   # bias for every model param in every successfull repetition
#
#   # ML
#   meanML <- rowMeans(ML, na.rm=T)
#   var_nreps[[i]][["ML"]] <- (ML - meanML)^2
#
#   # FIML
#   meanFIML <- rowMeans(FIML, na.rm=T)
#   var_nreps[[i]][["FIML"]] <- (FIML - meanFIML)^2
#
#   # semLASSO
#   meansemLASSO <- rowMeans(semLASSO, na.rm=T)
#   var_nreps[[i]][["semLASSO"]] <- (semLASSO - meansemLASSO)^2
#
#   # semRidge
#   meansemRidge <- rowMeans(semRidge, na.rm=T)
#   var_nreps[[i]][["semRidge"]] <- (semRidge - meansemRidge)^2
#
#   # covShrink
#   meancovShrink <- rowMeans(covShrink, na.rm=T)
#   var_nreps[[i]][["covShrink"]] <- (covShrink - meancovShrink)^2
#
#   # covLASSO
#   meancovLASSO <- rowMeans(covLASSO, na.rm=T)
#   var_nreps[[i]][["covLASSO"]] <- (covLASSO - meancovLASSO)^2
#
#   # covRidge
#   meancovRidge <- rowMeans(covRidge, na.rm=T)
#   var_nreps[[i]][["covRidge"]] <- (covRidge - meancovRidge)^2
#
#   # .. and mean across all params in model
#   var["ML", i] <- paste0(round(mean(rowSums(var_nreps[[i]][["ML"]], na.rm=T)/(used["ML", i]*100 - 1)), 2), " (", used["ML", i], ")")
#   var["FIML", i] <- paste0(round(mean(rowSums(var_nreps[[i]][["FIML"]], na.rm=T)/(used["FIML", i]*100 - 1)), 2), " (", used["FIML", i], ")")
#   var["semLASSO", i] <- paste0(round(mean(rowSums(var_nreps[[i]][["semLASSO"]], na.rm=T)/(used["semLASSO", i]*100 - 1)), 2), " (", used["semLASSO", i], ")")
#   var["semRidge", i] <- paste0(round(mean(rowSums(var_nreps[[i]][["semRidge"]], na.rm=T)/(used["semRidge", i]*100 - 1)), 2), " (", used["semRidge", i], ")")
#   var["covShrink", i] <- paste0(round(mean(rowSums(var_nreps[[i]][["covShrink"]], na.rm=T)/(used["covShrink", i]*100 - 1)), 2), " (", used["covShrink", i], ")")
#   var["covLASSO", i] <- paste0(round(mean(rowSums(var_nreps[[i]][["covLASSO"]], na.rm=T)/(used["covLASSO", i]*100 - 1)), 2), " (", used["covLASSO", i], ")")
#   var["covRidge", i] <- paste0(round(mean(rowSums(var_nreps[[i]][["covRidge"]], na.rm=T)/(used["covRidge", i]*100 - 1)), 2), " (", used["covRidge", i], ")")
# }
#
#
# # NA for conditions with zero convergence [used]
# for (i in 1:7){
#   for (j in 1:18){
#     if (var[i,j] == "0 (0)"){
#       var[i,j] <- ""
#     }
#   }
# }
#
# xtable(var)
#
# saveRDS(var_nreps, file = "ana/objects/var_nreps")
# saveRDS(var, file = "ana/objects/var")

var_nreps <- readRDS("ana/objects/var_nreps")
var <- readRDS("ana/objects/var")

### PD  ####################################################################################################################

# pdAll <- as.data.frame(matrix(ncol=2, nrow=18))
# colnames(pdAll) <- c("TRUE", "FALSE")
# for (i in 1:18){
#   cond <- selectedConds[i]
#   pd <- c()
#   for (nrep in 1:100){
#     data <- get(paste0("analyseList_cond_", cond))[[nrep]][["data"]][["dataS"]]
#     N <- nrow(data)
#     S <- ((N-1)/N)*cov(data)
#     pd <- append(pd, matrixcalc::is.positive.definite(S))
#   }
#   pdTab <- table(pd)
#   if (all(pd)){ # all true
#     pdAll[i, 1] <- 1
#     pdAll[i, 2] <- 0
#   } else if (length(pdTab) == 1){ # all false
#     pdAll[i, 1] <- 0
#     pdAll[i, 2] <- 1
#   } else {
#     pdAll[i, 1] <- unname(pdTab[2])
#     pdAll[i, 2] <- unname(pdTab[1])
#   }
# }
# 5,6, 11,12, 17,18 sind pd

### COMP TIME  ####################################################################################################################

# comptime for ALL repetitions (i.e., even if not converged) --> not used in paper
# compTime <- as.data.frame(matrix(nrow=7, ncol=18))
# rownames(compTime) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
# colnames(compTime) <- 1:18
# for (i in 1:18){
#   cond <- selectedConds[i]
#   ML <- c()
#   FIML <- c()
#   semLASSO <- c()
#   semRidge <- c()
#   covShrink <- c()
#   covLASSO <- c()
#   covRidge <- c()
#   for (nrep in 1:100){
#     ML <- append(ML, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["sem"]])
#     FIML <- append(FIML, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["fiml"]])
#     semLASSO <- append(semLASSO, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["semLASSO"]])
#     semRidge <- append(semRidge, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["semRidge"]])
#     covShrink <- append(covShrink, allConds_covShrink[[i]][[nrep]][["times"]][["all_covShrink"]])
#     covLASSO <- append(covLASSO, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["all_covLASSO"]])
#     covRidge <- append(covRidge, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["all_covRidge"]])
#   }
#   compTime["ML", i] <- round(mean(ML), 2)
#   compTime["FIML", i] <- round(mean(FIML), 2)
#   compTime["semLASSO", i] <- round(mean(semLASSO), 2)
#   compTime["semRidge", i] <- round(mean(semRidge), 2)
#   compTime["covShrink", i] <- round(mean(covShrink), 2)
#   compTime["covLASSO", i] <- round(mean(covLASSO), 2)
#   compTime["covRidge", i] <- round(mean(covRidge), 2)
# }

### comptime for all VALID repetitions
# compTimeConv <- as.data.frame(matrix(nrow=7, ncol=18))
# rownames(compTimeConv) <- c("ML", "FIML", "semLASSO", "semRidge", "covShrink", "covLASSO", "covRidge")
# colnames(compTimeConv) <- 1:18
# for (i in 1:18){
#   cond <- selectedConds[i]
#   ML <- c()
#   FIML <- c()
#   semLASSO <- c()
#   semRidge <- c()
#   covShrink <- c()
#   covLASSO <- c()
#   covRidge <- c()
#   for (nrep in 1:100){
#     if (allUsed[[i]][["ML"]][[nrep]]){ # if model converged
#       ML <- append(ML, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["sem"]])
#     } else {
#       ML <- append(ML, NA)
#     }
#     if (allUsed[[i]][["FIML"]][[nrep]]){ # if model converged
#       FIML <- append(FIML, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["fiml"]])
#     } else {
#       FIML <- append(FIML, NA)
#     }
#     if (allUsed[[i]][["semLASSO"]][[nrep]]){ # if model converged
#       semLASSO <- append(semLASSO, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["semLASSO"]])
#     } else {
#       semLASSO <- append(semLASSO, NA)
#     }
#     if (allUsed[[i]][["semRidge"]][[nrep]]){ # if model converged
#       semRidge <- append(semRidge, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["semRidge"]])
#     } else {
#       semRidge <- append(semRidge, NA)
#     }
#     if (allUsed[[i]][["covShrink"]][[nrep]]){ # if model converged
#       covShrink <- append(covShrink, allConds_covShrink[[i]][[nrep]][["times"]][["all_covShrink"]])
#     } else {
#       covShrink <- append(covShrink, NA)
#     }
#     if (allUsed[[i]][["covLASSO"]][[nrep]]){ # if model converged
#       covLASSO <- append(covLASSO, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["all_covLASSO"]])
#     } else {
#       covLASSO <- append(covLASSO, NA)
#     }
#     if (allUsed[[i]][["covRidge"]][[nrep]]){ # if model converged
#       covRidge <- append(covRidge, get(paste0("analyseList_cond_", cond))[[nrep]][["times"]][["all_covRidge"]])
#     } else {
#       covRidge <- append(covRidge, NA)
#     }
#   }
#   compTimeConv["ML", i] <- paste0(round(mean(ML, na.rm=TRUE), 2), " (", used["ML", i], ")")
#   compTimeConv["FIML", i] <- paste0(round(mean(FIML, na.rm=TRUE), 2), " (", used["FIML", i], ")")
#   compTimeConv["semLASSO", i] <- paste0(round(mean(semLASSO, na.rm=TRUE), 2), " (", used["semLASSO", i], ")")
#   compTimeConv["semRidge", i] <- paste0(round(mean(semRidge, na.rm=TRUE), 2), " (", used["semRidge", i], ")")
#   compTimeConv["covShrink", i] <- paste0(round(mean(covShrink, na.rm=TRUE), 2), " (", used["covShrink", i], ")")
#   compTimeConv["covLASSO", i] <- paste0(round(mean(covLASSO, na.rm=TRUE), 2), " (", used["covLASSO", i], ")")
#   compTimeConv["covRidge", i] <- paste0(round(mean(covRidge, na.rm=TRUE), 2), " (", used["covRidge", i], ")")
# }
#
# # NA for conditions with zero convergence [used]
# for (i in 1:7){
#   for (j in 1:18){
#     if (compTimeConv[i,j] == "NaN (0)"){
#       compTimeConv[i,j] <- ""
#     }
#   }
# }
#
# xtable(compTimeConv)
# saveRDS(compTimeConv, file = "ana/objects/compTimeConv")

compTimeConv <- readRDS("ana/objects/compTimeConv")

### TARGET MATRIX COVSHRINK  ####################################################################################################################

# target <- as.data.frame(matrix(nrow=3, ncol=18))
# rownames(target) <- c("Scaled Identity Matrix", "Identity Matrix", "Diagonal Matrix")
# colnames(target) <- 1:18
# for (i in 1:18){
#   cond <- selectedConds[i]
#   equal <- 0
#   identity <- 0
#   unequal<- 0
#   for (nrep in 1:100){
#     if (allConds_covShrink[[i]][[nrep]][["covShrink"]][["cov"]][["TargetType"]] == "equal"){
#       equal <- equal + 1
#     } else if (allConds_covShrink[[i]][[nrep]][["covShrink"]][["cov"]][["TargetType"]] == "identity") {
#       identity <- identity + 1
#     } else if (allConds_covShrink[[i]][[nrep]][["covShrink"]][["cov"]][["TargetType"]] == "unequal") {
#       unequal <- unequal + 1
#     }
#     }
#   target["Scaled Identity Matrix", i] <- equal
#   target["Identity Matrix", i] <- identity
#   target["Diagonal Matrix", i] <- unequal
# }
#
# xtable(target)
#
# saveRDS(target, file = "ana/objects/target")

target <- readRDS("ana/objects/target")

### PLOTS  ####################################################################################################################
### raw estims ###

## old rawEstims (wiht NAs) --> not used
# rawEstims <- c()
# eins <- c(1:6)
# zwei <- c(7:12)
# drei <- c(13:18)
#
# for (i in 1:18){
#   cond <- selectedConds[i]
#   if (i %in% eins){ # insert 112 - 22  = 90 NAs
#     p_pred <- 10
#     npars <- 22
#     for (nrep in 1:100){
#
#       if (allUsed[[i]][["ML"]][[nrep]]){
#         ML <- c(get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x[1:12], # factor loadings + reg weights
#                 rep(NA, 90),
#                 get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x[13:npars]) # resid vars
#       } else {
#         ML <- rep(NA, 112)
#       }
#
#       if (allUsed[[i]][["FIML"]][[nrep]]){
#         FIML <- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+1):(p_pred+6)], # factor loadings
#                   unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:(p_pred)], # pred weights
#                   rep(NA, 90),
#                   unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+7):npars]) # resid vars
#       } else {
#         FIML <- rep(NA, 112)
#       }
#
#       if (allUsed[[i]][["semLASSO"]][[nrep]]){
#         semLASSO <- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])[1:12],
#                       rep(NA, 90),
#                       unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])[13:npars])
#       } else {
#         semLASSO <- rep(NA, 112)
#       }
#
#       if (allUsed[[i]][["semRidge"]][[nrep]]){
#         semRidge <- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])[1:12],
#                       rep(NA, 90),
#                       unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])[13:npars])
#       } else {
#         semRidge <- rep(NA, 112)
#       }
#
#       covShrink <- c(allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x[1:12],
#                      rep(NA, 90),
#                      allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x[13:npars])
#
#       if (allUsed[[i]][["covLASSO"]][[nrep]]){
#         covLASSO <- c(get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x[1:12],
#                       rep(NA, 90),
#                       get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x[13:npars])
#       } else {
#         covLASSO <- rep(NA, 112)
#       }
#
#       if (allUsed[[i]][["covRidge"]][[nrep]]){
#         covRidge <- c(get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x[1:12],
#                       rep(NA, 90),
#                       get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x[13:npars])
#       } else {
#         covRidge <- rep(NA, 112)
#       }
#
#       rawEstims <- append(rawEstims, c(ML, FIML, semLASSO, semRidge, covShrink, covLASSO, covRidge))
#     }
#   }
#   if (i %in% zwei){ # insert 112 - 62 = 50 NAs
#     p_pred <- 50
#     npars <- 62
#     for (nrep in 1:100){
#
#       if (allUsed[[i]][["ML"]][[nrep]]){
#         ML <- c(get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x[1:12], # factor loadings + reg weights
#                 rep(NA, 50),
#                 get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x[13:npars]) # resid vars
#       } else {
#         ML <- rep(NA, 112)
#       }
#
#       if (allUsed[[i]][["FIML"]][[nrep]]){
#         FIML<- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+1):(p_pred+6)], # factor loadings
#                  unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:(p_pred)], # pred weights
#                  rep(NA, 50),
#                  unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+7):npars]) # resid vars
#       } else {
#         FIML <- rep(NA, 112)
#       }
#
#       if (allUsed[[i]][["semLASSO"]][[nrep]]){
#         semLASSO <- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])[1:12],
#                       rep(NA, 50),
#                       unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])[13:npars])
#       } else {
#         semLASSO <- rep(NA, 112)
#       }
#
#       if (allUsed[[i]][["semRidge"]][[nrep]]){
#         semRidge <- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])[1:12],
#                       rep(NA, 50),
#                       unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])[13:npars])
#       } else {
#         semRidge <- rep(NA, 112)
#       }
#
#       covShrink <- c(allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x[1:12],
#                      rep(NA, 50),
#                      allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x[13:npars])
#
#       if (allUsed[[i]][["covLASSO"]][[nrep]]){
#         covLASSO <- c(get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x[1:12],
#                       rep(NA, 50),
#                       get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x[13:npars])
#       } else {
#         covLASSO <- rep(NA, 112)
#       }
#
#       if (allUsed[[i]][["covRidge"]][[nrep]]){
#         covRidge <- c(get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x[1:12],
#                       rep(NA, 50),
#                       get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x[13:npars])
#       } else {
#         covRidge <- rep(NA, 112)
#       }
#
#       rawEstims <- append(rawEstims, c(ML, FIML, semLASSO, semRidge, covShrink, covLASSO, covRidge))
#     }
#   }
#   if (i %in% drei){
#     p_pred <- 100
#     npars <- 112
#     for (nrep in 1:100){
#       ML <- get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x
#       FIML<- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+1):(p_pred+6)], # factor loadings
#                unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:(p_pred)], # pred weights
#                unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+7):npars]) # resid vars
#       semLASSO <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])
#       semRidge <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])
#       covShrink <- allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x
#       covLASSO <- get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x
#       covRidge <- get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x
#
#       rawEstims <- append(rawEstims, c(ML, FIML, semLASSO, semRidge, covShrink, covLASSO, covRidge))
#     }
#   }
# }
# saveRDS(rawEstims, "ana/objects/rawEstims")
# gc()
#
# rawEstims <- readRDS("ana/objects/rawEstims") # methods hintereinander (each 112 params) in reps (100) in conds (18)


#### new param estims (without NAs)
# rawEstims <- c()
# eins <- c(1:6)
# zwei <- c(7:12)
# drei <- c(13:18)
#
# for (i in 1:18){
#   cond <- selectedConds[i]
#     if (i %in% eins){ # insert 112 - 22  = 90 NAs
#       p_pred <- 10
#       npars <- 22
#       for (nrep in 1:100){
#         # ordering of params: factor loadings, reg weights, resid vars (only fiml has to be reordered)
#         if (allUsed[[i]][["ML"]][[nrep]]){
#           ML <- get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x
#         } else {
#           ML <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["FIML"]][[nrep]]){
#           FIML <- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+1):(p_pred+6)], # factor loadings
#                    unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:(p_pred)], # pred weights
#                    unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+7):npars]) # resid vars
#         } else {
#           FIML <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["semLASSO"]][[nrep]]){
#           semLASSO <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])
#         } else {
#           semLASSO <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["semRidge"]][[nrep]]){
#           semRidge <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])
#         } else {
#           semRidge <- rep(NA, npars)
#         }
#
#         covShrink <- allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x
#
#         if (allUsed[[i]][["covLASSO"]][[nrep]]){
#           covLASSO <- get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x
#         } else {
#           covLASSO <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["covRidge"]][[nrep]]){
#           covRidge <- get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x
#         } else {
#           covRidge <- rep(NA, npars)
#         }
#
#       rawEstims <- append(rawEstims, c(ML, FIML, semLASSO, semRidge, covShrink, covLASSO, covRidge)) # 22*7*100*6 = 92400
#       }
#     }
#     if (i %in% zwei){ # insert 112 - 62 = 50 NAs
#       p_pred <- 50
#       npars <- 62
#       for (nrep in 1:100){
#
#         if (allUsed[[i]][["ML"]][[nrep]]){
#           ML <- get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x
#         } else {
#           ML <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["FIML"]][[nrep]]){
#           FIML<- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+1):(p_pred+6)], # factor loadings
#                    unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:(p_pred)], # pred weights
#                    unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+7):npars]) # resid vars
#         } else {
#           FIML <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["semLASSO"]][[nrep]]){
#           semLASSO <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])
#         } else {
#           semLASSO <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["semRidge"]][[nrep]]){
#           semRidge <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])
#         } else {
#           semRidge <- rep(NA, npars)
#         }
#
#         covShrink <- allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x
#
#         if (allUsed[[i]][["covLASSO"]][[nrep]]){
#           covLASSO <- get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x
#         } else {
#           covLASSO <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["covRidge"]][[nrep]]){
#           covRidge <- get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x
#         } else {
#           covRidge <- rep(NA, npars)
#         }
#
#       rawEstims <- append(rawEstims, c(ML, FIML, semLASSO, semRidge, covShrink, covLASSO, covRidge)) # 62*7*100*6 = 260400
#       }
#     }
#     if (i %in% drei){
#       p_pred <- 100
#       npars <- 112
#
#       for (nrep in 1:100){
#
#         if (allUsed[[i]][["ML"]][[nrep]]){
#           ML <- get(paste0("analyseList_cond_", cond))[[nrep]][["sem"]]@Fit@x
#         } else {
#           ML <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["FIML"]][[nrep]]){
#           FIML<- c(unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+1):(p_pred+6)], # factor loadings
#                    unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[1:(p_pred)], # pred weights
#                    unname(get(paste0("analyseList_cond_", cond))[[nrep]][["fiml"]]@output[["estimate"]])[(p_pred+7):npars]) # resid vars
#         } else {
#           FIML <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["semLASSO"]][[nrep]]){
#           semLASSO <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semLASSO"]][["final_pars"]])
#         } else {
#           semLASSO <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["semRidge"]][[nrep]]){
#           semRidge <- unname(get(paste0("analyseList_cond_", cond))[[nrep]][["semRidge"]][["final_pars"]])
#         } else {
#           semRidge <- rep(NA, npars)
#         }
#
#         covShrink <- allConds_covShrink[[i]][[nrep]][["covShrink"]][["sem"]]@Fit@x
#
#         if (allUsed[[i]][["covLASSO"]][[nrep]]){
#           covLASSO <- get(paste0("analyseList_cond_", cond))[[nrep]][["covLASSO"]][["sem"]]@Fit@x
#         } else {
#           covLASSO <- rep(NA, npars)
#         }
#
#         if (allUsed[[i]][["covRidge"]][[nrep]]){
#           covRidge <- get(paste0("analyseList_cond_", cond))[[nrep]][["covRidge"]][["sem"]]@Fit@x
#         } else {
#           covRidge <- rep(NA, npars)
#         }
#
#       rawEstims <- append(rawEstims, c(ML, FIML, semLASSO, semRidge, covShrink, covLASSO, covRidge)) # 112*7*100*6 = 470400
#       }
#     }
# }
# saveRDS(rawEstims, "ana/objects/rawEstims") # total rows: 823200
# gc()

rawEstims <- readRDS("ana/objects/rawEstims") # methods hintereinander (22, 62 ans 112 params, resepctively) in reps (100) in conds (18)

# Method <-factor(c(rep(c("M", "F", "sL", "sR", "cS", "cL", "cR"), each=22, times=6*100),
#                   rep(c("M", "F", "sL", "sR", "cS", "cL", "cR"), each=62, times=6*100),
#                   rep(c("M", "F", "sL", "sR", "cS", "cL", "cR"), each=112, times=6*100)),
#                 levels=c("M", "F", "sL", "sR", "cS", "cL", "cR"), ordered=T)
#
# N <- c(rep("N = 10", 22*7*2*100), # 1,2
#        rep("N = 16", 22*7*2*100), # 3,4
#        rep("N = 50", 22*7*2*100), # 5,6
#        rep("N = 50", 62*7*2*100), # 7,8
#        rep("N = 56", 62*7*2*100), # 9,10
#        rep("N = 100", 62*7*2*100), # 11,12
#        rep("N = 100", 112*7*2*100), # 13,14
#        rep("N = 106", 112*7*2*100), # 15, 16
#        rep("N = 200", 112*7*2*100)) # 17,18
# N <- factor(N, levels=c("N = 10", "N = 16", "N = 50", "N = 56", "N = 100", "N = 106", "N = 200"), ordered=T)
#
# p <- c(rep(16, 22*7*6*100), # 1,2,3,4,5,6
#        rep(56, 62*7*6*100), # 7,8,9,10,11,12
#        rep(106, 112*7*6*100)) # 13,14,15,16,17,18
#
# zeros <- c(rep(c("40%", "70%"), times=3, each=22*100*7),
#            rep(c("40%", "70%"), times=3, each=62*100*7),
#            rep(c("40%", "70%"), times=3, each=112*100*7))
# zeros <- factor(zeros, levels=c("40%", "70%"), ordered=T)
#
# cond <- c(rep(c(1:6), each=22*7*100),
#           rep(c(7:12), each=62*7*100),
#           rep(c(13:18), each=112*7*100))
# cond <- factor(cond, levels=as.character(1:18), ordered=T)
#
# rep <- c(rep(c(1:100), times=6, each=22*7), ### hier war fehler
#          rep(c(1:100), times=6, each=62*7),
#          rep(c(1:100), times=6, each=112*7))
#
# params <- c(rep(c(rep("1", 6), rep("2", 10), rep("3", 6)), 7*6*100), # 1,2,3,4,5,6 --> 22 params
#             rep(c(rep("1", 6), rep("2", 50), rep("3", 6)), 7*6*100), # 7,8,9,10,11,12 --> 62 params
#             rep(c(rep("1", 6), rep("2", 100), rep("3", 6)), 7*6*100)) # 13,14,15,16,17,18 --> 112 params
#
# paramValue <- c(rep(
#   c(rep(c(c("load 1", "load 0.8", "load 0.8", "load 0.8", "load 0.5", "load 0.5"), rep("reg 0", 4), rep("reg 0.2", 2), rep("reg 0.5", 2), rep("reg 0.8", 2), rep("resid 1", 6)), 100*7), # 40% p=16 (=112) * 100 * 7
#   rep(c(c("load 1", "load 0.8", "load 0.8", "load 0.8", "load 0.5", "load 0.5"), rep("reg 0", 7), rep("reg 0.2", 1), rep("reg 0.5", 1), rep("reg 0.8", 1), rep("resid 1", 6)), 100*7)), # 70% p=16 (=112) * 100 * 7
#                  3), # 6 conditions in total with p = 16 (raw1)
#   rep(
#     c(rep(c(c("load 1", "load 0.8", "load 0.8", "load 0.8", "load 0.5", "load 0.5"), rep("reg 0", 20), rep("reg 0.2", 10), rep("reg 0.5", 10), rep("reg 0.8", 10), rep("resid 1", 6)), 100*7), # 40% p=56 (=112) * 100 * 7
#       rep(c(c("load 1", "load 0.8", "load 0.8", "load 0.8", "load 0.5", "load 0.5"), rep("reg 0", 35), rep("reg 0.2", 5), rep("reg 0.5", 5), rep("reg 0.8", 5), rep("resid 1", 6)), 100*7)), # 70% p=56 (=112) * 100 * 7
#     3),
#   rep(
#     c(rep(c(c("load 1", "load 0.8", "load 0.8", "load 0.8", "load 0.5", "load 0.5"), rep("reg 0", 40), rep("reg 0.2", 20), rep("reg 0.5", 20), rep("reg 0.8", 20), rep("resid 1", 6)), 100*7), # 40% p=106 (=112) * 100 * 7
#       rep(c(c("load 1", "load 0.8", "load 0.8", "load 0.8", "load 0.5", "load 0.5"), rep("reg 0", 70), rep("reg 0.2", 10), rep("reg 0.5", 10), rep("reg 0.8", 10),  rep("resid 1", 6)), 100*7)), # 70% p=106 (=112) * 100 * 7
#     3))
# paramValue <- factor(paramValue, levels = c("load 1", "load 0.8", "load 0.5", "reg 0", "reg 0.2", "reg 0.5", "reg 0.8", "resid 1"), ordered=T)
#
# raw <- data.frame(Estimate=rawEstims, Method, N, p, zeros, cond, rep, params, paramValue)
#
# saveRDS(raw, "ana/objects/raw")
# gc()

raw <- readRDS("ana/objects/raw")

# auf drei facet grids aufteilen
# nrow(raw)/3
# 470400
raw1 <- raw[1:92400,] # p=16, nrow=92400
#raw1 <- raw1[!is.na(raw1$Estimate),] # remove all NA estimates
raw2 <- raw[92401:(92400+260400),] # p=56 , nrow=260400
#raw2 <- raw2[!is.na(raw2$Estimate),] # remove all NA estimates
raw3 <- raw[(92400+260400+1):823200,] # p=16 , nrow= 470400
#raw3 <- raw3[!is.na(raw3$Estimate),] # remove all NA estimates


library(ggplot2)
library(jtools) # for apa theme
library(patchwork) # for fusing two plots into one panel

# identify default colors
# library(scales)
# hue_pal()(4)

# Einzelüberschriften fur raw1, raw2, und raw3:

### raw1 ###

o1 <- ggplot(data=raw1, aes(x=Method, y=Estimate)) +
  geom_jitter(aes(col=paramValue), shape=1, size=0.5) +
  scale_color_manual(values=c("#3288bd", #load 1
                              "#66c2a5", # load 0.8
                              "#abdda4", # load 0.5
                              "#fee08b", # reg 0
                              "#fdae61", # reg 0.2
                              "#f46d43", # reg 0.5
                              "#d53e4f", # reg 0.8
                              "#000000" # resid 1
  )) +
  facet_grid(. ~ N + zeros + cond ) +
  theme_apa() + scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0)) + ggtitle(label=c("Panel A: Conditions with p = 16")) + ylab("Parameter Estimate")

# reingezoomt:

u1 <- ggplot(data=raw1, aes(x=Method, y=Estimate)) +
  geom_hline(yintercept=c(0.9975), col="#3288bd", size=1) + # load 1
  geom_hline(yintercept=c(1.0025), col="#000000") + # resid 1
  geom_hline(yintercept=c(0.7975), col="#66c2a5", size=1) + # load 0.8
  geom_hline(yintercept=c(0.8025), col="#d53e4f") + # reg 0.8
  geom_hline(yintercept=c(0.4975), col="#abdda4", size=1) + # load 0.5
  geom_hline(yintercept=c(0.5025), col="#f46d43") + # reg 0.5
  geom_hline(yintercept=c(0.2), col="#fdae61") + # reg 0.2
  geom_hline(yintercept=c(0), col="#fee08b") + # reg 0
  geom_jitter(aes(col=paramValue), #alpha=0.3,
              shape=1, size=0.5) +
  scale_color_manual(values=c("#3288bd", #load 1
                              "#66c2a5", # load 0.8
                              "#abdda4", # load 0.5
                              "#fee08b", # reg 0
                              "#fdae61", # reg 0.2
                              "#f46d43", # reg 0.5
                              "#d53e4f", # reg 0.8
                              "#000000" # resid 1
  )) +
  facet_grid(. ~ N + zeros + cond) +
  theme_apa() + scale_y_continuous(limits=c(-2, 2), expand = c(0, 0), breaks = c(-2:2)) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_blank()) + ylab("Parameter Estimate")

o1 / u1

### raw2 ###


o2 <- ggplot(data=raw2, aes(x=Method, y=Estimate)) +
  geom_jitter(aes(col=paramValue), shape=1, size=0.5) +
  scale_color_manual(values=c("#3288bd", #load 1
                              "#66c2a5", # load 0.8
                              "#abdda4", # load 0.5
                              "#fee08b", # reg 0
                              "#fdae61", # reg 0.2
                              "#f46d43", # reg 0.5
                              "#d53e4f", # reg 0.8
                              "#000000" # resid 1
  )) +
  facet_grid(. ~ N + zeros + cond ) +
  theme_apa() + scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0)) + ggtitle(label=c("Panel B: Conditions with p = 56")) + ylab("Parameter Estimate")

# reingezoomt:

u2 <- ggplot(data=raw2, aes(x=Method, y=Estimate)) +
  geom_hline(yintercept=c(0.9975), col="#3288bd", size=1) + # load 1
  geom_hline(yintercept=c(1.0025), col="#000000") + # resid 1
  geom_hline(yintercept=c(0.7975), col="#66c2a5", size=1) + # load 0.8
  geom_hline(yintercept=c(0.8025), col="#d53e4f") + # reg 0.8
  geom_hline(yintercept=c(0.4975), col="#abdda4", size=1) + # load 0.5
  geom_hline(yintercept=c(0.5025), col="#f46d43") + # reg 0.5
  geom_hline(yintercept=c(0.2), col="#fdae61") + # reg 0.2
  geom_hline(yintercept=c(0), col="#fee08b") + # reg 0
  geom_jitter(aes(col=paramValue), #alpha=0.3,
              shape=1, size=0.5) +
  scale_color_manual(values=c("#3288bd", #load 1
                       "#66c2a5", # load 0.8
                       "#abdda4", # load 0.5
                       "#fee08b", # reg 0
                       "#fdae61", # reg 0.2
                       "#f46d43", # reg 0.5
                       "#d53e4f", # reg 0.8
                       "#000000" # resid 1
                       )) +
  facet_grid(. ~ N + zeros + cond ) +
  theme_apa() + scale_y_continuous(limits=c(-2, 2), expand = c(0, 0), breaks = c(-2:2)) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_blank()) + ylab("Parameter Estimate")

o2 / u2


### raw3 ###


o3 <- ggplot(data=raw3, aes(x=Method, y=Estimate)) +
  geom_jitter(aes(col=paramValue), shape=1, size=0.5) +
  scale_color_manual(values=c("#3288bd", #load 1
                              "#66c2a5", # load 0.8
                              "#abdda4", # load 0.5
                              "#fee08b", # reg 0
                              "#fdae61", # reg 0.2
                              "#f46d43", # reg 0.5
                              "#d53e4f", # reg 0.8
                              "#000000" # resid 1
  )) +
  facet_grid(. ~ N + zeros + cond ) +
  theme_apa() + scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0)) + ggtitle(label=c("Panel C: Conditions with p = 106")) + ylab("Parameter Estimate")

# reingezoomt:

u3 <- ggplot(data=raw3, aes(x=Method, y=Estimate)) +
  geom_hline(yintercept=c(0.9975), col="#3288bd", size=1) + # load 1
  geom_hline(yintercept=c(1.0025), col="#000000") + # resid 1
  geom_hline(yintercept=c(0.7975), col="#66c2a5", size=1) + # load 0.8
  geom_hline(yintercept=c(0.8025), col="#d53e4f") + # reg 0.8
  geom_hline(yintercept=c(0.4975), col="#abdda4", size=1) + # load 0.5
  geom_hline(yintercept=c(0.5025), col="#f46d43") + # reg 0.5
  geom_hline(yintercept=c(0.2), col="#fdae61") + # reg 0.2
  geom_hline(yintercept=c(0), col="#fee08b") + # reg 0
  geom_jitter(aes(col=paramValue), #alpha=0.3,
              shape=1, size=0.5) +
  scale_color_manual(values=c("#3288bd", #load 1
                              "#66c2a5", # load 0.8
                              "#abdda4", # load 0.5
                              "#fee08b", # reg 0
                              "#fdae61", # reg 0.2
                              "#f46d43", # reg 0.5
                              "#d53e4f", # reg 0.8
                              "#000000" # resid 1
  )) +
  facet_grid(. ~ N + zeros + cond ) +
  theme_apa() + scale_y_continuous(limits=c(-2, 2), expand = c(0, 0), breaks = c(-2:2)) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_blank()) + ylab("Parameter Estimate")

o3 / u3

