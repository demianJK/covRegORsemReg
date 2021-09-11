### data analysis
# run script in multiple R sessions (i.e., no parallel computation package used)

rm(list=ls()) # clear workspace

## for users:
# if(!require(devtools))install.packages("devtools")
# devtools::install_github("demianJK/covRegORsemREG")
# library(covRegORsemReg)

## for me:
# when Project opened
devtools::load_all() # load covRegORsemReg

cond <- 33 # MANUELL festlegen (1:36)
nreps <- 100

# load data
dataList <- readRDS(paste0("data/dataList_cond_", cond))

# initialize output lists
analyseList <- vector(mode = "list", length = nreps)

for (rep in 1:nreps){ # number of repetitions
  analyseList[[rep]] <- readRDS(paste0("ana/single/cond_", cond, "/rep_", rep))
  dat <- dataList[[rep]]
  analyseList[[rep]] <- analyse(dat)
  saveRDS(analyseList[[rep]], file = paste0("ana/single/cond_", cond, "/rep_", rep)) # sanve single reps
}

if (rep == nreps){
  saveRDS(analyseList, file = paste0("ana/analyseList_cond_", cond)) # save all reps
}

