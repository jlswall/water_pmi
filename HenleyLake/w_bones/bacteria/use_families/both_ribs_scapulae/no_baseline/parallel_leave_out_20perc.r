library("tidyverse")
library("randomForest")
# library("parallel")

library("future.apply")
plan(multisession, workers=10)


# ##################################################
# Are we dealing with phlya, orders, or families?
taxalevel <- "families"

# Read in combined rib/scapula dataset (families taxa).
allT <- read_csv("../bones_combine_rib_scapula_massaged.csv")

# We're not using the baseline observations (degdays==0), so we remove them.
allT <- allT %>% filter(degdays > 0)
# ##################################################



# ##################################################
# Put the data in wide format; remove days, subj, and rare taxa.

# Move back to wide format.
wideT <- allT %>%
  filter(taxon!="Rare") %>%
  select(degdays, sampleName, taxon, fracBySample) %>%
  spread(taxon, fracBySample) %>%
  select(-sampleName)

rm(allT)
# ##################################################



# ##################################################
# Try random forests for regression using "degdays" as the response
# variable.

# #########
# How many predictors?  (All columns except response: "degdays").
numPredictors <- ncol(wideT) - 1

# Try different numbers of bootstrap samples.
numBtSampsVec <- c(600, 1500, 2100, 3000)

# Try different values for mtry (which represents how many variables
# can be chosen from at each split of the tree).
numVarSplitVec <- seq(2, numPredictors, by=1)

# Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)


# ###########################
# Do cross-validation over and over, leaving out a different 20% of
# the observations each time.

# Number of times to do cross-validation.
numCVs <- 1000
# How many observations to reserve for testing each time.
numLeaveOut <- round(0.20 * nrow(wideT))

# For matrix to hold cross-validation results.
cvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
cvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)
# origUnitsqrtcvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
# origUnitsqrtcvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)
# sqrtcvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
# sqrtcvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)



# #########################################
# Set up function for fitting random forest model using original
# units.
origUnitsF <- function(x, jCombo){
  rf <- randomForest(degdays ~ . , data=x$trainT,
    mtry=combos[jCombo, "numVarSplit"], ntree=combos[jCombo, "numBtSamps"],
    importance=T)
  return(predict(rf, newdata=x$validT))
}

# Set up function for fitting random forest model using square root
# units.
# sqrtUnitsF <- function(x, jCombo){
#   sqrtrf <- randomForest(sqrt(degdays) ~ . , data=x$trainT,
#     mtry=combos[jCombo, "numVarSplit"], ntree=combos[jCombo, "numBtSamps"],
#     importance=T)
#   return(predict(sqrtrf, newdata=x$validT))
# }
# #########################################


# #########################################
# Get set up for cross-validation.
crossvalidL <- vector("list", numCVs)
set.seed(6495244)
for (i in 1:numCVs){
  lvOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)
  trainT <- wideT[-lvOut,]
  validT <- wideT[lvOut,]
  crossvalidL[[i]] <- list(trainT=trainT, validT=validT)
}
rm(i, lvOut, trainT, validT)


# Try using lapply to fit the random forests.
origFitL <- vector("list", nrow(combos))
for (j in 1:nrow(combos)){
  # origFitL[[j]] <- mclapply(crossvalidL, mc.cores=4, origUnitsF, jCombo=j)
    origFitL[[j]] <- future_lapply(crossvalidL, FUN=origUnitsF, jCombo=j,
      future.seed=as.integer(9347622))
  if (j %% 2 == 0)
    print(paste0("In orig units, finished combo number ", j))
}    
rm(j)

# sqrtFitL <- vector("list", nrow(combos))
# for (j in 1:nrow(combos)){
#   sqrtFitL[[j]] <- mclapply(crossvalidL, mc.cores=4, sqrtUnitsF, jCombo=j)
#   if (j %% 2 == 0)
#     print(paste0("In sqrt units, finished combo number ", j))
# }
# rm(j)
# #########################################


# #########################################
# Now, calculate the various summary statistics for each cross-validation.

for (i in 1:numCVs){

  # Get the validation set for this run from the list.
  validT <- crossvalidL[[i]][["validT"]]

  # Calculate SSTotal for the cross-validation set.  If the validation
  # observations are all from the same day, then SSTot is 0.  This will only
  # happen if validations subset is at least as small as the number of
  # observations per day.
  SSTot <- sum( (validT$degdays-mean(validT$degdays))^2 )
  if (SSTot == 0)
    warning(paste("In validation subset", i, "SSTot is 0.  All degdays are",
      unique(validT$degdays), sep=" "))
  
  for (j in 1:nrow(combos)){
     
    # Calculate the MSE and error fraction of the SS Total for the
    # validation data in the original units.
    resid <- origFitL[[j]][[i]] - validT$degdays
    cvMSE[j,i] <- mean(resid^2)
    # If the validation observations are all from the same days, then SSTot is
    # 0. Avoid divide by 0.
    if (SSTot > 0)
      cvErrFrac[j,i] <- sum(resid^2)/SSTot
    else
      cvErrFrac[j,i] <- NA
    rm(resid)
  
    # Calculate the MSE and error fraction of the SS Total for the
    # validation data in the original units.
    # sqrtUnitResid <- sqrtFitL[[j]][[i]] - sqrt(validT$degdays)
    # origUnitResid <- sqrtFitL[[j]][[i]]^2 - validT$degdays
    # sqrtcvMSE[j,i] <- mean(sqrtUnitResid^2)
    # sqrtcvErrFrac[j,i] <- sum(sqrtUnitResid^2)/sum( ( sqrt(validT$degdays) - mean(sqrt(validT$degdays)) )^2 )
    # origUnitsqrtcvMSE[j,i] <- mean(origUnitResid^2)
    # origUnitsqrtcvErrFrac[j,i] <- sum(origUnitResid^2)/SSTot
    # rm(sqrtUnitResid, origUnitResid)
  }
}
rm(i, j, validT, SSTot)
# #########################################



combos$avgcvMSE <- apply(cvMSE, 1, mean)
# There may be NAs here if SSTot was 0.  See comments above.
combos$avgcvErrFrac <- apply(cvErrFrac, 1, mean, na.rm=T)

# combos$avgsqrtcvMSE <- apply(sqrtcvMSE, 1, mean)
# combos$avgsqrtcvErrFrac <- apply(sqrtcvErrFrac, 1, mean)
# combos$avgorigUnitsqrtcvMSE <- apply(origUnitsqrtcvMSE, 1, mean)
# combos$avgorigUnitsqrtcvErrFrac <- apply(origUnitsqrtcvErrFrac, 1, mean)

write_csv(combos, file="parallel_leave_out_20perc.csv")


ggplot(data=combos, aes(x=numBtSamps, y=avgcvMSE,
  color=as.factor(numVarSplit))) + geom_line()
# X11()
# ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvMSE,
#   color=as.factor(numVarSplit))) + geom_line()
# # X11()
# ggplot(data=combos, aes(x=numBtSamps, y=avgorigUnitsqrtcvMSE,
#   color=as.factor(numVarSplit))) + geom_line()


ggplot(data=combos, aes(x=numBtSamps, y=avgcvErrFrac,
  color=as.factor(numVarSplit))) + geom_line()
# X11()
# ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvErrFrac,
#   color=as.factor(numVarSplit))) + geom_line()
# X11()
# ggplot(data=combos, aes(x=numBtSamps, y=avgorigUnitsqrtcvErrFrac,
#   color=as.factor(numVarSplit))) + geom_line()
# ####################
