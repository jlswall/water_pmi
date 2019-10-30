library("tidyverse")
library("randomForest")
library("figdim")
library("parallel")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("../../", taxalevel, "_massaged.csv"))
## ##################################################


## ##################################################
## Are we dealing with ribs, scapulae, or water observations?
obstype <- "Rib"

## Filter the data to just that type.
taxaT <- allT %>% filter(type==obstype)
rm(allT)
## ##################################################



## ##################################################
## Put the data in wide format; remove rare taxa.

## Move back to wide format.  Leave sampleName in so that we can
## distinguish data belonging to individual ribs.
wideT <- taxaT %>%
  filter(taxon!="Rare") %>%
  select(degdays, sampleName, taxon, fracBySample) %>%
  spread(taxon, fracBySample)

## Make variable denoting individual ribs.
wideT$ribnum <- substring(wideT$sampleName, first=3, last=4)

rm(taxaT)
## ##################################################



## ##################################################
## From earlier experiments, we figured out these parameters work best
## for the random forest model.

## Number of bootstrap samples.
numBtSamps <- 3000

## Repeated cross-validation runs (1000 of them), leaving out 20% of
## the observations at a time, indicated that the number of variables
## to consider at each split is about 14 (also 15 is very close)
## for the response variable in the original units.
numVarSplit <- 14
## ##################################################



## ##################################################
## Set up these training and valiation datasets, plus matrices to hold
## cross-valiation results.

## We exclude in turn each combination of individual rib and degree
## day. If there were no missing observations, this would be 20 degree
## days x 5 individual ribs = 100 combos.  However, we have no
## observations on 18 different combos of rib and degdays.
excludeCombos <- wideT %>% select(degdays, ribnum) %>% distinct()
                 
## How many times do we want to run each exclusion combo?
numRunsEachCombo <- 1
excludeT <- NULL
for (i in 1:numRunsEachCombo)
    excludeT <- bind_rows(excludeT, excludeCombos)
rm(excludeCombos)
## Ensure ordering is by degdays and then ribnum.
excludeT <- excludeT %>% arrange(degdays, ribnum)

## Set up the training and validation datasets corresponding to each
## combo.
numCVs <- nrow(excludeT)
crossvalidL <- vector("list", numCVs)
for (i in 1:numCVs){
  lvOut <- (wideT$ribnum==pull(excludeT[i,"ribnum"])) | (wideT$degdays==pull(excludeT[i,"degdays"]))
  trainT <- wideT[!lvOut,] %>% select(-ribnum, -sampleName)
  validT <- wideT[lvOut,]
  crossvalidL[[i]] <- list(trainT=trainT, validT=validT)
}
rm(i, lvOut, trainT, validT)
## #########################################



## #########################################
## Set up function for fitting random forest model using original
## units.
origUnitsF <- function(x, mtry, ntree){
  rf <- randomForest(degdays ~ ., data=x$trainT, mtry=mtry, ntree=ntree, importance=T)
  return(predict(rf, newdata=x$validT))
}

## Set random seed for reproducibility.
set.seed(466267)

## Try using lapply to fit the random forests.
origFitL <- mclapply(crossvalidL, mc.cores=6, origUnitsF, mtry=numVarSplit, ntree=numBtSamps)
## #########################################



## #########################################
## Collect the residuals, making a note about which day and
## individual were left out.

## Set up vectors to hold cross-validation results.
cvMSE <- rep(NA, numCVs)
cvErrFrac <- rep(NA, numCVs)

residDF <- NULL
for (i in 1:numCVs){

  ## Get the validation set for this run from the list.
  validT <- crossvalidL[[i]][["validT"]]

  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (validT$degdays-mean(validT$degdays))^2 )

  ## Calculate the residuals for this validation set.
  resid <- validT$degdays - origFitL[[i]]

  ## Overall cross-validations statistics, using all residuals.
  cvMSE[i] <- mean(resid^2)
  cvErrFrac[i] <- sum(resid^2)/SSTot
  
  ## ## Find residuals which correspond to the degree day which was left
  ## ## out in the training dataset.
  ## dayOfInterest <- excludeMat[i,"degdays"]
  ## whichItems <- validT$degdays==dayOfInterest
  ## residOfInterest <- resid[whichItems]

  ## Build a data frame with these residuals, along with the day and
  ## individual that were left out in this validation.
  iresidDF <- data.frame(dayOmit=pull(excludeT[i,"degdays"]),
                         ribnumOmit=pull(excludeT[i,"ribnum"]),
                         ribnumactual=validT$ribnum,
                         yactual=validT$degdays,
                         yhat=origFitL[[i]],
                         resid=resid)
  ## Add this data.frame to the what we've already collected.
  residDF <- rbind(residDF, iresidDF)
}
rm(i, validT, resid, iresidDF)

## Write this info out.
write.csv(residDF, file="resids_leave_out_one_rib_and_one_day.csv", row.names=FALSE)
## #########################################



## #########################################
## We're interested in the error we're likely to get with the model in
## "regular use".  That means that we are interested the prediction
## the model makes for a rib and time slot that we're never
## observed.  That's only ONE prediction of interest per model run.
## So, we look at the root of the mean squared errors for the n
## predictions of interest in these n runs.

myresids <- residDF %>%
  filter((ribnumactual==ribnumOmit) & (dayOmit==yactual)) %>%
  pull(resid)
sqrt(mean(myresids^2))
## This is about 741-744, whether I use 1, 10, or 100 runs per
## combo.
## #########################################



## #########################################
## Make plot showing the residuals associated with days which were
## completely left out of the model.

ggplot(residDF %>%
       filter(dayOmit==yactual) %>%
       filter(subjOmit==subjactual),
       aes(x=yactual, y=resid)) +
  geom_point() +
  ## geom_point(aes(col=subjOmit)) +
  geom_hline(yintercept=0) +
  labs(x="Actual degree days", y="Error (actual - estimated)")
ggsave(filename="leave_out_one_subj_and_one_day_residuals.pdf", height=3.5, width=3.5, units="in")

ggplot(residDF, aes(x=yactual, y=resid)) +
  facet_wrap(~subjOmit) +
  geom_point(aes(col=subjOmit)) +
  labs(x="Actual degree day", y="Residual")
## #########################################

