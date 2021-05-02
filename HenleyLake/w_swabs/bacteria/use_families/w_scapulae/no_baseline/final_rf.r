library("tidyverse")
library("randomForest")
library("figdim")
library("parallel")


# ##################################################
# Are we dealing with phlya, orders, or families?
taxalevel <- "families"

# Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("../../../", taxalevel, "_massaged.csv"))
# ##################################################



# ##################################################
# Are we dealing with ribs, scapulae, or water observations?
obstype <- "Scapula"

# Filter the data to just that type.
taxaT <- allT %>% filter((type==obstype) & (degdays>0))
rm(allT)
# ##################################################



# ##################################################
# Put the data in wide format; remove days, subj, and rare taxa.

# Move back to wide format.
wideT <- taxaT %>%
  filter(taxon!="Rare") %>%
  select(degdays, sampleName, taxon, fracBySample) %>%
  spread(taxon, fracBySample) %>%
  select(-sampleName)

# rm(taxaT)  ## Use to make plot of influential taxa at finish.
# ##################################################



# ##################################################
# From earlier experiments, we figured out these parameters work best
# for the random forest model.

# Number of bootstrap samples.
numBtSamps <- 3000

# Repeated cross-validation runs (1000 of them), leaving out 20% of the
# observations at a time, indicated that the number of variables to consider at
# each split is about 19 (but 22 and 21 are close) for the response variable in
# the original units.
numVarSplit <- 19
## ##################################################



## ##################################################
## Fit the random forest using all the data (do not hold any back for
## cross-validation).

## ###########################
## Set up function for fitting random forest model using full dataset.
fitFullDataF <- function(x, mtry, ntree){

  rf <- randomForest(degdays ~ . , data=x, mtry=mtry, ntree=ntree, importance=T)

  ## Order the taxa according to decreasing values of the scaled
  ## importance %IncMSE.  Return as a tibble, with %IncMSE column
  ## renamed to PercIncMSE.
  importanceT <- as.tibble(importance(rf), rownames="taxa") %>%
        rename(PercIncMSE=`%IncMSE`) %>%
        arrange(desc(PercIncMSE))

  myresults <- list(predicted=rf[["predicted"]], y=rf[["y"]],
        importanceT=importanceT)

  return(myresults)
}
## ###########################


## ###########################
## Set up list as input for the function, which is just the original
## dataset repeated over and over.
numRepeat <- 1000
repDataL <- vector("list", numRepeat)
for (i in 1:numRepeat)
  repDataL[[i]] <- wideT
## ###########################


## ###########################
## Now, fit random forests to the full dataset over and over.

set.seed(387045)
fullResultsL <- mclapply(repDataL, mc.cores=6, fitFullDataF, mtry=numVarSplit,
  ntree=numBtSamps)

## Calculate the RMSE and pseudo-Rsquared for these runs with the full
## dataset.
fullRMSE <- rep(NA, 1000)
fullRsq <- rep(NA, 1000)
fullImportanceT <- NULL
for (i in 1:length(fullResultsL)){

  ## Get results from run i.
  iTmp <- fullResultsL[[i]]

  ## Find residuals:
  resids <- iTmp$predicted - iTmp$y

  ## Calculate the RMSE and pseudo-Rsquared:
  fullRMSE[i] <- sqrt( mean( resids^2 ) )
  fullRsq[i] <- 1.0 - ( sum(resids^2)/sum( (iTmp$y - mean(iTmp$y))^2 ) )

  ## Store measures of importance in a long tibble.
  fullImportanceT <- rbind(fullImportanceT, iTmp$importanceT)
}
rm(i, iTmp, resids)

## See summary of %IncMSE (measure of importance) over all model runs.
fullImportanceT %>%
  group_by(taxa) %>%
  summarize(meanPercIncMSE=mean(PercIncMSE),
    lbPercIncMSE=quantile(PercIncMSE, 0.025),
    ubPercIncMSE=quantile(PercIncMSE, 0.975)
    ) %>%
  arrange(desc(meanPercIncMSE))

## Get summary statistics for report.
c(mean(fullRMSE), 1.96*sd(fullRMSE))
## RMSE: 521.95793  12.97127
c(mean(fullRsq), 1.96*sd(fullRsq))
## Rsq: 0.827283672 0.008583558

write_csv(data.frame(fullRMSE, fullRsq),
    file="cvstats_w_full_dataset_final_params.csv")
rm(fullRMSE, fullRsq)
# ###########################
# ##################################################



# ##################################################
# Fit the final random forest with all the data (no cross-validation).

set.seed(2920675)

# Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(degdays ~ . , data=wideT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

## init.fig.dimen(file=paste0("orig_units_all_data_families_imp_plot.pdf"), width=8, height=6)
## varImpPlot(rf, main="Importance of eukaryotic family-level taxa (all time steps)")
## dev.off()


# Find residuals:
resids <- rf$predicted - wideT$degdays

# Print out RMSE:
sqrt( mean( resids^2 ) )
# RMSE: 526.454

# Estimate of explained variance, which R documentation calls "pseudo
# R-squared"
1 - ( sum(resids^2)/sum( (wideT$degdays - mean(wideT$degdays))^2 ) )
# Expl. frac.: 0.8243236

# Save the fitted model so that we can re-create graphics and summary
# statistics without running it again.
save(rf, file="families_scapulae_rfmodel.RData")
# ##################################################



## ## ##################################################
## ## Make graph of just IncNodePurity alone.

## ## Get the top "n" (whether 8, 10, whatever) influential taxa.
## n <- 8

## ## Turn importance measures into a tibble, sorted by IncNodePurity in
## ## increasing order.
## importanceT <- importance(rf) %>%
##   as.data.frame() %>% as_tibble() %>%
##   rownames_to_column("family") %>%
##   arrange(IncNodePurity)
## ## Turn family names into factors, so that we can make the bar chart
## ## with the bars in decreasing order.
## importanceT$family <- factor(importanceT$family, levels=importanceT$family)
## ggplot(importanceT %>% top_n(n, wt=IncNodePurity),
##        aes(x=family, y=IncNodePurity)) +
##   coord_flip() +
##   geom_col() +
##   labs(x="Eukaryotic family-level taxa", y="Decrease in node impurity")
## ggsave(filename="orig_units_all_data_families_IncNodePurity_barchart.pdf", height=2.5, width=4.5, units="in")
## ## ##################################################



# ##################################################
# Make graph of just %IncMSE alone.

# Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 15

# Turn importance measures into a tibble, sorted by IncNodePurity in
# increasing order.
importanceT <- importance(rf) %>%
  as.data.frame() %>%
  rownames_to_column("family") %>%
   as_tibble() %>%
  arrange(`%IncMSE`)
# Turn family names into factors, so that we can make the bar chart
# with the bars in decreasing order.
importanceT$family <- factor(importanceT$family, levels=importanceT$family)
ggplot(importanceT %>% top_n(n, wt=`%IncMSE`),
       aes(x=family, y=`%IncMSE`)) +
  coord_flip() +
  geom_col() +
  labs(x="Scapulae: family-level taxa", y="Mean % increase in MSE when excluded")
ggsave(filename="families_scapula_PercIncMSE_barchart.pdf", height=4.5, width=6,
  units="in")
# ##################################################



# ##################################################
# Make scatter plots of the percentages over time (by taxa) for the
# top n taxa in terms of %IncMSE.

# Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 6

# Save the names of the families that are in the top 10 in
# terms of %IncMSE.
topChoices <- as.character(importanceT %>%
  arrange(desc(`%IncMSE`)) %>%
  pull(family))[1:n]

# Find the percentages for these taxa.
chooseT <- taxaT %>%
  filter(taxon %in% topChoices)
chooseT$taxon <- factor(chooseT$taxon, levels=topChoices)


# From sampleName variable, extract the scapulanumber.  Include the scapula
# number is the scatterplot, so that we can check whether one scapula has
# frequently unusual data associated with it.
chooseT$scapnum <- substring(chooseT$sampleName, first=3, last=4)


ggplot(chooseT, aes(degdays, fracBySample)) +
  geom_point(aes(color=scapnum)) +
  labs(x="Degree days", y="Fraction", color="Scapula") +
  theme(legend.title=element_text(size=rel(0.8)),
    legend.text=element_text(size=rel(0.8))) + 
  # Allow diff. y-scales across panels.
    facet_wrap(~taxon, ncol=3, scales="free_y") 
  # facet_wrap(~taxon, ncol=3)  ## Keep y-scales same across panels.
ggsave("infl_scapula_family_scatter.pdf", width=8, height=4, units="in")
# ##################################################
