library("tidyverse")
library("figdim")
library("randomForest")


## ##################################################
## Read in dataset with family-level taxa.

## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("../", taxalevel, "_massaged.csv"))

## All the family-level taxa have a prefix of "f__" attached to the
## family name.  We want to remove this prefix.  The only taxon that
## is not affected is the "Rare" category, because it's the only one
## which doesn't have the prefix.
allT$taxon <- str_remove(allT$taxon, "f__")

## Filter data to just ribs, and then store.  Then, do the same for scapulae.
ribT <- allT %>% filter(type=="Rib")
scapT <- allT %>% filter(type=="Scapula")

## For both rib and scapula data, remove "Rare" taxa category and put
## into the wide format which was used by the random forest model.
ribwideT <- ribT %>%
  filter(taxon!="Rare") %>%
  select(degdays, sampleName, taxon, fracBySample) %>%
  spread(taxon, fracBySample) %>%
  select(-sampleName)
scapwideT <- scapT %>%
  filter(taxon!="Rare") %>%
  select(degdays, sampleName, taxon, fracBySample) %>%
  spread(taxon, fracBySample) %>%
  select(-sampleName)

rm(allT)
## ##################################################



## ##################################################
## Read in the final fitted models for ribs and scapulae.

load("w_ribs/families_ribs_rfmodel.RData")
## The object was named "rf", but we don't want to mix up the model
## for ribs with the one for scapulae.  So we re-name it.
ribRF <- rf
rm(rf)

load("w_scapulae/families_scapulae_rfmodel.RData")
## As with the model for ribs, we give it a more specific name.
scapRF <- rf
rm(rf)
## ##################################################




## ##################################################
## Make six-panel (3 rows by 2 columns) figure for use in publication,
## with rib plots in one column and scapula in the other.


## ########################
## Ribs: Make graph of just %IncMSE alone.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 10

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
importanceT <- importance(ribRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  arrange(`%IncMSE`)
## Remove the "f__" from the family taxon names.
importanceT$family <- str_remove(importanceT$family, "f__")

## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
importanceT$family <- factor(importanceT$family, levels=importanceT$family)

r1c1Panel <- ggplot(importanceT %>% top_n(n, wt=`%IncMSE`),
                  aes(x=family, y=`%IncMSE`)) +
  theme_minimal() +
  coord_flip() +
  geom_col() +
  labs(x="Family-level Taxa for Ribs", y="Mean % Decrease in MSE When Excluded from Random Forest Model")
## ########################


## ########################
## Scapulae: Make graph of just %IncMSE alone.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 10

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
importanceT <- importance(scapRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  arrange(`%IncMSE`)
## Remove the "f__" from the family taxon names.
importanceT$family <- str_remove(importanceT$family, "f__")

## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
importanceT$family <- factor(importanceT$family, levels=importanceT$family)

r1c2Panel <- ggplot(importanceT %>% top_n(n, wt=`%IncMSE`),
                  aes(x=family, y=`%IncMSE`)) +
  theme_minimal() +
  coord_flip() +
  geom_col() +
  labs(x="Family-level Taxa for Scapulae", y="Mean % Decrease in MSE When Excluded from Random Forest Model")
## ########################


## ########################
## Ribs: Show line plot of changing relative abundance for the top 5 taxa.

## For the rib family-level taxa, it looks like the first 4 taxa are
## the most influential.  However, for consistency with the way we did
## this plot for the Forger et al (2019) paper, I'm going to draw the
## top 5 taxa, as measured by %IncMSE.
n <- 5

## Turn importance measures into a tibble, sorted by %IncMSE in
## increasing order.
importanceT <- importance(ribRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  arrange(`%IncMSE`)
## Remove the "f__" from the family taxon names.
importanceT$family <- str_remove(importanceT$family, "f__")


## Save the names of the families that are in the top 10 in
## terms of %IncMSE.
topChoices <- as.character(importanceT %>% arrange(desc(`%IncMSE`)) %>% pull(family))[1:n]

## Find the percentages for these taxa.
chooseT <- ribT %>%
  filter(taxon %in% topChoices)

## Average the value across samples for each taxa and each day.
summTopT <- chooseT %>%
  group_by(taxon, sampleName) %>%
  summarize(meanPercByDay=100*mean(fracBySample), medianPercByDay=100*median(fracBySample))


## #####
## As in Forger et al (2019), we want a plot of average relative
## abundance vs. time for these five influential taxa.  The x-axis had
## the time steps evenly spaced (not reflecting actual time passage),
## with each tick mark labeled with the day/degreeday.  To do this,
## but yet keep days in order, we need to build a new factor variable
## of the form day/degree day, with ordered levels.

## Read in dates and ADD for each sample.
infoT <- read_csv("../sampling_info.csv")
## Using the dates, calculate how many days have passed since the
## first day of the study (2016-11-17).
days <- as.vector(infoT$date - min(infoT$date))

infoT$dayADD <- factor(with(infoT, paste(degdays, days, sep="/")))
summTopT$dayADD <- factor(with(summTopT, paste(degdays, days, sep="/")), levels=orderedLevels)
rm(orderedLevels)
## #####

dev.new(width=4.5, height=4)
trPanel <- ggplot(summTopT, aes(x=dayADD, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon)) +
  scale_y_continuous(limits=c(0, 100), expand=c(0,0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5),
        legend.position=c(0.95, 0.98),
        legend.justification=c("right", "top"),
        legend.title=element_blank(),
        legend.key.size=unit(0.5, 'lines'),
        legend.background=element_rect(fill="white")) +
  labs(x="Accumulated Degree Days/Days", y="Relative Abundance")## tag="A")
## ########################


## ########################
## This panel is based on Fig. 1c from Pechal et al. (2015).  It
## is predicted vs. actual ADD.

## Make a tibble of actual and predicted values for each observation.
predvactT <- as.tibble(data.frame(predicted=rf$predicted, actual=rf$y))
Rsq <- with(predvactT, round(cor(actual, predicted)^2, 2))
## RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(resids^2)), 2)  
blPanel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=1700, hjust=0, label=paste("R^2  ==", Rsq), parse=T) +
  annotate("text", x=50, y=1600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, max(as.vector(predvactT))), y=c(0, max(as.vector(predvactT)))) +
  labs(x="Actual Accumulated Degree Days", y="Predicted Accumulated Degree Days")
## ########################


## ########################
## This panel shows the same results as the panel above, but on the
## log scale.

## Make new columns with log, base 10.  We have 6 observations made at
## day/ADD 0/0.  Since log10(0) is undefined, we make these values NA.
predvactT$logactual <- with(predvactT, ifelse(actual>0, log10(actual), NA))
predvactT$logpredicted <- with(predvactT, ifelse(predicted>0, log10(predicted), NA))
minAxisLmt <- min(c(predvactT$logactual, predvactT$logpredicted), na.rm=T)
maxAxisLmt <- max(c(predvactT$logactual, predvactT$logpredicted), na.rm=T)

Rsq <- with(predvactT, round(cor(actual, predicted)^2, 2))
## RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(resids^2)), 2)  
blPanel <- ggplot(predvactT, aes(x=logactual, y=logpredicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  ## annotate("text", x=50, y=1700, hjust=0, label=paste("R^2  ==", Rsq), parse=T) +
  ## annotate("text", x=50, y=1600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(minAxisLmt, maxAxisLmt), y=c(minAxisLmt, maxAxisLmt)) +
  labs(x="Log 10 of Actual Accumulated Degree Days", y="Log 10 of Predicted Accumulated Degree Days")
## ########################


## ########################
library("cowplot")
plot_grid(tlPanel, trPanel, blPanel, labels=c("a", "b", "c", "d"), rel_widths=c(1.125, 1), rel_heights=c(1, 1))
ggsave(file="relative_abundance_Rsq_rmse.pdf", width=8.5, height=4, units="in")
## ########################
## ##################################################
