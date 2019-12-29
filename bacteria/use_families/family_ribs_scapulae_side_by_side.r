library("tidyverse")
library("figdim")
library("randomForest")
library("scales")  ## For hue_pal() function.


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
## Find influential taxa for ribs and scapulae, and set colors for
## these to be used consistently across plot panels.


## ########################
## Get top n influential taxa for ribs.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 10

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
ribimportT <- importance(ribRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  top_n(n, wt=`%IncMSE`) %>%
  arrange(`%IncMSE`)
## Remove the "f__" from the family taxon names.
ribimportT$family <- str_remove(ribimportT$family, "f__")

## WORKING HERE!
## We need to remove the following line, leaving it just before we
## make the horizontail bar chart.  Then, we can make the code below
## more readable, since we'll be able to leave out the "as.character"
## parts.
## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
ribimportT$family <- factor(ribimportT$family, levels=ribimportT$family)
## ########################


## ########################
## Get top n influential taxa for scapulae.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 10

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
scapimportT <- importance(scapRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  top_n(n, wt=`%IncMSE`) %>%
  arrange(`%IncMSE`)
## Remove the "f__" from the family taxon names.
scapimportT$family <- str_remove(scapimportT$family, "f__")

## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
scapimportT$family <- factor(scapimportT$family, levels=scapimportT$family)
## ########################
## ##################################################



## ##################################################
## Make four-panel (2 rows by 2 columns) figure for use in publication,
## with rib plots in one column and scapula in the other.


## ########################
## Set up colors so that we can keep them consistent among the 4
## panels. For the taxa displayed in relative abundance plots, use the
## same colors for the bar charts showing the mean % decrease in MSE.

## Will assign colors to "m" most influential taxa for both ribs and
## scapula; some taxa may be influential to both.
m <- 5
mostinfl <- unique( c(
    as.character(ribimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family)),
    as.character(scapimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family))
) )

## The most influential taxa will get colors.  Taxa which are less
## influential, but still appear in the bar chart (top n, but not top
## m), will have bars plotted in gray.
infltaxaColors <- c(hue_pal()(length(mostinfl)))
names(infltaxaColors) <- mostinfl

graytaxaNames <- c(as.character(ribimportT$family), as.character(scapimportT$family))[!c(as.character(ribimportT$family), as.character(scapimportT$family)) %in% mostinfl]
graytaxaColors <- rep("#999999", length(graytaxaNames))
names(graytaxaColors) <- graytaxaNames

## Combine gray and colored taxa into one vector of colors.
taxaColors <- c(infltaxaColors, graytaxaColors)
## ########################


## ########################
## Ribs: Make graph of just %IncMSE alone.

## ## Get the top "n" (whether 8, 10, whatever) influential taxa.
## n <- 10

## ## Turn importance measures into a tibble, sorted by IncNodePurity in
## ## increasing order.
## importanceT <- importance(ribRF) %>%
##   as.data.frame() %>% 
##   rownames_to_column("family") %>%
##   as_tibble() %>%
##   arrange(`%IncMSE`)
## ## Remove the "f__" from the family taxon names.
## importanceT$family <- str_remove(importanceT$family, "f__")

## ## Turn family names into factors, so that we can make the bar chart
## ## with the bars in decreasing order.
## importanceT$family <- factor(importanceT$family, levels=importanceT$family)


## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
ribimportT$family <- factor(ribimportT$family, levels=ribimportT$family)

r1c1Panel <- ggplot(ribimportT, aes(x=family, y=`%IncMSE`, fill=family)) +
  theme_minimal() +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  labs(x=NULL, y="Mean % Decrease in MSE") +
  scale_fill_manual(values=taxaColors)
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
  labs(x=NULL, y="Mean % Decrease in MSE")
## ########################


## ########################
## Ribs: Show line plot of changing relative abundance for the top 5 taxa.

## #####
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
## #####


## The following code was used to label the x-axis with a combo of ADD
## and actual number of days elapsed.
## ## #####
## ## Read in dates and ADD for each sample from a separate CSV file.
## infoT <- read_csv("../sampling_info.csv")

## ## Calculate actual number of days elapsed since the first day of the
## ## study (2016-11-17).
## infoT$days <- as.vector(infoT$date - min(infoT$date))

## ## Build a variable of form "ADD/days" to use with plots.  Make it an
## ## ordered factor to maintain chronological order.
## ADDday <- with(infoT, paste(degdays, days, sep="/"))
## infoT$ADDday <- factor(ADDday, levels=unique(ADDday))

## ## Merge this info with the taxa percentages.
## chooseT <- chooseT %>%
##   left_join(infoT %>% select(sampleName, ADDday, degdays, days))
## ## #####


## #####
## Average the value across samples for each taxa and each day.
## Remove day 0 from figures.
summTopT <- chooseT %>%
  filter(degdays > 0) %>%
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample), medianPercByDay=100*median(fracBySample))
## #####


## #####
## Draw plot of average relative abundance vs. time for these five
## influential taxa.
## dev.new(width=4.5, height=4)
r2c1Panel <- ggplot(summTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon)) +
  scale_y_continuous(limits=c(0, 50), expand=c(0,0)) +
  theme_minimal() +
  theme(##axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5),
        legend.position=c(0.95, 0.98),
        legend.justification=c("right", "top"),
        legend.title=element_blank(),
        legend.key.size=unit(0.5, 'lines'),
        legend.background=element_rect(fill="white")) +
  labs(x="Accumulated Degree Days", y="Relative Abundance (Ribs)")## tag="A")
## ########################


## ########################
## Scapulae: Show line plot of changing relative abundance for the top 5 taxa.

## #####
## For the scapulae family-level taxa, it looks like the first taxa is
## most influential, with the next 6 gradually decreasing in
## importance.  Then, there's another break between the 7th and 8th
## taxa. However, for consistency with the way we did this plot for
## the Forger et al (2019) paper, I'm going to draw the top 5 taxa, as
## measured by %IncMSE.
n <- 5

## Turn importance measures into a tibble, sorted by %IncMSE in
## increasing order.
importanceT <- importance(scapRF) %>%
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
chooseT <- scapT %>%
  filter(taxon %in% topChoices)
## #####


## The following code was used to label the x-axis with a combo of ADD
## and actual number of days elapsed.
## ## #####
## ## Read in dates and ADD for each sample from a separate CSV file.
## infoT <- read_csv("../sampling_info.csv")

## ## Calculate actual number of days elapsed since the first day of the
## ## study (2016-11-17).
## infoT$days <- as.vector(infoT$date - min(infoT$date))

## ## Build a variable of form "ADD/days" to use with plots.  Make it an
## ## ordered factor to maintain chronological order.
## ADDday <- with(infoT, paste(degdays, days, sep="/"))
## infoT$ADDday <- factor(ADDday, levels=unique(ADDday))

## ## Merge this info with the taxa percentages.
## chooseT <- chooseT %>%
##   left_join(infoT %>% select(sampleName, ADDday, degdays, days))
## ## #####


## #####
## Average the value across samples for each taxa and each day.
## Remove day 0 from figures.
summTopT <- chooseT %>%
  filter(degdays > 0) %>%
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample), medianPercByDay=100*median(fracBySample))
## #####


## #####
## Draw plot of average relative abundance vs. time for these five
## influential taxa.
## dev.new(width=4.5, height=4)
r2c2Panel <- ggplot(summTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon)) +
  scale_y_continuous(limits=c(0, 50), expand=c(0,0)) +
  theme_minimal() +
  theme(## axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5),
        legend.position=c(0.95, 0.98),
        legend.justification=c("right", "top"),
        legend.title=element_blank(),
        legend.key.size=unit(0.5, 'lines'),
        legend.background=element_rect(fill="white")) +
  labs(x="Accumulated Degree Days", y="Relative Abundance (Scapulae)")## tag="A")
## ########################

## ########################
library("cowplot")
plot_grid(r1c1Panel, r1c2Panel, r2c1Panel, r2c2Panel, labels=c("a", "b", "c", "d"), nrow=2)

ggsave(file="ribs_scapula_family_4panels.pdf", height=10, width=7.5, units="in")
## ########################
## ##################################################



## ##################################################
## Make two-panel figure showing predicted vs. actual ADD for ribs and
## scapulae.


## ########################
## Ribs: Predicted vs. actual ADD

## This panel is based on Fig. 6c from Forger et al. (2019).  That
## figure was based on Fig. 1c from Pechal et al. (2015).

## Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=ribRF$predicted, actual=ribRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
Rsq <- with(predvactT, round(cor(actual, predicted)^2, 2))
## RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
r3c1Panel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=4925, hjust=0, label=paste("R^2  ==", Rsq), parse=T) +
  annotate("text", x=50, y=4600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, max(as.vector(predvactT))), y=c(0, max(as.vector(predvactT)))) +
  labs(x="Actual Accumulated Degree Days (Ribs)", y="Predicted Accumulated Degree Days (Ribs)")
## ########################


## ########################
## Scapulae: Predicted vs. actual ADD

## This panel is based on Fig. 6c from Forger et al. (2019).  That
## figure was based on Fig. 1c from Pechal et al. (2015).

## Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=scapRF$predicted, actual=scapRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
Rsq <- with(predvactT, round(cor(actual, predicted)^2, 2))
## RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
r3c2Panel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=4925, hjust=0, label=paste("R^2  ==", Rsq), parse=T) +
  annotate("text", x=50, y=4600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, max(as.vector(predvactT))), y=c(0, max(as.vector(predvactT)))) +
  labs(x="Actual Accumulated Degree Days (Scapulae)", y="Predicted Accumulated Degree Days (Scapulae)")
## ########################


## ########################
library("cowplot")
plot_grid(r3c1Panel, r3c2Panel, labels=c("a", "b"), nrow=1)

ggsave(file="ribs_scapula_predicted_vs_actual_ADD_family.pdf", height=4, width=7.5, units="in")
## ########################
## ##################################################

## ##################################################
