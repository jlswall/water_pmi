library("tidyverse")
library("figdim")
library("randomForest")
library("scales")  ## For hue_pal() function.
library("cowplot")
library("ggpubr")


## ##################################################
## Read in dataset with family-level taxa.

## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0(taxalevel, "_massaged.csv"))

## ## All the family-level taxa have a prefix of "f__" attached to the
## ## family name.  We want to remove this prefix.  The only taxon that
## ## is not affected is the "Rare" category, because it's the only one
## ## which doesn't have the prefix.
## allT$taxon <- str_remove(allT$taxon, "f__")

## Filter data to just ribs, and then store.  Then, do the same for scapulae.
ribT <- allT %>% filter(type=="Rib")
scapT <- allT %>% filter(type=="Scapula")
## Change "_unclassified" to "_uncl" in the family taxon names to make
## the names shorter to display on plots.
ribT$taxon <- str_replace(ribT$taxon, "_unclassified", "_uncl")
scapT$taxon <- str_replace(scapT$taxon, "_unclassified", "_uncl")

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
## Change "_unclassified" to "_uncl" in the family taxon names to make
## the names shorter to display on plots.
ribimportT$family <- str_replace(ribimportT$family, "_unclassified", "_uncl")
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
## Change "_unclassified" to "_uncl" in the family taxon names to make
## the names shorter to display on plots.
scapimportT$family <- str_replace(scapimportT$family, "_unclassified", "_uncl")
## ########################
## ##################################################



## ##################################################
## Find mean relative abundance across time for top m influential taxa
## for both ribs and scapulae.  This will allow us to set consistent
## y-axes for rib plot and scapula plot.

## Set number of top taxa to plot.
m <- 5

## #####
## For ribs:

## Identify top m taxa, pull their percentages.
topChoicesRib <- as.character(ribimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family))
## Find the percentages for these taxa for each ADD.
chooseT <- ribT %>%
  filter(taxon %in% topChoicesRib)

## Average the value across samples for each taxa and each day.
summTopRibT <- chooseT %>%
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample), medianPercByDay=100*median(fracBySample))
rm(chooseT)
## #####


## #####
## For scapulae:

## Identify top m taxa, pull their percentages.
topChoicesScap <- as.character(scapimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family))
## Find the percentages for these taxa for each ADD.
chooseT <- scapT %>%
  filter(taxon %in% topChoicesScap)

## Average the value across samples for each taxa and each day.
summTopScapT <- chooseT %>%
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample), medianPercByDay=100*median(fracBySample))
rm(chooseT)
## #####


## Find the highest mean relative abundance for both ribs and
## scapulae.  This will allow us specifiy consistent y-axes for the
## relative abundance plots.
relAbundMax <- ceiling( max( c(summTopRibT %>% pull(meanPercByDay), summTopScapT %>% pull(meanPercByDay)) ) )
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
mostinfl <- sort(unique(c(summTopRibT$taxon, summTopScapT$taxon)))
## m <- 5
## mostinfl <- sort( unique( c(
##     ribimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family),
##     scapimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family)
## ) ) )

## The most influential taxa will get colors.  Taxa which are less
## influential, and don not appear in the top "m" taxa for either ribs
## or scapulae, will have bars plotted in gray.
infltaxaColors <- c(hue_pal()(length(mostinfl)))
names(infltaxaColors) <- mostinfl

graytaxaNames <- c(ribimportT$family, scapimportT$family)[!c(ribimportT$family, scapimportT$family) %in% mostinfl]
graytaxaColors <- rep("#999999", length(graytaxaNames))
names(graytaxaColors) <- graytaxaNames

## Combine gray and colored taxa into one vector of colors.
taxaColors <- c(infltaxaColors, graytaxaColors)
## ########################


## ########################
## Set up range for the bars in the bar chart so that axes can be
## consistent between bar charts for ribs and scapulae.

barMax <- ceiling(max(c(ribimportT %>% pull(`%IncMSE`), scapimportT %>% pull(`%IncMSE`))))
## ########################


## ########################
## Ribs: Make graph of just %IncMSE alone.

## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
ribimportT$family <- factor(ribimportT$family, levels=ribimportT$family)

ribbarPanel <- ggplot(ribimportT, aes(x=family, y=`%IncMSE`, fill=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  ##  labs(x=NULL, y="Mean % Decrease in MSE") +
  labs(x=NULL, y="") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors)
## ########################


## ########################
## Scapulae: Make graph of just %IncMSE alone.

## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
scapimportT$family <- factor(scapimportT$family, levels=scapimportT$family)

scapbarPanel <- ggplot(scapimportT, aes(x=family, y=`%IncMSE`, fill=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  labs(x=NULL, y="Mean % Decrease in MSE") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors)
## ########################


## ########################
## Ribs: Show line plot of changing relative abundance for the top m taxa.

## Draw plot of average relative abundance vs. time (ADD) for these
## top m influential taxa.
riblinePanel <- ggplot(summTopRibT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_y_continuous(limits=c(0, relAbundMax), expand=c(0,0)) +
  theme_minimal() +
  ##  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  labs(x="", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
  scale_color_manual(values=taxaColors)
## ########################


## ########################
## Scapulae: Show line plot of changing relative abundance for the top m taxa.

## Draw plot of average relative abundance vs. time (ADD) for these
## top m influential taxa.
scaplinePanel <- ggplot(summTopScapT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_y_continuous(limits=c(0, relAbundMax), expand=c(0,0)) +
  theme_minimal() +
  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
  scale_color_manual(values=taxaColors)
## ########################


## ########################
## First row shows info about rib taxa (bar chart, relative abundance graphs).
plotrow1 <- plot_grid(ribbarPanel, riblinePanel, rel_widths=c(1, 2), nrow=1)
plotrow1 <- annotate_figure(plotrow1, left=text_grob("Ribs", face="bold", rot=90, size=14))

## First row shows info about rib taxa (bar chart, relative abundance graphs).
plotrow2 <- plot_grid(scapbarPanel, scaplinePanel, rel_widths=c(1, 2), nrow=1)
plotrow2 <- annotate_figure(plotrow2, left=text_grob("Scapulae", face="bold", rot=90, size=14))

## Put the rows together to make 4-panel figure.
plot_grid(plotrow1, plotrow2, nrow=2)

ggsave(file="hl_rib_scapula_family_4panels.pdf", height=5, width=7.5, units="in")
## ########################
## ##################################################



## ##################################################
## Make two-panel figure showing predicted vs. actual ADD for ribs and
## scapulae.
## This panel is based on Fig. 6c from Forger et al. (2019).  That
## figure was based on Fig. 1c from Pechal et al. (2015).


## ########################
## Ribs: Predicted vs. actual ADD

## Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=ribRF$predicted, actual=ribRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
## Force Rsq to be printed to 2 decmial places.
Rsq <- with(predvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
## RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
ribscatterPanel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=4925, hjust=0, label=paste("R^2  ==", deparse(Rsq)), parse=T) +
  annotate("text", x=50, y=4600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, max(as.vector(predvactT))), y=c(0, max(as.vector(predvactT)))) +
  theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
  labs(x="Actual Accumulated Degree Days", y="Predicted Accumulated Degree Days")

ribscatterPanel <- annotate_figure(ribscatterPanel, top=text_grob("Ribs", face="bold", size=14, vjust=1))
## ########################


## ########################
## Scapulae: Predicted vs. actual ADD

## Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=scapRF$predicted, actual=scapRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
## Force Rsq to be printed to 2 decmial places.
Rsq <- with(predvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
## RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
scapscatterPanel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=4925, hjust=0, label=paste("R^2  ==", deparse(Rsq)), parse=T) +
  annotate("text", x=50, y=4600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, max(as.vector(predvactT))), y=c(0, max(as.vector(predvactT)))) +
  theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
  ##  labs(x="Actual Accumulated Degree Days", y="Predicted Accumulated Degree Days")
  labs(x="Actual Accumulated Degree Days", y="")

scapscatterPanel <- annotate_figure(scapscatterPanel, top=text_grob("Scapulae", face="bold", size=14, vjust=1))
## ########################


## ########################
plot_grid(ribscatterPanel, scapscatterPanel, nrow=1)

ggsave(file="hl_rib_scapula_family_predicted_vs_actual_ADD.pdf", height=4, width=7.5, units="in")
## ########################
## ##################################################
