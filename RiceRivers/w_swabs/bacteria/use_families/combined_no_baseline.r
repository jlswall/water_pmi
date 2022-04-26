library("tidyverse")
library("figdim")
library("randomForest")
library("scales")  # For hue_pal() function.
library("cowplot")
library("ggpubr")


# ##################################################
# Read in dataset with family-level taxa.

# Are we dealing with phlya, orders, or families?
taxalevel <- "families"

# Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("../", taxalevel, "_massaged.csv"))

# All the family-level taxa have a prefix of "f__" attached to the
# family name.  We want to remove this prefix.  The only taxon that
# is not affected is the "Rare" category, because it's the only one
# which doesn't have the prefix.
allT$taxon <- str_remove(allT$taxon, "f__")

# Filter data to just ribs, and then store.  Then, do the same for scapulae.
ribT <- allT %>% filter( (type=="Rib") & (degdays > 0) )
scapT <- allT %>% filter( (type=="Scapula") & (degdays > 0) )
bothT <- allT %>% filter(degdays > 0)

rm(allT)
# ##################################################



# ##################################################
# Read in the final fitted models for ribs, scapulae, and for ribs/scapulae
# combined.

load("w_ribs/no_baseline/families_ribs_rfmodel.RData")
# The object was named "rf", but we don't want to mix up the model
# for ribs with the one for scapulae.  So we re-name it.
ribRF <- rf
rm(rf)

load("w_scapulae/no_baseline/families_scapulae_rfmodel.RData")
# As with the model for ribs, we give it a more specific name.
scapRF <- rf
rm(rf)

load("both_ribs_scapulae/no_baseline/families_combined_rfmodel.RData")
# As with the previous models, we give it a more specific name.
bothRF <- rf
rm(rf)
# ##################################################



# ##################################################
# Find influential taxa for ribs and scapulae, and set colors for
# these to be used consistently across plot panels.

# ########################
# Buiid function to extract the top n (8, 10, whatever) most influential taxa
# (measured according to "%IncMSE") from each model, saving the taxa name and
# the measures of influence.
extractInflTaxa <- function(modelRF, n){

  importT <- importance(modelRF) %>%
    as.data.frame() %>% 
    rownames_to_column("family") %>%
    as_tibble() %>%
    top_n(n, wt=`%IncMSE`) %>%
    arrange(`%IncMSE`)

    return(importT)
}
# ########################


# ########################
# Get top 10 influential taxa for ribs, scapulae, and ribs/scapulae combined.
ribimportT <- extractInflTaxa(ribRF, n=10)
scapimportT <- extractInflTaxa(scapRF, n=10)
bothimportT <- extractInflTaxa(bothRF, n=10)

# Remove the "f__" from the family taxon names.
ribimportT$family <- str_remove(ribimportT$family, "f__")
scapimportT$family <- str_remove(scapimportT$family, "f__")
bothimportT$family <- str_remove(bothimportT$family, "f__")
# ########################
# ##################################################




# ##################################################
# Make four-panel (2 rows by 2 columns) figure for use in publication,
# with rib plots in one column and scapula in the other.


# ########################
# Set up colors so that we can keep them consistent among the 4
# panels. For the taxa displayed in relative abundance plots, use the
# same colors for the bar charts showing the mean % decrease in MSE.

# Read in colors assigned to most influential taxa to keep colors the same
# across plots.
colorpaletteDF <- read.csv("../../../color_table_infl_taxa.csv")
assignedColors <- colorpaletteDF$taxaColors
names(assignedColors) <- colorpaletteDF$"taxon"
rm(colorpaletteDF)

# Identify "m" most influential taxa for both ribs and scapula; some taxa may be
# influential to both.
m <- 5
mostinfl <- sort( unique( c(
    ribimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family),
    scapimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family),
    bothimportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family)
) ) )


# The most influential taxa will get the previously determined colors.  Taxa
# which are less influential, and do not appear in the top "m" taxa for either
# ribs or scapulae, will have bars plotted in white with black outlines.
taxaColors <- assignedColors[mostinfl]


# # The most influential taxa will get the previously determined colors.  Taxa
# # which are less influential, and do not appear in the top "m" taxa for either
# # ribs or scapulae, will have bars plotted in gray.
# infltaxaColors <- assignedColors[mostinfl]

# graytaxaNames <- c(ribimportT$family, scapimportT$family)[!c(ribimportT$family,
#   scapimportT$family) %in% mostinfl]
# graytaxaColors <- rep("#999999", length(graytaxaNames))
# names(graytaxaColors) <- graytaxaNames

# # Combine gray and colored taxa into one vector of colors.
# taxaColors <- c(infltaxaColors, graytaxaColors)
# ########################


# ########################
# Set up range for the bars in the bar chart so that axes can be
# consistent between bar charts for ribs and scapulae.

barMax <- ceiling( max( c( ribimportT %>% pull(`%IncMSE`),
                           scapimportT %>% pull(`%IncMSE`),
                           bothimportT %>% pull(`%IncMSE`)
                          )))
# ########################


# ########################
# Ribs: Make graph of just %IncMSE alone.

# Turn family names into factors, so that we can make the bar chart
# with the bars in decreasing order.
ribimportT$family <- factor(ribimportT$family, levels=ribimportT$family)

ribbarPanel <- ggplot(ribimportT, aes(x=family, y=`%IncMSE`, fill=family, color=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  #  labs(x=NULL, y="Mean % Decrease in MSE") +
  labs(x=NULL, y="") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors, na.value = "white") + 
  scale_color_manual(values=taxaColors, na.value = "black")


# ########################


# ########################
# Scapulae: Make graph of just %IncMSE alone.

# Turn family names into factors, so that we can make the bar chart
# with the bars in decreasing order.
scapimportT$family <- factor(scapimportT$family, levels=scapimportT$family)

scapbarPanel <- ggplot(scapimportT, aes(x=family, y=`%IncMSE`, fill=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  labs(x=NULL, y="Mean % Decrease in MSE") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors)
# ########################


# ########################
# Combined ribs/scapulae: Make graph of just %IncMSE alone.

# Turn family names into factors, so that we can make the bar chart
# with the bars in decreasing order.
bothimportT$family <- factor(bothimportT$family, levels=bothimportT$family)

bothbarPanel <- ggplot(bothimportT, aes(x=family, y=`%IncMSE`, fill=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  labs(x=NULL, y="Mean % Decrease in MSE") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors)
# ########################



# ########################
# Build function to find sample averages by degday for each of the influential
# taxa.

calcAvgByTaxaDay <- function(importT, countsT, m){

  # Save the names of the families that are in the top m in terms of %IncMSE.
  topChoices <- as.character(importT %>% top_n(m, wt=`%IncMSE`) %>% pull(family))

  # Find the percentages for these taxa.
  chooseT <- countsT %>%
    filter(taxon %in% topChoices)

  # Average the value across samples for each taxa and each day.
  summTopT <- chooseT %>%
    group_by(taxon, degdays) %>%
    summarize(meanPercByDay=100*mean(fracBySample),
      medianPercByDay=100*median(fracBySample))

  return(summTopT)
}

# Call this function for ribs, scapula, and combined ribs/scapulae.  Set the
# upper y-axis limit for the line plots to accomodate the highest average.
ribSummTopT <- calcAvgByTaxaDay(ribimportT, ribT)
scapSummTopT <- calcAvgByTaxaDay(scapimportT, scapT)
bothSummTopT <- calcAvgByTaxaDay(bothimportT, bothT)

lineMax <- max(ribSummTopT$meanPercByDay,
              scapSummTopT$meanPercByDay,
              bothSummTopT$meanPercByDay)
# ########################



# ########################
# Ribs: Show line plot of changing relative abundance for the top 5 taxa.

# #####
# For the rib family-level taxa, it looks like the first 4 taxa are
# the most influential.  However, for consistency with the way we did
# this plot for the Forger et al (2019) paper, I'm going to draw the
# top m taxa, as measured by %IncMSE.

# Save the names of the families that are in the top m in terms of
# %IncMSE.
topChoices <- as.character(ribimportT %>% top_n(m, wt=`%IncMSE`) %>%
  pull(family))

# Find the percentages for these taxa.
chooseT <- ribT %>%
  filter(taxon %in% topChoices)
# #####


# #####
# Average the value across samples for each taxa and each day.
# Remove day 0 from figures.
summTopT <- chooseT %>%
  filter(degdays > 0) %>%  # For analyses which don't use baseline obs
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample),
    medianPercByDay=100*median(fracBySample))
# #####


# #####
# Draw plot of average relative abundance vs. time for these five
# influential taxa.
# dev.new(width=4.5, height=4)
riblinePanel <- ggplot(summTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_y_continuous(limits=c(0, 25), expand=c(0,0)) +
  theme_minimal() +
  #  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  labs(x="", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  scale_color_manual(values=taxaColors)
# ########################


# ########################
# Scapulae: Show line plot of changing relative abundance for the top 5 taxa.

# #####
# For the scapulae family-level taxa, it looks like the first taxa is
# most influential, with the next 6 gradually decreasing in
# importance.  Then, there's another break between the 7th and 8th
# taxa. However, for consistency with the way we did this plot for
# the Forger et al (2019) paper, I'm going to draw the top m taxa, as
# measured by %IncMSE.

# Save the names of the families that are in the top 10 in
# terms of %IncMSE.
topChoices <- as.character(scapimportT %>% top_n(m, wt=`%IncMSE`) %>%
  pull(family))

# Find the percentages for these taxa.
chooseT <- scapT %>%
  filter(taxon %in% topChoices)
# #####


# #####
# Average the value across samples for each taxa and each day.
# Remove day 0 from figures.
summTopT <- chooseT %>%
  filter(degdays > 0) %>% # For analyses which don't use baseline obs
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample),
    medianPercByDay=100*median(fracBySample))
# #####


# #####
# Draw plot of average relative abundance vs. time for these five
# influential taxa.
# dev.new(width=4.5, height=4)
scaplinePanel <- ggplot(summTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_y_continuous(limits=c(0, 25), expand=c(0,0)) +
  theme_minimal() +
  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  scale_color_manual(values=taxaColors)
# ########################


# ########################
# First row shows info about rib taxa (bar chart, relative abundance graphs).
plotrow1 <- plot_grid(ribbarPanel, riblinePanel, rel_widths=c(1, 2), nrow=1)
plotrow1 <- annotate_figure(plotrow1, left=text_grob("Ribs", face="bold",
  rot=90, size=14))

# First row shows info about rib taxa (bar chart, relative abundance graphs).
plotrow2 <- plot_grid(scapbarPanel, scaplinePanel, rel_widths=c(1, 2), nrow=1)
plotrow2 <- annotate_figure(plotrow2, left=text_grob("Scapulae", face="bold",
  rot=90, size=14))

# Put the rows together to make 4-panel figure.
plot_grid(plotrow1, plotrow2, nrow=2)

ggsave(file="rr_swabs_rib_scapula_family_no_baseline_4panels.pdf", height=5,
  width=7.5, units="in")
# ########################
# ##################################################



# ##################################################
# Make two-panel figure showing predicted vs. actual ADD for ribs and
# scapulae.
# This panel is based on Fig. 6c from Forger et al. (2019).  That
# figure was based on Fig. 1c from Pechal et al. (2015).


# ########################
# Ribs: Predicted vs. actual ADD

# Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=ribRF$predicted, actual=ribRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
# Force Rsq to be printed to 2 decmial places.
Rsq <- with(predvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
ribscatterPanel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=3600, hjust=0, label=paste("R^2  ==", deparse(Rsq)),
    parse=T) +
  annotate("text", x=50, y=3300, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, max(as.vector(predvactT))), y=c(0, max(as.vector(predvactT)))) +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  labs(x="Actual Accumulated Degree Days",
    y="Predicted Accumulated Degree Days")

ribscatterPanel <- annotate_figure(ribscatterPanel, top=text_grob("Ribs",
  face="bold", size=14, vjust=1))
# ########################


# ########################
# Scapulae: Predicted vs. actual ADD

# Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=scapRF$predicted, actual=scapRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
# Force Rsq to be printed to 2 decmial places.
Rsq <- with(predvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
scapscatterPanel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=3600, hjust=0, label=paste("R^2  ==", deparse(Rsq)),
    parse=T) +
  annotate("text", x=50, y=3300, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, max(as.vector(predvactT))), y=c(0, max(as.vector(predvactT)))) +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  #  labs(x="Actual Accumulated Degree Days", y="Predicted Accumulated Degree Days")
  labs(x="Actual Accumulated Degree Days", y="")

scapscatterPanel <- annotate_figure(scapscatterPanel, top=text_grob("Scapulae",
  face="bold", size=14, vjust=1))
# ########################


# ########################
plot_grid(ribscatterPanel, scapscatterPanel, nrow=1)

ggsave(file="rr_swabs_rib_scapula_family_no_baseline_predicted_vs_actual_ADD.pdf",
  height=4, width=7.5, units="in")
# ########################
# ##################################################
