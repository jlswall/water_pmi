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
# ########################


# ########################
# Set up range for the bars in the bar chart so that axes can be
# consistent between bar charts.

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

scapbarPanel <- ggplot(scapimportT, aes(x=family, y=`%IncMSE`, fill=family, color=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  # labs(x=NULL, y="Mean % Decrease in MSE") +
  labs(x=NULL, y="") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors, na.value = "white") + 
  scale_color_manual(values=taxaColors, na.value = "black")
# ########################


# ########################
# Combined ribs/scapulae: Make graph of just %IncMSE alone.

# Turn family names into factors, so that we can make the bar chart
# with the bars in decreasing order.
bothimportT$family <- factor(bothimportT$family, levels=bothimportT$family)

bothbarPanel <- ggplot(bothimportT, aes(x=family, y=`%IncMSE`, fill=family, color=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  labs(x=NULL, y="Mean % Decrease in MSE") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors, na.value = "white") + 
  scale_color_manual(values=taxaColors, na.value = "black")
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
      medianPercByDay=100*median(fracBySample)) %>%
    ungroup()

  return(summTopT)
}

# Call this function for ribs, scapula, and combined ribs/scapulae.  Set the
# upper y-axis limit for the line plots to accomodate the highest average.
ribSummTopT <- calcAvgByTaxaDay(ribimportT, ribT, m)
scapSummTopT <- calcAvgByTaxaDay(scapimportT, scapT, m)
bothSummTopT <- calcAvgByTaxaDay(bothimportT, bothT, m)

lineMax <- max(ribSummTopT$meanPercByDay,
              scapSummTopT$meanPercByDay,
              bothSummTopT$meanPercByDay)
# ########################



# ########################
# Show line plots of changing relative abundance for the top m taxa.

# #####
# Draw plot of average relative abundance vs. time for the top m influential
# taxa for ribs.
riblinePanel <- ggplot(ribSummTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_y_continuous(limits=c(0, lineMax), expand=c(0,0)) +
  theme_minimal() +
  #  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  labs(x="", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  scale_color_manual(values=taxaColors)
# #####


# #####
# Draw plot of average relative abundance vs. time for the top m influential
# taxa for scapulae.
scaplinePanel <- ggplot(scapSummTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_y_continuous(limits=c(0, lineMax), expand=c(0,0)) +
  theme_minimal() +
  # labs(x="Accumulated Degree Days", y="Relative Abundance") +
  labs(x="", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  scale_color_manual(values=taxaColors)
# #####


# #####
# Draw plot of average relative abundance vs. time for the top m influential
# taxa for ribs and scapulae combined.
bothlinePanel <- ggplot(bothSummTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_y_continuous(limits=c(0, lineMax), expand=c(0,0)) +
  theme_minimal() +
  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  scale_color_manual(values=taxaColors)
# #####
# ########################


# ########################
# First row shows info about rib taxa (bar chart, relative abundance graphs).
plotrow1 <- plot_grid(ribbarPanel, riblinePanel, rel_widths=c(1, 2), nrow=1)
plotrow1 <- annotate_figure(plotrow1, left=text_grob("Ribs", face="bold",
  rot=90, size=14))

# Next show info about scapulae taxa (bar chart, relative abundance graphs).
plotrow2 <- plot_grid(scapbarPanel, scaplinePanel, rel_widths=c(1, 2), nrow=1)
plotrow2 <- annotate_figure(plotrow2, left=text_grob("Scapulae", face="bold",
  rot=90, size=14))

# Next show info about scapulae taxa (bar chart, relative abundance graphs).
plotrow3 <- plot_grid(bothbarPanel, bothlinePanel, rel_widths=c(1, 2), nrow=1)
plotrow3 <- annotate_figure(plotrow3, left=text_grob("Ribs & scapulae", face="bold",
  rot=90, size=14))

# Put the rows together to make figure.
plot_grid(plotrow1, plotrow2, plotrow3, nrow=3)

ggsave(file="rr_combined_family_no_baseline_6panels.pdf", height=7.5,
  width=7.5, units="in")
# ########################
# ##################################################



# ##################################################
# Make three-panel figure showing predicted vs. actual ADD for ribs, scapulae,
# and combined ribs/scapulae.


# ########################
# Build tibble for each of the models, containing predicted and actual columns.
# Then, for scatterplots with consistent axes, calculate the range of values.


ribpredvactT <- with(ribRF, as_tibble(data.frame(actual=y, predicted=predicted,
  resids=y-predicted)))
scappredvactT <- with(scapRF, as_tibble(data.frame(actual=y, predicted=predicted,
  resids=y-predicted)))
bothpredvactT <- with(bothRF, as_tibble(data.frame(actual=y, predicted=predicted,
  resids=y-predicted)))

# Find maximum value across all models.  This allows us to set axis range to
# be the same for all plots.
maxAxis <- max(ribpredvactT, scappredvactT, bothpredvactT)
# ########################


# ########################
# Ribs: Predicted vs. actual ADD

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(ribpredvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(ribpredvactT$resids^2)), 2)
ribscatterPanel <- ggplot(ribpredvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  # annotate("text", x=50, y=3600, hjust=0, label=paste("R^2  ==", deparse(Rsq)),
  #  parse=T) +
  lims(x=c(0, maxAxis), y=c(0, maxAxis)) +
  annotate("text", x=0.02*maxAxis, y=0.95*maxAxis, hjust=0,
    label=paste("R^2  ==", deparse(Rsq)), parse=T) +
  annotate("text", x=0.02*maxAxis, y=0.88*maxAxis, hjust=0,
    label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  labs(x="Actual Accumulated Degree Days",
    y="Predicted Accumulated Degree Days")

ribscatterPanel <- annotate_figure(ribscatterPanel, top=text_grob("Ribs",
  face="bold", size=14, vjust=1))
# ########################


# ########################
# Scapulae: Predicted vs. actual ADD

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(scappredvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(scappredvactT$resids^2)), 2)  
scapscatterPanel <- ggplot(scappredvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  lims(x=c(0, maxAxis), y=c(0, maxAxis)) +
  annotate("text", x=0.02*maxAxis, y=0.95*maxAxis, hjust=0,
    label=paste("R^2  ==", deparse(Rsq)), parse=T) +
  annotate("text", x=0.02*maxAxis, y=0.88*maxAxis, hjust=0,
    label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  #  labs(x="Actual Accumulated Degree Days", y="Predicted Accumulated Degree Days")
  labs(x="Actual Accumulated Degree Days", y="")

scapscatterPanel <- annotate_figure(scapscatterPanel, top=text_grob("Scapulae",
  face="bold", size=14, vjust=1))
# ########################


# ########################
# Both ribs and scapulae: Predicted vs. actual ADD

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(bothpredvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(bothpredvactT$resids^2)), 2)  
bothscatterPanel <- ggplot(bothpredvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  lims(x=c(0, maxAxis), y=c(0, maxAxis)) +
  annotate("text", x=0.02*maxAxis, y=0.95*maxAxis, hjust=0,
    label=paste("R^2  ==", deparse(Rsq)), parse=T) +
  annotate("text", x=0.02*maxAxis, y=0.88*maxAxis, hjust=0,
    label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  labs(x="Actual Accumulated Degree Days", y="Predicted Accumulated Degree Days")

bothscatterPanel <- annotate_figure(bothscatterPanel,
  top=text_grob("Combined ribs and scapulae", face="bold", size=14, vjust=1))
# ########################


# ########################
topRowPlots <- plot_grid(ribscatterPanel, scapscatterPanel, nrow=1)
plot_grid(topRowPlots, bothscatterPanel, nrow=2)

ggsave(file="rr_combined_family_no_baseline_predicted_vs_actual_ADD.pdf",
  height=7, width=7, units="in")
# ########################
# ##################################################
