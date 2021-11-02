library("tidyverse")
library("figdim")
library("randomForest")
library("scales")  # For hue_pal() function.
library("cowplot")
library("ggpubr")


# ##################################################
# Read in dataset with bones' family-level taxa.

# Are we dealing with phlya, orders, or families?
taxalevel <- "families"

# Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("w_bones/bacteria/", taxalevel, "_massaged.csv"))

# All the family-level taxa have a prefix of "f__" attached to the
# family name.  We want to remove this prefix.  The only taxon that
# is not affected is the "Rare" category, because it's the only one
# which doesn't have the prefix.
allT$taxon <- str_remove(allT$taxon, "f__")

# Filter data to just scapulae and no baseline observations.
bonesT <- allT %>% filter( (type=="Scapula") & (degdays > 0))

# Remove "Rare" taxa category and put into the wide format which was used by the
# random forest model.
boneswideT <- bonesT %>%
  filter(taxon!="Rare") %>%
  select(degdays, sampleName, taxon, fracBySample) %>%
  spread(taxon, fracBySample) %>%
  select(-sampleName)

rm(allT)
# ##################################################


# ##################################################
# Read in dataset with swabs' family-level taxa.

# Are we dealing with phlya, orders, or families?
taxalevel <- "families"

# Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("w_swabs/bacteria/", taxalevel, "_massaged.csv"))

# All the family-level taxa have a prefix of "f__" attached to the
# family name.  We want to remove this prefix.  The only taxon that
# is not affected is the "Rare" category, because it's the only one
# which doesn't have the prefix.
allT$taxon <- str_remove(allT$taxon, "f__")

# Filter data to just ribs and no baseline observations.
swabsT <- allT %>% filter( (type=="Scapula") & (degdays > 0))

# For both datasets, remove "Rare" taxa category and put into the wide format
# which was used by the random forest model.
swabswideT <- swabsT %>%
  filter(taxon!="Rare") %>%
  select(degdays, sampleName, taxon, fracBySample) %>%
  spread(taxon, fracBySample) %>%
  select(-sampleName)

rm(allT)
# ##################################################



# ##################################################
# Read in the final fitted models for both.

load("w_bones/bacteria/use_families/w_scapulae/no_baseline/families_scapulae_rfmodel.RData")
# The object was named "rf", but we don't want to mix up the 2 models being
# compared.  So we re-name it.
wBonesRF <- rf
rm(rf)

load("w_swabs/bacteria/use_families/w_scapulae/no_baseline/families_scapulae_rfmodel.RData")
# As before, we give the saved object a more specific name.
wSwabsRF <- rf
rm(rf)
# ##################################################



# ##################################################
# Find influential taxa for both analyses, and set colors for these to be used
# consistently across plot panels.


# ########################
# Get top n influential taxa for first analysis.

# Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 10

# Turn importance measures into a tibble, sorted by IncNodePurity in
# increasing order.
wBonesImportT <- importance(wBonesRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  top_n(n, wt=`%IncMSE`) %>%
  arrange(`%IncMSE`)
# Remove the "f__" from the family taxon names.
wBonesImportT$family <- str_remove(wBonesImportT$family, "f__")
# ########################


# ########################
# Get top n influential taxa for second analysis.

# Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 10

# Turn importance measures into a tibble, sorted by IncNodePurity in
# increasing order.
wSwabsImportT <- importance(wSwabsRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  top_n(n, wt=`%IncMSE`) %>%
  arrange(`%IncMSE`)
# Remove the "f__" from the family taxon names.
wSwabsImportT$family <- str_remove(wSwabsImportT$family, "f__")
# ########################
# ##################################################



# ##################################################
# Make four-panel (2 rows by 2 columns) figure, with analysis 1 results in one
# column and analysis 2 results in the other.


# ########################
# Set up colors so that we can keep them consistent among the 4
# panels. For the taxa displayed in relative abundance plots, use the
# same colors for the bar charts showing the mean % decrease in MSE.

# Read in colors assigned to most influential taxa to keep colors the same
# across plots.
colorpaletteDF <- read.csv("color_table_infl_taxa.csv")
assignedColors <- colorpaletteDF$taxaColors
names(assignedColors) <- colorpaletteDF$"taxon"
rm(colorpaletteDF)

# Identify "m" most influential taxa for both analyses; some taxa may be
# influential to both.
m <- 5
mostinfl <- sort( unique( c(
    wBonesImportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family),
    wSwabsImportT %>% top_n(m, wt=`%IncMSE`) %>% pull(family)
) ) )

# The most influential taxa will get the previously determined colors.  Taxa
# which are less influential, and do not appear in the top "m" taxa for either
# analysis 1 or 2, will have bars plotted in gray.
infltaxaColors <- assignedColors[mostinfl]

graytaxaNames <- c(wBonesImportT$family, wSwabsImportT$family)[!c(wBonesImportT$family, wSwabsImportT$family) %in% mostinfl]
graytaxaColors <- rep("#999999", length(graytaxaNames))
names(graytaxaColors) <- graytaxaNames

# Combine gray and colored taxa into one vector of colors.
taxaColors <- c(infltaxaColors, graytaxaColors)
# ########################


# ########################
# Set up range for the bars in the bar chart so that axes can be
# consistent between bar charts for both analyses.

barMax <- ceiling(max(c(wBonesImportT %>% pull(`%IncMSE`),
 wSwabsImportT %>% pull(`%IncMSE`))))
# ########################


# ########################
# Analysis 1 (with bones): Make graph of just %IncMSE alone.

# Turn family names into factors, so that we can make the bar chart
# with the bars in decreasing order.
wBonesImportT$family <- factor(wBonesImportT$family,
  levels=wBonesImportT$family)

wBonesbarPanel <- ggplot(wBonesImportT, aes(x=family, y=`%IncMSE`,
    fill=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  #  labs(x=NULL, y="Mean % Decrease in MSE") +
  labs(x=NULL, y="") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors)
# ########################


# ########################
# Analysis 2 (with swabs): Make graph of just %IncMSE alone.

# Turn family names into factors, so that we can make the bar chart
# with the bars in decreasing order.
wSwabsImportT$family <- factor(wSwabsImportT$family,
  levels=wSwabsImportT$family)

wSwabsbarPanel <- ggplot(wSwabsImportT, aes(x=family, y=`%IncMSE`,
    fill=family)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0, barMax), expand=c(0,0)) +
  coord_flip() +
  geom_col(show.legend=FALSE) +
  labs(x=NULL, y="Mean % Decrease in MSE") +
  theme(axis.title.x = element_text(size=10)) +
  scale_fill_manual(values=taxaColors)
# ########################


# ########################
# Analysis 1: Show line plot of changing relative abundance for the top 5 taxa.

# #####
# For consistency with the way we did this plot for previous analyses, I'm going
# to draw the top m taxa, as measured by %IncMSE.

# Save the names of the families that are in the top m in terms of
# %IncMSE.
topChoices <- as.character(wBonesImportT %>% top_n(m, wt=`%IncMSE`) %>%
  pull(family))

# Find the percentages for these taxa.
chooseT <- bonesT %>%
  filter(taxon %in% topChoices)
# #####


# #####
# Average the value across samples for each taxa and each day.
summTopT <- chooseT %>%
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample),
    medianPercByDay=100*median(fracBySample))
# #####


# #####
# Draw plot of average relative abundance vs. time for these five
# influential taxa.
# dev.new(width=4.5, height=4)
wBonesLinePanel <- ggplot(summTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
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
# Analysis 2: Show line plot of changing relative abundance for the top 5 taxa.

# #####
# For consistency with the way we did this plot for previous analyses, I'm going
# to draw the top m taxa, as measured by %IncMSE.

# Save the names of the families that are in the top 10 in
# terms of %IncMSE.
topChoices <- as.character(wSwabsImportT %>% top_n(m, wt=`%IncMSE`) %>%
  pull(family))

# Find the percentages for these taxa.
chooseT <- swabsT %>%
  filter(taxon %in% topChoices)
# #####


# #####
# Average the value across samples for each taxa and each day.
summTopT <- chooseT %>%
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample),
    medianPercByDay=100*median(fracBySample))
# #####


# #####
# Draw plot of average relative abundance vs. time for these five
# influential taxa.
# dev.new(width=4.5, height=4)
wSwabsLinePanel <- ggplot(summTopT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_y_continuous(limits=c(0, 25), expand=c(0,0)) +
  theme_minimal() +
  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  scale_color_manual(values=taxaColors)
# ########################


# ########################
# First row shows info about taxa in first analysis (bar chart, relative
# abundance graphs).
plotrow1 <- plot_grid(wBonesbarPanel, wBonesLinePanel, rel_widths=c(1, 2),
  nrow=1)
plotrow1 <- annotate_figure(plotrow1,
  left=text_grob("Using bones, without baseline obs", face="bold", rot=90,
      size=10))

# Second row shows info about taxa in second analysis (bar chart, relative
# abundance graphs).
plotrow2 <- plot_grid(wSwabsbarPanel, wSwabsLinePanel, rel_widths=c(1, 2),
  nrow=1)
plotrow2 <- annotate_figure(plotrow2,
  left=text_grob("Using swabs, without baseline obs", face="bold", rot=90,
      size=10))

# Put the rows together to make 4-panel figure.
plot_grid(plotrow1, plotrow2, nrow=2)

ggsave(file="rr_scapula_bones_vs_swabs_no_baseline_4panels.pdf", height=5, width=7.5,
  units="in")
# ########################
# ##################################################



# ##################################################
# Make two-panel figure showing predicted vs. actual ADD for the 2 analyses
# (with bones and with swabs) 

# To keep the axes consistent on the following 2 plots, we figure out the
# largest value that will need to be plotted on either plot.
maxaxis <- max(c(wBonesRF$predicted, wBonesRF$y, wSwabsRF$predicted,
  wSwabsRF$y))

# ########################
# Analysis 1 (with bones): Predicted vs. actual ADD

# Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=wBonesRF$predicted,
  actual=wBonesRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
# Force Rsq to be printed to 2 decmial places.
Rsq <- with(predvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
wBonesScatterPanel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=5700, hjust=0, label=paste("R^2  ==", deparse(Rsq)),
    parse=T) +
  annotate("text", x=50, y=5300, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, maxaxis), y=c(0, maxaxis)) +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  labs(x="Actual Accumulated Degree Days",
    y="Predicted Accumulated Degree Days")

wBonesScatterPanel <- annotate_figure(wBonesScatterPanel,
  top=text_grob("Using bones, without baseline obs", face="bold", size=14, vjust=1))
# ########################


# ########################
# Analysis 2 (with swabs): Predicted vs. actual ADD

# Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=wSwabsRF$predicted, actual=wSwabsRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
# Force Rsq to be printed to 2 decmial places.
Rsq <- with(predvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
wSwabsScatterPanel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=5700, hjust=0, label=paste("R^2  ==", deparse(Rsq)),
    parse=T) +
  annotate("text", x=50, y=5300, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, maxaxis), y=c(0, maxaxis)) +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  labs(x="Actual Accumulated Degree Days", y="")

wSwabsScatterPanel <- annotate_figure(wSwabsScatterPanel,
  top=text_grob("Using swabs, without baseline obs", face="bold", size=14, vjust=1))
# ########################


# ########################
plot_grid(wBonesScatterPanel, wSwabsScatterPanel, nrow=1)

ggsave(file="rr_scapula_bones_vs_swabs_no_baseline_predicted_vs_actual_ADD.pdf",
  height=4, width=7.5, units="in")
# ########################
# ##################################################
