library("tidyverse")
library("randomForest")
library("cowplot")
library("ggpubr")


# Are we dealing with phlya, orders, or families?
taxalevel <- "families"

# Are we dealing with Henley Lake or Rice Rivers data?
location <- "RiceRivers"

# Get the top "n" (whether 8, 10, whatever) influential taxa.
# We've set up enough unique colors to for the top 6 infl. taxa from each model.
nInflTaxa <- 6


# ##################################################
# Read in dataset with bones' family-level taxa.

# Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("../", location, "/w_bones/bacteria/", taxalevel,
                        "_massaged.csv"))

# All the family-level taxa have a prefix of "f__" attached to the
# family name.  We want to remove this prefix.  The only taxon that
# is not affected is the "Rare" category, because it's the only one
# which doesn't have the prefix.
allT$taxon <- str_remove(allT$taxon, "f__")

# Filter data to ribs and scapulae (no water, etc.).
bonesT <- allT %>% filter(type=="Rib" | type=="Scapula")

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

# Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("../", location, "/w_swabs/bacteria/", taxalevel,
                        "_massaged.csv"))

# All the family-level taxa have a prefix of "f__" attached to the
# family name.  We want to remove this prefix.  The only taxon that
# is not affected is the "Rare" category, because it's the only one
# which doesn't have the prefix.
allT$taxon <- str_remove(allT$taxon, "f__")

# Filter data to ribs and scapulae (no water, etc.).
swabsT <- allT %>% filter(type=="Rib" | type=="Scapula")

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

load(paste0("../", location, "/w_bones/bacteria/use_families/both_ribs_scapulae/w_baseline/families_combined_rfmodel.RData"))
# The object was named "rf", but we don't want to mix up the 2 models being
# compared.  So we re-name it.
wBonesRF <- rf
rm(rf)

load(paste0("../", location, "/w_swabs/bacteria/use_families/both_ribs_scapulae/w_baseline/families_combined_rfmodel.RData"))
# As before, we give the saved object a more specific name.
wSwabsRF <- rf
rm(rf)
# ##################################################



# ##################################################
# Find influential taxa for both analyses, and set colors for these to be used
# consistently across plot panels.


# ########################
# Get top n influential taxa for first analysis.

# Turn importance measures into a tibble, sorted by IncNodePurity in
# increasing order.
wBonesImportT <- importance(wBonesRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  top_n(nInflTaxa, wt=`%IncMSE`) %>%
  arrange(`%IncMSE`)
# Remove the "f__" from the family taxon names.
wBonesImportT$family <- str_remove(wBonesImportT$family, "f__")
# ########################


# ########################
# Get top n influential taxa for second analysis.

# Turn importance measures into a tibble, sorted by IncNodePurity in
# increasing order.
wSwabsImportT <- importance(wSwabsRF) %>%
  as.data.frame() %>% 
  rownames_to_column("family") %>%
  as_tibble() %>%
  top_n(nInflTaxa, wt=`%IncMSE`) %>%
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
colorpaletteDF <- read.csv("color_table_infl_taxa_ribs_scap.csv")
assignedColors <- colorpaletteDF$taxaColors
names(assignedColors) <- colorpaletteDF$"taxon"
rm(colorpaletteDF)
# ########################



# ########################
# Bar plots of influential taxa (left side of panel of figures).


# #####
# Set up range for the bars in the bar chart so that axes can be
# consistent between bar charts for both analyses.

barMax <- ceiling(max(c(wBonesImportT %>% pull(`%IncMSE`),
                        wSwabsImportT %>% pull(`%IncMSE`))))
# #####


# #####
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
  scale_fill_manual(values=assignedColors)
# #####


# #####
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
  scale_fill_manual(values=assignedColors)
# #####
# ########################



# ########################
# Analyses 1 and 2: Show line plots of changing relative abundance for the top taxa.

# ###########
# Get numbers for analysis 1 (with bones):

# Save the names of the families that are in the top nInflTaxa in terms of
# %IncMSE.
topChoices <- as.character(wBonesImportT %>% top_n(nInflTaxa, wt=`%IncMSE`) %>%
                             pull(family))

# Find the percentages for these taxa.
chooseT <- bonesT %>%
  filter(taxon %in% topChoices)

# Average the value across samples for each taxa and each day.
summBonesT <- chooseT %>%
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample),
            medianPercByDay=100*median(fracBySample))

rm(topChoices, chooseT)
# ###########


# ###########
# Get numbers for analysis 2 (with swabs):

# Save the names of the families that are in the top nInflTaxa in
# terms of %IncMSE.
topChoices <- as.character(wSwabsImportT %>% top_n(nInflTaxa, wt=`%IncMSE`) %>%
                             pull(family))

# Find the percentages for these taxa.
chooseT <- swabsT %>%
  filter(taxon %in% topChoices)

# Average the value across samples for each taxa and each day.
summSwabsT <- chooseT %>%
  group_by(taxon, degdays) %>%
  summarize(meanPercByDay=100*mean(fracBySample),
            medianPercByDay=100*median(fracBySample))

rm(topChoices, chooseT)
# ###########


# ###########
# Find the maximum value for the x- and y-axes for setting consistent axes.

xMax <- max(c(summBonesT$degdays, summSwabsT$degdays))
yMax <- max(c(summBonesT$meanPercByDay, summSwabsT$meanPercByDay))
# ###########


# ###########
# Draw plot of average relative abundance vs. time for these influential taxa.
wBonesLinePanel <- ggplot(summBonesT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(linewidth=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_x_continuous(limits=c(0, xMax)) +
  scale_y_continuous(limits=c(0, yMax), expand=c(0,0)) +
  theme_minimal() +
  #  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  labs(x="", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
  scale_color_manual(values=assignedColors)


# #####
# Draw plot of average relative abundance vs. time for these five
# influential taxa.
# dev.new(width=4.5, height=4)
wSwabsLinePanel <- ggplot(summSwabsT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(linewidth=1.25, aes(color=taxon), show.legend=FALSE) +
  scale_x_continuous(limits=c(0, xMax)) +
  scale_y_continuous(limits=c(0, yMax), expand=c(0,0)) +
  theme_minimal() +
  labs(x="Accumulated Degree Days", y="Relative Abundance") +
  theme(axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
  scale_color_manual(values=assignedColors)
# ########################




# ########################
# Build 4-panel figure with a title in its own "row" at the top.

# For the title.
topTitle <- ggdraw() + draw_label("Ribs & scapulae", fontface="bold")

# First row of panel shows info about taxa in second analysis (bar chart,
# relative abundance graphs).
plotrow1 <- plot_grid(wBonesbarPanel, wBonesLinePanel, rel_widths=c(1, 2),
  nrow=1)
plotrow1 <- annotate_figure(plotrow1,
  left=text_grob("Bones", face="bold", rot=90, size=10))

# Second row shows info about taxa in second analysis (bar chart, relative
# abundance graphs).
plotrow2 <- plot_grid(wSwabsbarPanel, wSwabsLinePanel, rel_widths=c(1, 2),
  nrow=1)
plotrow2 <- annotate_figure(plotrow2,
  left=text_grob("Swabs", face="bold", rot=90, size=10))

# Put the rows together to make 4-panel figure.
plot_grid(topTitle, plotrow1, plotrow2, nrow=3, rel_heights=c(0.1, 1, 1))

ggsave(file="rr_indic_taxa_combined_rib_scap_bones_vs_swabs_with_baseline.pdf", height=5,
       width=7.5, units="in")
# ########################
# ##################################################
