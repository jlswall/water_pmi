library("tidyverse")
library("figdim")
library("randomForest")
library("scales")  # For hue_pal() function.
library("cowplot")
library("ggpubr")



# ##################################################
# Read in the final fitted models.

# The objects containing the fitted models were names "rf" when they were
# originally run (one at a time), but we need to give them re-name them here so
# they don't get mixed up.

# Models using SWABS: ribs, scapulae, and ribs/scapulae combined.
#   Ribs
load("w_swabs/bacteria/use_families/w_ribs/w_baseline/families_ribs_rfmodel.RData")
swabRibRF <- rf
rm(rf)
#   Scapulae
load("w_swabs/bacteria/use_families/w_scapulae/w_baseline/families_scapulae_rfmodel.RData")
swabScapRF <- rf
rm(rf)
#   Combined (ribs and scapulae)
load("w_swabs/bacteria/use_families/both_ribs_scapulae/w_baseline/families_combined_rfmodel.RData")
swabCombRF <- rf
rm(rf)


# Models using BONES: ribs, scapulae, and ribs/scapulae combined.
#   Ribs
load("w_bones/bacteria/use_families/w_ribs/w_baseline/families_ribs_rfmodel.RData")
boneRibRF <- rf
rm(rf)
#   Scapulae
load("w_bones/bacteria/use_families/w_scapulae/w_baseline/families_scapulae_rfmodel.RData")
boneScapRF <- rf
rm(rf)
#   Combined (ribs and scapulae)
load("w_bones/bacteria/use_families/both_ribs_scapulae/w_baseline/families_combined_rfmodel.RData")
boneCombRF <- rf
rm(rf)
# ##################################################



# ##################################################
# Build tibble for each of the models, containing predicted and actual columns.

### REFER TO FILE HenleyLake/w_swabs/bacteria/use_families/combined_w_baseline.r

# Models using SWABS:
swabRibT <- with(swabRibRF,
                 as_tibble(data.frame(type="Rib", actual=y, predicted=predicted,
                                      resids=y-predicted)
                           )
                 )
swabScapT <- with(swabScapRF,
                  as_tibble(data.frame(type="Scapula", actual=y,
                                       predicted=predicted, resids=y-predicted)
                            )
                  )
# For the combined model, we read in actual and predicted values for which
# the type (rib or scapula) was saved for each sample.
swabCombT <- read_csv("w_swabs/bacteria/use_families/both_ribs_scapulae/w_baseline/predicted_actual_w_type.csv")
swabCombT <- swabCombT %>% mutate(resids = actual - predicted)


# Models using BONES:
boneRibT <- with(boneRibRF,
                 as_tibble(data.frame(type="Rib", actual=y, predicted=predicted,
                                      resids=y-predicted)
                 )
)
boneScapT <- with(boneScapRF,
                  as_tibble(data.frame(type="Scapula", actual=y,
                                       predicted=predicted, resids=y-predicted)
                  )
)
# For the combined model, we read in actual and predicted values for which
# the type (rib or scapula) was saved for each sample.
boneCombT <- read_csv("w_bones/bacteria/use_families/both_ribs_scapulae/w_baseline/predicted_actual_w_type.csv")
boneCombT <- boneCombT %>% mutate(resids = actual - predicted)
# ##################################################



# ##################################################
# Make six-panel figure showing predicted vs. actual ADD.  Swabs in left column,
# bones in right column.  Rows are ribs, scapulae, and ribs/scapulae combined.

# To keep the axes consistent on all plots, we figure out the largest value that
# will need to be plotted on the x-axis or y-axis of any of the plots.
maxAxis <- max(swabRibT %>% select("actual", "predicted"),
               swabScapT %>% select("actual", "predicted"),
               swabCombT %>% select("actual", "predicted"),
               boneRibT %>% select("actual", "predicted"),
               boneScapT %>% select("actual", "predicted"),
               boneCombT %>% select("actual", "predicted")
               )


# Assign different colors and symbols to ribs and scapula and use these for each
# scatterplot.
typeColors <- hue_pal()(2)  # Based on way ggplot assigns colors for categories
names(typeColors) <- c("Rib", "Scapula")
typeSymbols <- c("circle", "square")
names(typeSymbols) <- c("Rib", "Scapula")


### WORKING HERE - NEED TO MAKE SYMBOLS BIGGER.


# ########################
# Swabs, ribs

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(swabRibT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(swabRibT$resids^2)), 2)
ribscatterPanel <- ggplot(swabRibT, aes(x=actual, y=predicted)) +
  geom_point(aes(shape=type, color=type), show.legend=F) +
  scale_color_manual(values=typeColors) +
  scale_shape_manual(values=typeSymbols) +
  geom_abline(slope=1, intercept=0) +
  lims(x=c(0, maxAxis), y=c(0, maxAxis)) +
  annotate("text", x=0.02*maxAxis, y=0.95*maxAxis, hjust=0,
           label=paste("R^2  ==", deparse(Rsq)), parse=T) +
  annotate("text", x=0.02*maxAxis, y=0.88*maxAxis, hjust=0,
           label=paste("RMSE = ", RMSE)) +
  coord_fixed(ratio=1) +
  theme_bw() +
  theme(axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold")) +
  labs(x="Actual Accumulated Degree Days",
       y="Predicted Accumulated Degree Days",
       title="Swabs, ribs")
# ########################      



# ########################
# Swabs from ribs: Predicted vs. actual ADD

# Make a tibble of actual and predicted values for each observation.
predvactT <- as_tibble(data.frame(predicted=swabsRibsRF$predicted,
  actual=swabsRibsRF$y))
predvactT$resids <- with(predvactT, actual - predicted)
# Force Rsq to be printed to 2 decmial places.
Rsq <- with(predvactT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(predvactT$resids^2)), 2)  
swabsRibsScatter <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=5000, hjust=0, label=paste("R^2  ==", deparse(Rsq)),
    parse=T) +
  annotate("text", x=50, y=4600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, maxaxis), y=c(0, maxaxis)) +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  labs(x="Actual Accumulated Degree Days",
    y="Predicted Accumulated Degree Days")

wBonesScatterPanel <- annotate_figure(wBonesScatterPanel,
  top=text_grob("Using bones, with baseline obs", face="bold", size=14, vjust=1))
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
  annotate("text", x=50, y=5000, hjust=0, label=paste("R^2  ==", deparse(Rsq)),
    parse=T) +
  annotate("text", x=50, y=4600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, maxaxis), y=c(0, maxaxis)) +
  theme(axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10)) +
  labs(x="Actual Accumulated Degree Days", y="")

wSwabsScatterPanel <- annotate_figure(wSwabsScatterPanel,
  top=text_grob("Using swabs, with baseline obs", face="bold", size=14, vjust=1))
# ########################


# ########################
plot_grid(wBonesScatterPanel, wSwabsScatterPanel, nrow=1)

ggsave(file="hl_scapula_bones_vs_swabs_with_baseline_predicted_vs_actual_ADD.pdf",
  height=4, width=7.5, units="in")
# ########################
# ##################################################
