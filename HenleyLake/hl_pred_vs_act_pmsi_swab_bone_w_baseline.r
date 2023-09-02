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



# ########################
# Swabs, ribs

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(swabRibT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(swabRibT$resids^2)), 2)
swabRibScatter <- ggplot(swabRibT, aes(x=actual, y=predicted)) +
  geom_point(aes(shape=type, color=type), size=1.75, show.legend=F) +
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
  theme(plot.margin=unit(c(t=1, r=1, b=1, l=1), "pt"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold")) +
  labs(x="Actual Accumulated Degree Days",
       y="Predicted Accumulated Degree Days")
       # title="Swabs from ribs")
# ########################      


# ########################
# Swabs, scapulae

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(swabScapT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(swabScapT$resids^2)), 2)
swabScapScatter <- ggplot(swabScapT, aes(x=actual, y=predicted)) +
  geom_point(aes(shape=type, color=type), size=1.75, show.legend=F) +
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
  theme(plot.margin=unit(c(t=1, r=1, b=1, l=1), "pt"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold")) +
  labs(x="Actual Accumulated Degree Days",
       y="Predicted Accumulated Degree Days")
       # title="Swabs from scapulae")
# ########################      


# ########################
# Swabs, ribs and scapulae (combined)

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(swabCombT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(swabCombT$resids^2)), 2)
swabCombScatter <- ggplot(swabCombT, aes(x=actual, y=predicted)) +
  geom_point(aes(shape=type, color=type), size=1.75, show.legend=F) +
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
  theme(plot.margin=unit(c(t=1, r=1, b=1, l=1), "pt"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold")) +
  labs(x="Actual Accumulated Degree Days",
       y="Predicted Accumulated Degree Days")
       # title="Swabs from ribs and scapulae")
# ########################      



# ########################
# Bones, ribs

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(boneRibT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(boneRibT$resids^2)), 2)
boneRibScatter <- ggplot(boneRibT, aes(x=actual, y=predicted)) +
  geom_point(aes(shape=type, color=type), size=1.75, show.legend=F) +
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
  theme(plot.margin=unit(c(t=1, r=1, b=1, l=1), "pt"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold")) +
  labs(x="Actual Accumulated Degree Days",
       y="Predicted Accumulated Degree Days")
  #      title="Ribs bones")
# ########################      


# ########################
# Bones, scapulae

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(boneScapT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(boneScapT$resids^2)), 2)
boneScapScatter <- ggplot(boneScapT, aes(x=actual, y=predicted)) +
  geom_point(aes(shape=type, color=type), size=1.75, show.legend=F) +
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
  theme(plot.margin=unit(c(t=1, r=1, b=1, l=1), "pt"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold")) +
  labs(x="Actual Accumulated Degree Days",
       y="Predicted Accumulated Degree Days")
       # title="Scapulae bones")
# ########################      


# ########################
# Bones, ribs and scapulae (combined)

# Force Rsq to be printed to 2 decmial places.
Rsq <- with(boneCombT, format(round(cor(actual, predicted)^2, 2), nsmall=2))
# RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(boneCombT$resids^2)), 2)
boneCombScatter <- ggplot(boneCombT, aes(x=actual, y=predicted)) +
  geom_point(aes(shape=type, color=type), size=1.75, show.legend=F) +
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
  theme(plot.margin=unit(c(t=1, r=1, b=1, l=1), "pt"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold")) +
  labs(x="Actual Accumulated Degree Days",
       y="Predicted Accumulated Degree Days")
       # title="Ribs and scapulae, bones")
# ########################      



# ########################
# To see default plot margin sizes: theme_get()$plot.margin
swabTitle <- ggdraw() + draw_label("Swabs", fontface="bold")
boneTitle <- ggdraw() + draw_label("Bones", fontface="bold")

ribTitle <- ggdraw() + draw_label("Ribs", fontface="bold", angle=90)
scapTitle <- ggdraw() + draw_label("Scapulae", fontface="bold", angle=90)
combTitle <- ggdraw() + draw_label("Ribs & scapulae", fontface="bold", angle=90)

topTitleRow <- plot_grid("", swabTitle, boneTitle,
                         nrow=1, ncol=3, rel_widths=c(0.1, 1, 1))
ribRow <- plot_grid(ribTitle, swabRibScatter, boneRibScatter,
                     nrow=1, ncol=3, rel_widths=c(0.1, 1, 1))
scapRow <- plot_grid(scapTitle, swabScapScatter, boneScapScatter,
                     nrow=1, ncol=3, rel_widths=c(0.1, 1, 1))
combRow <- plot_grid(combTitle, swabCombScatter, boneCombScatter,
                     nrow=1, ncol=3, rel_widths=c(0.1, 1, 1))

plot_grid(topTitleRow, ribRow, scapRow,combRow,
          nrow=4, ncol=1, rel_heights=c(0.15, 1, 1, 1))


ggsave(file="hl_scatter_pred_vs_act_pmsi_swab_bone_w_baseline.pdf", height=9,
       width=6.5, units="in")
# ########################
# ##################################################
