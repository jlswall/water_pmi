library("tidyverse")
library("figdim")
library("randomForest")
library("pals")  # For color palettes.
library("cowplot")
library("ggpubr")


#  THIS CODE HAS BEEN ALTERED SO THAT IT ONLY CONSIDERS INFLUENTIAL TAXA BASED
#  ON MODELS WHICH USE THE BASELINE OBSERVATIONS (ADD 0).



# ##################################################
# Get top influential taxa (based on %IncMSE) for all models which used bones.


# ########################
# Make a list of all the files containing all the models of interest.

# ##########
# For analyses using BONES

# For analyses using BONES and RIBS:
#   For ribs, Henley Lake and Rice Rivers, WITH baseline observations.
modelFiles <- c(
  "HenleyLake/w_bones/bacteria/use_families/w_ribs/w_baseline/families_ribs_rfmodel.RData",
  "RiceRivers/w_bones/bacteria/use_families/w_ribs/w_baseline/families_ribs_rfmodel.RData"
)
# #   For ribs, Henley Lake and Rice Rivers, WITHOUT baseline observations.
# modelFiles <- c(
#   modelFiles,
#   "HenleyLake/w_bones/bacteria/use_families/w_ribs/no_baseline/families_ribs_rfmodel.RData",
#   "RiceRivers/w_bones/bacteria/use_families/w_ribs/no_baseline/families_ribs_rfmodel.RData"
# )


# For analyses using BONES and SCAPULAE:
#   For scapulae, Henley Lake and Rice Rivers, WITH baseline observations.
modelFiles <- c(
  modelFiles,
  "HenleyLake/w_bones/bacteria/use_families/w_scapulae/w_baseline/families_scapulae_rfmodel.RData",
  "RiceRivers/w_bones/bacteria/use_families/w_scapulae/w_baseline/families_scapulae_rfmodel.RData"
)
# #   For scapulae, Henley Lake and Rice Rivers, WITHOUT baseline observations.
# modelFiles <- c(
#   modelFiles,
#   "HenleyLake/w_bones/bacteria/use_families/w_scapulae/no_baseline/families_scapulae_rfmodel.RData",
#   "RiceRivers/w_bones/bacteria/use_families/w_scapulae/no_baseline/families_scapulae_rfmodel.RData"
# )


# For analyses using BONES and COMBINED ribs and scapulae:
#   For combined r&s, Henley Lake and Rice Rivers, WITH baseline observations.
modelFiles <- c(
  modelFiles,
  "HenleyLake/w_bones/bacteria/use_families/both_ribs_scapulae/w_baseline/families_combined_rfmodel.RData",
  "RiceRivers/w_bones/bacteria/use_families/both_ribs_scapulae/w_baseline/families_combined_rfmodel.RData"
)
# #   For combined r&s, Henley Lake and Rice Rivers, WITHOUT baseline observations.
# modelFiles <- c(
#   modelFiles,
#   "HenleyLake/w_bones/bacteria/use_families/both_ribs_scapulae/no_baseline/families_combined_rfmodel.RData",
#   "RiceRivers/w_bones/bacteria/use_families/both_ribs_scapulae/no_baseline/families_combined_rfmodel.RData"
# )
# ##########



# ##########
# For analyses using SWABS

# For analyses using SWABS and RIBS:
#   For ribs, Henley Lake and Rice Rivers, WITH baseline observations.
modelFiles <- c(
  modelFiles,
  "HenleyLake/w_swabs/bacteria/use_families/w_ribs/w_baseline/families_ribs_rfmodel.RData",
  "RiceRivers/w_swabs/bacteria/use_families/w_ribs/w_baseline/families_ribs_rfmodel.RData"
)
# #   For ribs, Henley Lake and Rice Rivers, WITHOUT baseline observations.
# modelFiles <- c(
#   modelFiles,
#   "HenleyLake/w_swabs/bacteria/use_families/w_ribs/no_baseline/families_ribs_rfmodel.RData",
#   "RiceRivers/w_swabs/bacteria/use_families/w_ribs/no_baseline/families_ribs_rfmodel.RData"
# )


# For analyses using SWABS and SCAPULAE:
#   For scapulae, Henley Lake and Rice Rivers, WITH baseline observations.
modelFiles <- c(
  modelFiles,
  "HenleyLake/w_swabs/bacteria/use_families/w_scapulae/w_baseline/families_scapulae_rfmodel.RData",
  "RiceRivers/w_swabs/bacteria/use_families/w_scapulae/w_baseline/families_scapulae_rfmodel.RData"
)
# #   For scapulae, Henley Lake and Rice Rivers, WITHOUT baseline observations.
# modelFiles <- c(
#   modelFiles,
#   "HenleyLake/w_swabs/bacteria/use_families/w_scapulae/no_baseline/families_scapulae_rfmodel.RData",
#   "RiceRivers/w_swabs/bacteria/use_families/w_scapulae/no_baseline/families_scapulae_rfmodel.RData"
# )


# For analyses using SWABS and COMBINED ribs and scapulae:
#   For combined r&s, Henley Lake and Rice Rivers, WITH baseline observations.
modelFiles <- c(
  modelFiles,
  "HenleyLake/w_swabs/bacteria/use_families/both_ribs_scapulae/w_baseline/families_combined_rfmodel.RData",
  "RiceRivers/w_swabs/bacteria/use_families/both_ribs_scapulae/w_baseline/families_combined_rfmodel.RData"
)
# #   For combined r&s, Henley Lake and Rice Rivers, WITHOUT baseline observations.
# modelFiles <- c(
#   modelFiles,
#   "HenleyLake/w_swabs/bacteria/use_families/both_ribs_scapulae/no_baseline/families_combined_rfmodel.RData",
#   "RiceRivers/w_swabs/bacteria/use_families/both_ribs_scapulae/no_baseline/families_combined_rfmodel.RData"
# )
# ##########
# ########################



# ########################
# For each model, read in top n most influential taxa, based on "%IncMSE".

# Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 6

# We store the list of taxa name as we go through the loop in this object:
infltaxa <- NULL

for (iFile in modelFiles){
  
  # Load file containing each fitted model in turn.
  load(iFile)

  # The fitted model is named "rf" in each file.
  toptaxa <- importance(rf) %>%
    as.data.frame() %>% 
    rownames_to_column("family") %>%
    as_tibble() %>%
    top_n(n, wt=`%IncMSE`) %>%
    pull(family)

  # Append these taxa to the list for all the models which came before.
  infltaxa <- c(infltaxa, toptaxa)
}

# Many of the influential taxa are influential in more than one analysis, so we
# remove the duplicates.  Also, all the family-level taxa have a prefix of "f__"
# attached to the family name.  We remove this prefix.
infltaxa <- str_remove(unique(infltaxa), "f__")

rm(rf, toptaxa, iFile)
# ##########
# ##################################################




# ##################################################
# Set up colors for each of the influential taxa.

# fullSetColors <- alphabet2()
# Now remove some colors which could be hard to distinguish.
# excludeColors <- c("amethyst", "honey", "iron", "jade", "pink", "quagmire", "sea")
# fullSetColors <- fullSetColors[-which(names(fullSetColors) %in% excludeColors)]
# taxaColors <- fullSetColors[1:length(infltaxa)]
# names(taxaColors) <- infltaxa

# Get n=length(infltaxa) colors from the "polychrome" palette in pals package.
taxaColors <- polychrome(n=length(infltaxa))
names(taxaColors) <- infltaxa


# Write these out in case we want to use them later in other plots to ensure
# consistent colors.
write.csv(as.data.frame(taxaColors) %>% rownames_to_column("taxon"),
  file="color_table_infl_taxa_ribs_scap.csv", row.names=F)
# ##################################################
